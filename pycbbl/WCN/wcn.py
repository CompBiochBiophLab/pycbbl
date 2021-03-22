from Bio.PDB import *
from collections import OrderedDict
import numpy as np
import re

class WCNObject:
    """
    Class that generates different types of WCN profiles for defined structures.

    Attributes
    ----------
    input_file : str
        Path to the input file structure.
    structure : Bio.PDB.Structure
        Biopython structure object.
    all_atoms : list
        List containing all the atom objects in the structure/
    calpha_atoms : list
        List containing all the atom objects corresponding to the alpha carbon atoms
        of each residue in the protein structure.
    residues : list
        List containing all the residues in the protein structure.
    chains : list
        List containing all the chains in the protein structure.
    bf_aa : numpy.ndarray
        B-factor profile
    wcn_ca : numpy.ndarray
        Alpha carbon atoms WCN profile.
    wcn_aa : numpy.ndarray
        All atoms WCN profile.
    wcn_pr : numpy.ndarray
        All atom WCN profile averaged by residue.

    Methods
    -------
    removeWaters()
        Exclude all water atoms from the WCN analysis.
    removeLigands()
        Exclude all ligand atoms from the WCN analysis.
    removeHydrogens()
        Exclude all hydrogen atoms from the WCN analysis.
    getBfactors()
        Reads the b-factor data of the input structure file.
    squareInverseDistanceMatrix(atoms)
        Calculate the inverse distance matrix for the set of atoms given.
    calculateCAlpha()
        Calculates an alpha carbon WCN profile.
    calculateAllAtom()
        Calculates an all atom WCN profile.
    calculateAveragePerResidue()
        Calculates an all atom WCN profile averaged by residue.
    """

    def __init__(self, input_file):
        """
        Initializes an instance of the WCNObject class from a PDB or CIF files containing
        a protein structure.

        Parameters
        ----------
        input_file : str
            Path to the structural file.
        """
        if input_file.endswith('.pdb'):
            parser = PDBParser()
        elif input_file.endswith('.cif'):
            parser = MMCIFParser()

        self.input_file = input_file
        self.structure = parser.get_structure(input_file, input_file)
        self.all_atoms = [a for a in self.structure.get_atoms()]
        self.calpha_atoms = [a for a in self.structure.get_atoms() if a.name == 'CA']
        self.residues = [r for r in self.structure.get_residues()]
        self.chains = [c for c in self.structure.get_chains()]

        #Check that all defined residues have a C-alpha atom
        for r in self.residues:
            if len([a for a in r if a.name == 'CA']) != 1:
                if 'H_' != r.get_full_id()[3][0][:2] and  r.get_full_id()[3][0] != 'W':
                    print('Warning! defined residue '+str(r.get_resname())+' '+str(r.get_full_id()[-1][1])+' of chain '+r.get_parent().id+' do not have a C-alpha atom')

        self.wcn_aa = None
        self.wcn_ca = None
        self.wcn_pr = None
        self.bf_aa = None

    def removeWaters(self):
        """
        Exclude all water atoms from the analysis. It removes all the water atoms
        from the set of all atoms contained in the attribute 'all_atoms'.
        """
        atomsToRemove = set()
        residuesToRemove = set()

        for a in self.all_atoms:
            if a.get_parent().get_full_id()[3][0] == 'W':
                atomsToRemove.add(a)
                residuesToRemove.add(a.get_parent())

        self.all_atoms = [a for a in self.all_atoms if a not in atomsToRemove]
        self.residues = [r for r in self.residues if r not in residuesToRemove]

    def removeLigands(self):
        """
        Exclude all hetero-atoms from the analysis. It removes all the atoms marked
        as heteroatoms ('HETATM') from the set of all atoms contained in the attribute
        'all_atoms'.
        """
        atomsToRemove = set()
        residuesToRemove = set()

        for a in self.all_atoms:
            if 'H_' in a.get_parent().get_full_id()[3][0]:
                atomsToRemove.add(a)
                residuesToRemove.add(a.get_parent())

        self.all_atoms = [a for a in self.all_atoms if a not in atomsToRemove]
        self.calpha_atoms = [a for a in self.all_atoms if a not in atomsToRemove and a.name == 'CA']
        self.residues = [r for r in self.residues if r not in residuesToRemove]

    def removeHydrogens(self):
        """
        Exclude all hydrogen atoms from the analysis. It removes all the hydrogen
        atoms from the set of all atoms contained in the attribute 'all_atoms'.
        """

        _hydrogen = re.compile("[123 ]*H.*")
        atomsToRemove = set()
        for a in self.all_atoms:
            if _hydrogen.match(a.name):
                atomsToRemove.add(a)
        self.all_atoms = [a for a in self.all_atoms if a not in atomsToRemove]

    def getBfactors(self, normalized=True):
        """
        Get B-factor column of values.

        Parameters
        ----------
        normalized : bool
            Normalize B-factor column profile?

        Returns
        -------
        bf_aa : np.ndarray
            B-factor column profile for all atoms.
        """

        if not isinstance(self.bf_aa, np.ndarray):

            self.bf_aa = np.array([atom.bfactor for atom in self.all_atoms])

            if normalized == True:
                self.bf_aa = _normalize(self.bf_aa)

        return self.bf_aa

    def squareInverseDistanceMatrix(self, atoms, CUDA=False):
        """
        Calculate the inverse distance matrix for all the atoms in the given set.

        Parameters
        ----------
        atoms : list
            List of atoms objects.
        CUDA : bool
            Whether to use CUDA to calculate the matrix.

        Returns
        -------
        distance_matrix : np.ndarray
            Squared inverse atom (pair-wise) distances array
        """

        # Cuda implementation
        # if CUDA == True:
        #     coords = [a.coord for a in atoms]
        #     size = len(coords)
        #     distance = 1/cp.sum((cp.array(coords*len(coords)).reshape(size,size,3)-cp.array(coords).reshape(size,1,3))**2,axis=2)
        #     cp.fill_diagonal(distance,0)
        #     distance_matrix = cp.asnumpy(distance)
        #
        #     return distance_matrix

        # else:
        distance_matrix = np.zeros((len(atoms), len(atoms)))
        for i in range(len(atoms)):
            for j in range(i+1,len(atoms)):
                distance_matrix[i][j] = np.linalg.norm(atoms[i].coord-atoms[j].coord)
                distance_matrix[i][j] = np.square(distance_matrix[i][j])
                distance_matrix[i][j] = np.reciprocal(distance_matrix[i][j])
                distance_matrix[j][i] = distance_matrix[i][j]

        return distance_matrix

    def calculateCAlpha(self, normalized=True, CUDA=False):
        """
        Calculates the WCN profile only for alpha carbon atoms.

        Parameters
        ----------
        normalized : bool
            Normalize WCN profile?
        CUDA : bool
            Whether to use CUDA to calculate the matrix (requires cupy).

        Returns
        -------
        wcn_ca : np.ndarray
            WCN profile for alpha carbon atoms
        """
        if CUDA == True:
             self.wcn_ca = self.squareInverseDistanceMatrix(self.calpha_atoms, CUDA=True)
        else:
            self.wcn_ca = self.squareInverseDistanceMatrix(self.calpha_atoms)

        self.wcn_ca = np.reciprocal(np.sum(self.wcn_ca, axis=0))

        if normalized == True:
            self.wcn_ca = _normalize(self.wcn_ca)

        return self.wcn_ca

    def calculateAllAtom(self, normalized=True, CUDA=False):
        """
        Calculates the WCN profile for all atoms.

        Parameters
        ----------
        normalized : bool
            Normalize WCN profile?
        CUDA : bool
            Whether to use CUDA to calculate the matrix (requires cupy).

        Returns
        -------
        wcn_aa : np.ndarray
            WCN profile for all atoms
        """

        if CUDA == True:
             self.wcn_aa = self.squareInverseDistanceMatrix(self.all_atoms, CUDA=True)
        else:
            self.wcn_aa = self.squareInverseDistanceMatrix(self.all_atoms)

        self.wcn_aa = np.reciprocal(np.sum(self.wcn_aa, axis=0))

        if normalized == True:
            self.wcn_aa = _normalize(self.wcn_aa)

        return self.wcn_aa

    def calculateAveragePerResidue(self, normalized=True, CUDA=False):
        """
        Calculates first an all atom WCN profile and then averages it by residue.

        Parameters
        ----------
        normalized : bool
            Normalize WCN profile?
        CUDA : bool
            Whether to use CUDA to calculate the matrix (requires cupy).

        Returns
        -------
        wcn_pr : np.ndarray
            All atom WCN profile averaged by residue.
        """

        if self.wcn_aa == None:
            self.calculateAllAtom(normalized=False, CUDA=CUDA)

        self.wcn_pr = np.zeros(len(self.residues))

        for j,r in enumerate(self.residues):
            indexes = []
            for i,a in enumerate(self.all_atoms):
                if a.get_parent() == r:
                    indexes.append(i)
            self.wcn_pr[j] = np.average(self.wcn_aa[indexes])

        if normalized == True:
            self.wcn_aa = _normalize(self.wcn_aa)
            self.wcn_pr = _normalize(self.wcn_pr)

        return self.wcn_pr

    def separateByChains(self):
        """
        Separate by chains all the so far calculated WCN profiles.
        """

        self.wcn_aa_chain = {}
        self.wcn_ca_chain = {}
        self.wcn_pr_chain = {}

        for c in self.chains:
            self.wcn_aa_chain[c.id] = None
            self.wcn_ca_chain[c.id] = None
            self.wcn_pr_chain[c.id] = None

        if isinstance(self.wcn_aa, np.ndarray):
            self.wcn_aa_chain = {}
            for i,a in enumerate(self.all_atoms):
                c_id = a.get_parent().get_parent().id
                if c_id not in self.wcn_aa_chain:
                    self.wcn_aa_chain[c_id] = []
                self.wcn_aa_chain[c_id].append(self.wcn_aa[i])
            self.wcn_aa_chain[c_id] = np.array(self.wcn_aa_chain[c_id])

            for c in self.chains:
                self.wcn_aa_chain[c.id] = np.array(self.wcn_ca_chain[c.id])

        if isinstance(self.wcn_ca, np.ndarray):
            self.wcn_ca_chain = {}
            self.residues_in_chain = {}
            for i,a in enumerate(self.calpha_atoms):
                c_id = a.get_parent().get_parent().id
                if c_id not in self.wcn_ca_chain:
                    self.wcn_ca_chain[c_id] = []
                if c_id not in self.residues_in_chain:
                    self.residues_in_chain[c_id] = []
                self.wcn_ca_chain[c_id].append(self.wcn_ca[i])
                self.residues_in_chain[c_id].append(a.get_parent())

            for c in self.residues_in_chain:
                self.wcn_ca_chain[c] = np.array(self.wcn_ca_chain[c])
                try:
                    self.residues_in_chain[c] = np.array(self.residues_in_chain[c])
                except:
                    continue


def _normalize(data):
    """
    Normalize the array of data.

    Parameters
    ----------
    data : np.ndarray
        Array of data to normalize.

    Returns
    -------
    data : np.ndarray
        Normalized array.
    """

    average = np.average(data)
    stdev = np.std(data)
    data = (data - average)/stdev

    return data
