import numpy as np
import re
import os
import shutil
from Bio.PDB import *
from collections import OrderedDict
from zipfile import ZipFile
from ..PDB import saveStructureToPDB, getChainSequence
from ..alignment import mafft
from ..MD import alignment
from scipy.spatial import distance_matrix

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
    chain_ids : list
        List containing all the chain ids in the protein structure.
    bf_aa : numpy.ndarray
        B-factor profile
        wcn_ca : numpy.ndarray
        Alpha carbon atoms WCN profile.
    wcn_aa : numpy.ndarray
        All atoms WCN profile.
    wcn_pr : numpy.ndarray
        All atom WCN profile averaged by residue.
    wcn_per_chain : dict
        All atom WCN profiles for structures represented only by selected subchains.

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
        self.chain_ids = [c.id for c in self.structure.get_chains()]

        #Check that all defined residues have a C-alpha atom
        for r in self.residues:
            if len([a for a in r if a.name == 'CA']) != 1:
                if 'H_' != r.get_full_id()[3][0][:2] and  r.get_full_id()[3][0] != 'W':
                    print('Warning! defined residue '+str(r.get_resname())+' '+str(r.get_full_id()[-1][1])+' of chain '+r.get_parent().id+' do not have a C-alpha atom')

        self.wcn_aa = None
        self.wcn_ca = None
        self.wcn_pr = None
        self.bf_aa = None
        self.bf_ca = None
        self.bf_pr = None
        self.wcn_per_chain = {}

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

        # Update atoms
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
        Get B-factor column of values and create all-atom, C-alpha and per-residue
        profiles. It only returns the all atom b-factor list, however all this can
        be accessed by the attributes: .bf_aa, .bf_ca, .bf_pr, respectively.

        Parameters
        ----------
        normalized : bool
            Normalize B-factor column profile (Z-score)?

        Returns
        -------
        bf_aa : np.ndarray
            B-factor column profile for all atoms.
        """

        if not isinstance(self.bf_aa, np.ndarray):

            # Store all atom
            self.bf_aa = np.array([atom.bfactor for atom in self.all_atoms])

            # Store alpha carbon
            self.bf_ca = np.array([atom.bfactor for atom in self.calpha_atoms])

            # Store average BF per residue
            self.bf_pr = np.zeros(len(self.residues))
            for j,r in enumerate(self.residues):
                indexes = []
                for i,a in enumerate(self.all_atoms):
                    if a.get_parent() == r:
                        indexes.append(i)
                self.bf_pr[j] = np.average(self.bf_aa[indexes])

            if normalized == True:
                self.bf_aa = _normalize(self.bf_aa)
                self.bf_ca = _normalize(self.bf_ca)
                self.bf_pr = _normalize(self.bf_pr)

        return self.bf_aa

    def squareInverseDistanceMatrix(self, atoms):
        """
        Calculate the inverse distance matrix for all the atoms in the given set.

        Parameters
        ----------
        atoms : list
            List of atoms objects.

        Returns
        -------
        distance_matrix : np.ndarray
            Squared inverse atom (pair-wise) distances array
        """

        atoms_xyz = np.array([atoms[i].coord for i in range(len(atoms))])
        dm = distance_matrix(atoms_xyz, atoms_xyz)
        dm = np.square(dm)
        dm = np.reciprocal(dm, out=dm, where=dm!=0)
        np.fill_diagonal(dm, 0)

        return dm

    def calculateCAlpha(self, normalized=True):
        """
        Calculates the WCN profile only for alpha carbon atoms.

        Parameters
        ----------
        normalized : bool
            Normalize WCN profile?

        Returns
        -------
        wcn_ca : np.ndarray
            WCN profile for alpha carbon atoms
        """

        self.wcn_ca = self.squareInverseDistanceMatrix(self.calpha_atoms)
        self.wcn_ca = np.reciprocal(np.sum(self.wcn_ca, axis=0))

        if normalized == True:
            self.wcn_ca = _normalize(self.wcn_ca)

        return self.wcn_ca

    def calculateAllAtom(self, normalized=True):
        """
        Calculates the WCN profile for all atoms.

        Parameters
        ----------
        normalized : bool
            Normalize WCN profile?

        Returns
        -------
        wcn_aa : np.ndarray
            WCN profile for all atoms
        """

        self.wcn_aa = self.squareInverseDistanceMatrix(self.all_atoms)

        self.wcn_aa = np.reciprocal(np.sum(self.wcn_aa, axis=0))

        if normalized == True:
            self.wcn_aa = _normalize(self.wcn_aa)

        return self.wcn_aa

    def calculateAveragePerResidue(self, normalized=True):
        """
        Calculates first an all atom WCN profile and then averages it by residue.

        Parameters
        ----------
        normalized : bool
            Normalize WCN profile?

        Returns
        -------
        wcn_pr : np.ndarray
            All atom WCN profile averaged by residue.
        """

        if isinstance(self.wcn_aa, type(None)):
            self.calculateAllAtom(normalized=normalized)

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

            for c in self.chains:
                self.wcn_aa_chain[c.id] = np.array(self.wcn_aa_chain[c.id])

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
                    self.residues_in_chain[c] = np.array(self.residues_in_chain[c], dtype=object)
                except:
                    continue

        if isinstance(self.wcn_pr, np.ndarray):
            self.wcn_pr_chain = {}
            self.residues_in_chain = {}
            for i,a in enumerate(self.calpha_atoms):
                c_id = a.get_parent().get_parent().id
                if c_id not in self.wcn_pr_chain:
                    self.wcn_pr_chain[c_id] = []
                if c_id not in self.residues_in_chain:
                    self.residues_in_chain[c_id] = []
                self.wcn_pr_chain[c_id].append(self.wcn_pr[i])
                self.residues_in_chain[c_id].append(a.get_parent())

            for c in self.residues_in_chain:
                self.wcn_pr_chain[c] = np.array(self.wcn_pr_chain[c])\

    def calculateAllAtomForChains(self, chain_ids, normalized=True, overwrite=False):

        if isinstance(chain_ids, str):
            chain_ids = [chain_ids]

        if not isinstance(chain_ids, list):
            raise ValueError('Chain IDs must be given as a str or (if more than one) as a list')

        if ''.join(chain_ids) in self.wcn_per_chain and not overwrite:
            print('Already calculated. Give overwrite=True to recalculate')

        for c in chain_ids:
            if c not in self.chain_ids:
                raise ValueError('Chain '+c+' is not present in the protein structure!')

        # Create atoms list
        atoms = []

        # Save atoms
        for r in self.residues:
            if r.get_parent().id in chain_ids:
                for a in r.get_atoms():
                    atoms.append(a)

        # Calculate all atom WCN
        self.wcn_per_chain[''.join(chain_ids)] = self.squareInverseDistanceMatrix(atoms)

        self.wcn_per_chain[''.join(chain_ids)] = np.reciprocal(np.sum(self.wcn_per_chain[''.join(chain_ids)], axis=0))

        if normalized == True:
            self.wcn_per_chain[''.join(chain_ids)] = _normalize(self.wcn_per_chain[''.join(chain_ids)])

    def dumpStructure(self, output_file, wcn_aa_as_bf=False, wcn_ca_as_bf=False,
                      wcn_pr_as_bf=False):
        """
        Write a structure with b-factors or WCN in the b-factors columns. The WCN values
        can be all atom, alpha carbon, or per-residue. In the last two cases the
        same WCN is given to all atoms inside a residue.
        """

        # Create chains objects
        chains = {}
        atoms = []
        residues = []

        # Create output struture objects
        structure = Structure.Structure(0)
        model = Model.Model(0)

        # Save atoms and residues
        for c in self.chains:
            chains[c.id] = Chain.Chain(c.id)
            for r in c.get_residues():
                chains[c.id].add(r)
                residues.append(r)
                for a in r.get_atoms():
                    atoms.append(a)

        if wcn_aa_as_bf:
            assert len(atoms) == len(self.wcn_aa)
            for i,wcn in enumerate(self.wcn_aa):
                atoms[i].bfactor = wcn

        elif wcn_ca_as_bf:
            assert len(residues) == len(self.wcn_ca)
            for i, wcn in enumerate(self.wcn_ca):
                for a in residues[i].get_atoms():
                    a.bfactor = wcn

        elif wcn_pr_as_bf:
            assert len(residues) == len(self.wcn_pr)
            for i, wcn in enumerate(self.wcn_pr):
                for a in residues[i].get_atoms():
                    a.bfactor = wcn

        else:
            assert len(atoms) == len(self.bf_aa)
            for i,bf in enumerate(self.bf_aa):
                atoms[i].bfactor = bf

        for c in chains:
            model.add(chains[c])
        structure.add(model)

        saveStructureToPDB(structure, output_file)

    def saveWCNProfilesToFile(self, file_name):

        directory = '/'.join(file_name.split('/')[:-1])
        file_name = file_name.split('/')[-1]

        os.chdir(directory)

        # Append a gz extension to given file name
        if file_name.endswith('.gz'):
            file_name = file_name.replace('.gz', '')

        # Store WCN into a single compressed file
        created_files = []

        # Store WCN CA profile
        if type(self.wcn_ca).__module__ == np.__name__:
            np.save(file_name+'_wcn_ca.npy', self.wcn_ca)
            created_files.append(file_name+'_wcn_ca.npy')

        # Store WCN AA profile
        if type(self.wcn_aa).__module__ == np.__name__:
            np.save(file_name+'_wcn_aa.npy', self.wcn_aa)
            created_files.append(file_name+'_wcn_aa.npy')

        # Store WCN PR profile
        if type(self.wcn_pr).__module__ == np.__name__:
            np.save(file_name+'_wcn_pr.npy', self.wcn_pr)
            created_files.append(file_name+'_wcn_pr.npy')

        # Store WCN per-chains profiles
        for c in self.wcn_per_chain:
            if type(self.wcn_per_chain[c]).__module__ == np.__name__:
                np.save(file_name+'_wcn_per_chain_'+c+'.npy', self.wcn_per_chain[c])
                created_files.append(file_name+'_wcn_per_chain_'+c+'.npy')

        # Create zip object
        zipObj = ZipFile(file_name+'.gz', 'w')

        # Add files to zip and then remove them
        for f in created_files:
            zipObj.write(f)
            os.remove(f)

        # Close Zip object
        zipObj.close()

        os.chdir('..'*len(directory.split()))


    def readWCNProfilesFromFile(self, wcn_file, verbose=True):

        directory = '/'.join(wcn_file.split('/')[:-1])
        wcn_file = wcn_file.split('/')[-1]

        # Append .gz to the file name
        if wcn_file.endswith('.gz'):
            wcn_file = wcn_file.replace('.gz', '')

        tmpdir = 'wcnTMP'

        # Extract all files to a folder
        with ZipFile(directory+'/'+wcn_file+'.gz', 'r') as zip_ref:
            zip_ref.extractall(tmpdir)

        # Check for WCN profile files and load them
        # Check AA
        wcn_aa_file = tmpdir+'/'+wcn_file+'_wcn_aa.npy'
        if os.path.exists(wcn_aa_file):
            if verbose:
                print('Found AA WCN profile')
            self.wcn_aa = np.load(wcn_aa_file)

        # Check CA
        wcn_ca_file = tmpdir+'/'+wcn_file+'_wcn_ca.npy'
        if os.path.exists(wcn_ca_file):
            if verbose:
                print('Found CA WCN profile')
            self.wcn_ca = np.load(wcn_ca_file)

        # Check average per residue
        wcn_pr_file = tmpdir+'/'+wcn_file+'_wcn_pr.npy'
        if os.path.exists(wcn_pr_file):
            if verbose:
                print('Found average-per-residue WCN profile')
            self.wcn_pr = np.load(wcn_pr_file)

        # Check per-chain
        for d in os.listdir(tmpdir):
            if '_wcn_per_chain_' in d:
                chains = d.split('_')[-1].replace('.npy','')
                if verbose:
                    print('Found per-chain WCN profile for chains: '+chains)
                self.wcn_per_chain[chains] = np.load(tmpdir+'/'+d)

        shutil.rmtree(tmpdir)

    def matchAtomsToTrajectory(self, trajectory, return_by_chains=True):
        """
        Matches the atoms between a trajectory and the WCNObject.structure. It assumes
        that the order of the chains between the objects is the same. The matching is
        based on a sequence alignemnt between the two objects, and then a matching between
        the atom names between aligned residues (it is not recommended to match
        hydrogens). The function returns two dictionares of indexes (traj_indexes, struture_indexes)
        to slice the respective per chain objects into arrays of common atoms. if
        return_by_chains is True, then only a single list is returned containing
        all the atoms (i.e., not separated by chains).

        Parameters
        ==========
        trajectory : md.Trajectory
            An mdtraj object containing the trajectory to match against
        return_by_chains : bool
            Return the indexes separated by chains

        Returns
        =======
        traj_indexes : dict (or list)
            A dictionary containing the trajectory atoms' indexes for each chain common
            to the WCN structure object. It is a list if return_by_chains is False.
        structure_indexes : dict (or list)
            A dictionary containing the structure atoms' indexes for each chain common
            to the trajectory object. It is a list if return_by_chains is False.
        """

        # Get the chains for the structure and the trajectory object
        chains = {}
        chains['structure'] = [c for c in self.chains]
        chains['trajectory'] = [*range(len([*trajectory.topology.chains]))]

        if len(chains['structure']) != len(chains['trajectory']):
            raise ValueError('The WCN Object and the trajectory have a different number of chains')

        sequences = {}
        traj_indexes = {}
        structure_indexes = {}

        # Iterate over all the chains in the structure and the trajectory
        for i,c in enumerate(chains['structure']):

            traj_indexes[c.id] = []
            structure_indexes[c.id] = []

            # Get the sequences for the structure and trajectory object
            sequences[c.id] = {}
            sequences[c.id]['trajectory'] = alignment.getTopologySequence(trajectory.topology, i)
            sequences[c.id]['structure'] = getChainSequence(c)

            # Do a sequence alignment between the two objects
            sa = mafft.multipleSequenceAlignment(sequences[c.id])

            # Match common residues between the WCN and trajectory objects
            i1 = 0
            i2 = 0
            t_res_indexes = []
            s_res_indexes = []
            for i in range(sa.get_alignment_length()):
                if sa[0].seq[i] != '-' and sa[1].seq[i] != '-':
                    t_res_indexes.append(i1)
                    s_res_indexes.append(i2)
                if sa[0].seq[i] != '-':
                    i1 += 1
                if sa[1].seq[i] != '-':
                    i2 += 1

            # Get common residues indexes as numpy arrays
            t_res_indexes = np.array(t_res_indexes)
            s_res_indexes = np.array(s_res_indexes)

            # Create list of residues in each object
            traj_residues = [r for r in trajectory.topology.residues]
            structure_residues = [r for r in c.get_residues()]

            # Match atoms by name
            for z in zip(t_res_indexes, s_res_indexes):
                tr = traj_residues[z[0]]
                sr = structure_residues[z[1]]
                t_atoms = [a for a in tr.atoms]
                s_atoms = [a for a in sr.get_atoms()]
                s_atoms_names = [a.name for a in sr.get_atoms()]

                for i,n in enumerate(t_atoms):
                    try:
                        j = s_atoms_names.index(n.name)
                    except:
                        j = None

                    # Save only indexes of matched atoms by name
                    if not isinstance(j, type(None)):
                        traj_indexes[c.id].append(t_atoms[i].index)
                        structure_indexes[c.id].append(self.all_atoms.index(s_atoms[j]))

        # Convert dictionaries into lists containing all indexes
        if not return_by_chains:
            ti = []
            si = []
            for c in traj_indexes:
                ti += traj_indexes[c]
                si += structure_indexes[c]
            traj_indexes = ti
            structure_indexes = si

        return traj_indexes, structure_indexes

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


def calculateRMSF(trajectory, c_alpha=False, normalize=True):

    # Align trajectory
    trajectory.superpose(trajectory[0])
    # Calculate RMSF
    rmsf = np.average(np.linalg.norm(trajectory.xyz-np.average(trajectory.xyz, axis=0), axis=2), axis=0)

    if c_alpha:
        atoms = []
        for a in trajectory.topology.atoms:
            if a.name == 'CA':
                atoms.append(a.index)
        atoms = np.array(atoms)
        rmsf = rmsf[atoms]

    #Normalize
    if normalize:
        rmsf = _normalize(rmsf)

    return rmsf
