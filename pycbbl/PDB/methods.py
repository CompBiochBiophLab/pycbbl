from Bio import PDB
from Bio import AlignIO
import mdtraj as md
import os
import re
import shutil
import numpy as np
from collections import OrderedDict

"""
Collection of general methods to work with PDB files.

Methods
-------
retrievePDBs()
    Download PDB files from the PDB database.

"""

def retrievePDBs(pdb_codes, pdb_directory='PDB'):
    """
    Download a set of pdb structures from the PDB database.

    Parameters
    ----------
    pdb_codes : list or str
        A list with the PDB codes to retrieve their files. Optionally a string
        with a unique code can be given.
    pdb_directory : str
        The name of the directory to store the retrieved files.

    Returns
    -------
    pdb_paths : list
        A list containing the paths to the retrieved PDB files.
    """

    pdb_paths = []

    pdbl = PDB.PDBList()

    # Create directory to store files
    if not os.path.exists(pdb_directory):
        os.mkdir(pdb_directory)

    # Convert to list if only a string is given
    if isinstance(pdb_codes, str):
        pdb_codes = [pdb_codes]

    if isinstance(pdb_codes, list):
        # Iterate pdb codes
        for code in pdb_codes:
            # Download PDB file
            if not os.path.exists(pdb_directory+'/'+code.upper()+'.pdb'):
                pdbl.retrieve_pdb_file(code, file_format='pdb', pdir=pdb_directory)
            else: # If file already exists
                print('Structure exists: '+code.upper()+'.pdb')
                pdb_paths.append(pdb_directory+'/'+code.upper()+'.pdb')

    for f in os.listdir(pdb_directory):
        if f.endswith('.ent'):
            # Rename file
            os.rename(pdb_directory+'/'+f,
            pdb_directory+'/'+f.replace('pdb','').upper().replace('.ENT','.pdb'))
            # Append path
            pdb_paths.append(pdb_directory+'/'+f.replace('pdb','').upper().replace('.ENT','.pdb'))

    # Remove unnecesary folders created by Bio.PDB method
    shutil.rmtree('obsolete')

    return pdb_paths

def getPDBpaths(pdb_folder, exclude=[], return_excluded=False):
    """
    Get all files with extension .pdb from a specific folder.

    Parameters
    ----------
    pdb_folder : str
        Path to the folder containing the PDB files
    exclude : str or list
        PDB file or name (wihout extension) to exclude from the analysis.

    Returns
    -------
    pdb_paths : dict
        Dictionary that matches the pdb codes to the path of the corresponding pdb
        files.
    """

    # Check input variables
    if isinstance(exclude, str):
        exclude = [exclude]
    if not isinstance(exclude, list):
        raise ValueError('PDB files to exclude from the analysis must be given as a list.')

    # Remove .pdb extension form excluded files
    exclude = [f.replace('.pdb','') for f in exclude]

    # Add '/' to end of pdb_folder if not given
    if not pdb_folder.endswith('/'):
        pdb_folder = pdb_folder+'/'

    pdb_paths = {}
    excluded = []

    # Store PDB paths
    for pdb in sorted(os.listdir(pdb_folder)):
        if pdb.endswith('.pdb'):
            name = pdb.replace('.pdb', '')
            if name not in exclude:
                pdb_paths[name] = pdb_folder+pdb
            elif return_excluded:
                excluded.append(name)

    if return_excluded:
        return pdb_paths, excluded
    else:
        return pdb_paths

def getChainSequence(chain):
    """
    Get the one-letter protein sequence of a Bio.PDB.Chain object.

    Parameters
    ----------
    chain : Bio.PDB.Chain
        Input chain to retrieve its sequence from.

    Returns
    -------
    sequence : str
        Sequence of the input protein chain.
    None
        If chain does not contain protein residues.
    """
    sequence = ''
    for r in chain:
        if r.id[0] == ' ': # Non heteroatom filter
            try:
                sequence += PDB.Polypeptide.three_to_one(r.resname)
            except:
                sequence += 'X'
    if sequence == '':
        return None
    else:
        return sequence

def chainsAsStructure(chains):
    """
    This method creates a new Structure object containing only the given chains.

    Parameters
    ----------
    chains : list or Bio.PDB.Chain
        Chain or chains to be added to the new structure object.

    Returns
    -------
    structure : Bio.PDB.Structure
    """

    if not isinstance(chains, list):
        chains = [chains]

    structure = PDB.Structure.Structure(0)
    model = PDB.Model.Model(0)
    for chain in chains:
        model.add(chain)
    structure.add(model)

    return structure

def renumberResidues(structure, by_chain=False):
    """
    Renumber residues in a structure object starting from one. Two methods are possible:
    if by_chain is set to True the renumbering is restarted at the begining of every
    chain.

    Parameters
    ----------
    structure : Bio.PDB.Structure
        Input structure object
    by_chain : bool
        Whether each chain should be renumerated from one.

    Returns
    -------

    structure_copy : Bio.PDB.Structure
        Renumbered copy of the input structure
    """

    count = 0
    structure_copy = PDB.Structure.Structure(0)
    model = PDB.Model.Model(0)
    auxiliar = PDB.Chain.Chain(0)

    for chain in structure.get_chains():
        new_chain = PDB.Chain.Chain(chain.id)
        if by_chain:
            count = 0
        for residue in chain.get_residues():
            count += 1
            residue.set_parent(auxiliar)
            residue.id = (residue.id[0], count, residue.id[2])
            new_chain.add(residue)
        model.add(new_chain)
    structure_copy.add(model)

    return structure_copy

class notHydrogen(PDB.Select):
    def accept_atom(self, atom):
        """
        Verify if atom is not Hydrogen.
        """
        _hydrogen = re.compile("[123 ]*H.*")
        name = atom.get_id()
        if _hydrogen.match(name):
            return 0
        else:
            return 1

class notWater(PDB.Select):
    def accept_residue(self, residue):
        """
        Verify if residue is water.
        """
        _restype = residue.id[0]
        if _restype == 'W':
            return 0
        else:
            return 1

class onlyProtein(PDB.Select):

    def accept_residue(self, residue):
        """
        Verify if residues are protein.
        """
        _restype = residue.id[0]
        if _restype != ' ':
            return 0
        else:
            return 1

def saveStructureToPDB(structure, output, remove_hydrogens=False, remove_water=False,
                        only_protein=False):
    """
    Saves a structure into a PDB file

    Parameters
    ----------
    structure : list or Bio.PDB.Structure
        Structure to save
    """

    io = PDB.PDBIO()
    io.set_structure(structure)

    selector = None
    if remove_hydrogens:
        selector = notHydrogen()
    elif remove_water:
        selector = notWater()
    elif only_protein:
        selector = onlyProtein()

    if selector != None:
        io.save(output, selector)
    else:
        io.save(output)

def readPDB(pdb_file):
    """
    Reads a PDB file with the Biopython PDB parser.

    Parameters
    ----------
    pdb_file : str
        path to the pdb file.

    Returns
    -------
    structure : list or Bio.PDB.Structure
        Structure objects
    """

    parser = PDB.PDBParser()
    name = pdb_file.split('/')[-1].split('.pdb')[0]
    structure = parser.get_structure(name, pdb_file)

    return structure

def PDBsToTrajectory(folder, return_filenames=False, verbose=False):
    """
    Method to create a trajectory (mdtraj object) from many PDB files inside a specific
    folder. It is mandatory that the PDB models contain the same number of atoms,
    otherwise the method will fail. The PDB names are sorted before reading.

    Parameters
    ----------
    folder : str
        path to the folder containing the pdb files.
    return_filenames : bool
        If True the function will return a 2-tuple with the trajectory as the first
        element and a list containing the file names as second element (no pdb extension).
    verbose : bool
        Print wich PDBs is being process.

    Returns
    -------
    trajectory : mdtraj.Trajectory
        Trajectory object containing the PDBs coordinates and topology.

    if return_filenames == True
        (trajectory, file_names): tuple
            Tuple containing the trajectory and a list with the file names.
    """

    topology = None
    trajectory = None
    file_names = []
    for f in sorted(os.listdir(folder)):
        if f.endswith('.pdb'):
            if verbose:
                print('Processing file: '+f, end = "\r")
            traj = md.load(folder+'/'+f)
            file_names.append(f.replace('.pdb', ''))
            if isinstance(topology, type(None)):
                topology = traj.topology
            if isinstance(trajectory, type(None)):
                trajectory = traj
            else:
                if trajectory.xyz.shape[1:] != traj.xyz.shape[1:]:
                    print('PDB file: '+folder+'/'+f+' has different number of atoms \
than the reference topology! Discarding it')
                    continue
                traj = md.Trajectory(traj.xyz, topology)
                trajectory = md.join((trajectory, traj))

    if return_filenames:
        return (trajectory, file_names)
    else:
        return trajectory
