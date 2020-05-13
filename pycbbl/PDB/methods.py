from Bio import PDB
import os
import shutil

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

def getChainSequence(chain):
    """
    Get the one-letter protein sequence of a Bio.PDB.Chain object.

    Attributes
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
