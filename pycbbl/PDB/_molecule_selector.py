import Bio.PDB as PDB
import os
import shutil

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

    pdb_list = []

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
                pdb_list.append(pdb_directory+'/'+code.upper()+'.pdb')

    # Rename files
    for f in os.listdir(pdb_directory):
        if f.endswith('.ent'):
            os.rename(pdb_directory+'/'+f,pdb_directory+'/'+f.replace('pdb','').upper().replace('.ENT','.pdb'))
            pdb_list.append(pdb_directory+'/'+f.replace('pdb','').upper().replace('.ENT','.pdb'))

    # Remove unnecesary folders created by Bio.PDB method
    shutil.rmtree('obsolete')

    return pdb_list


class moleculeSelector:

    def __init__(self, pdb_foler):
        self.pdb_paths = {}
        # for pdb in os.listdir(pdb_folder)
