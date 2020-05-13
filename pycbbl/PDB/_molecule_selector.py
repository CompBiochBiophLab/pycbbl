from pycbbl.PDB import getChainSequence
from collections import OrderedDict
import Bio.PDB as PDB
import os

class moleculeSelector:
    """
    Molecule selector is a class that allows to find chains for a group of similar
    proteins in a set of PDB files.

    Attributes
    ----------
    pdb_paths : OrderedDict
        Dictionary that matches the pdb codes to the path of the corresponding pdb
        files.
    excluded : list
        List containing the codes excluded from the analysis.
    structure : OrderedDict
        Dictionary containing Bio.PDB.Structure objects for each pdb file.
    sequence : OrderedDict()
        Dictionary containing the sequences of each PDB protein chain ID.

    Methods
    -------

    """
    def __init__(self, pdb_folder, exclude=[]):
        """
        Initialize molecule selector class. This stores the pdb paths, reads structures
        with the Bio.PDB.PDBParser and gather the sequences of each chain.

        Parameters
        ----------
        pdb_folder : str
            Path to the folder containing the PDB files
        exclude : str or list
            PDB code(s) to exclude from the analysis.
        """

        # Define attributes
        self.pdb_paths = OrderedDict()
        self.structure = OrderedDict()
        self.sequence = OrderedDict()
        self.excluded = []

        # Check input variables
        if isinstance(exclude, str):
            exclude = [exclude]
        if not isinstance(exclude, list):
            raise ValueError('Codes to exclude from the analysis must be given as a list.')

        # Store PDB paths
        for pdb in sorted(os.listdir(pdb_folder)):
            if pdb.endswith('.pdb'):
                code = pdb.replace('.pdb','')
                if code not in exclude:
                    self.pdb_paths[code] = pdb_folder+'/'+pdb
                else: # Exclude files
                    self.excluded.append(code)

        # Read structures
        parser = PDB.PDBParser()
        for pdb in self.pdb_paths:
            self.structure[pdb] = parser.get_structure(pdb, self.pdb_paths[pdb])

        # Generate sequences for each PDB protein chain
        for pdb in self.structure:
            self.sequence[pdb] = {}
            for chain in self.structure[pdb].get_chains():
                self.sequence[pdb][chain.id] = getChainSequence(chain)
