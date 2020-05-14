from pycbbl.PDB import getChainSequence
from pycbbl.PDB import blast
from pycbbl import clustering
from collections import OrderedDict
import matplotlib.pyplot as plt
import Bio.PDB as PDB
import numpy as np
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
    sequence : OrderedDict
        Dictionary containing the sequences of each PDB protein chain ID.
    pidMatrix : np.array
        Matrix of PID values for all sequences.
    sequence_matrix_map : dict
        Dictionary mapping (pdb_code, chain_id) tuples to the position in the comparison
        square matrices.
    matrix_sequence_map : dict
        Dictionary mapping the positions in the comparison square matrices to the
        (pdb_code, chain_id) tuples.
    pid_cluster = pycbbl.clustering.hierarchical.clusterByDistance
        Object that contains pid clustering informations
    pid_file = str
        File to write or read PIDs from.

    Methods
    -------
    calculatePIDMatrix()
        Calculate pairwise PID values for all sequences in the PDBs
    plotSequenceLengthDistribution()
        Plot the distribution of chain sequence lenghts

    """
    def __init__(self, pdb_folder, exclude=[], pid_file=None):
        """
        Initialize molecule selector class. This stores the pdb paths, reads structures
        with the Bio.PDB.PDBParser and gather the sequences of each chain.

        Parameters
        ----------
        pdb_folder : str
            Path to the folder containing the PDB files
        exclude : str or list
            PDB code(s) to exclude from the analysis.
        save_pids : str
            Path to file to write pid values, if file exists pid values will be read
            from it.
        """

        # Define attributes
        self.pdb_paths = OrderedDict()
        self.structure = OrderedDict()
        self.sequence = OrderedDict()
        self.sequence_matrix_map = {}
        self.matrix_sequence_map = {}
        self.excluded = []
        self.pid_file = pid_file

        # Define attributes holders (Initialized as None)
        self.pidMatrix = None
        self.pid_cluster = None

        # Check input variables
        if isinstance(exclude, str):
            exclude = [exclude]
        if not isinstance(exclude, list):
            raise ValueError('Codes to exclude from the analysis must be given as a list.')
        if not self.pid_file.endswith('.npy'):
            self.pid_file = pid_file+'.npy'
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
        count = 0
        for pdb in self.structure:
            self.sequence[pdb] = OrderedDict()
            for chain in self.structure[pdb].get_chains():
                self.sequence[pdb][chain.id] = getChainSequence(chain)
                # Store mapping function to find elements in comparison matrix
                self.sequence_matrix_map[(pdb, chain.id)] = count
                self.matrix_sequence_map[(count)] = (pdb, chain.id)
                count += 1

    def calculatePIDMatrix(self, force_calculation=False):
        """
        Calculate matrix with pairwise comparison of PID values for all sequences.

        Parameters
        ----------
        force_calculation : bool (False)
            Force PID matrix recomputation even if it was already calculated.

        Returns
        -------
        pidMatrix : numpy.array
            Matrix of PID values for all sequences.
        """

        #Check whether matrix was already calculated
        if not isinstance(self.pidMatrix, np.ndarray) or force_calculation:

            # Iterate sequences using the matrix_sequence mapper
            sequences = []
            for index in self.matrix_sequence_map:
                (pdb, chain) = self.matrix_sequence_map[index]
                sequences.append(self.sequence[pdb][chain])

            # Check if pid_file was given
            if self.pid_file != None:
                # Check if pid_file exists
                if os.path.exists(self.pid_file):
                    # Load pid matrix from file
                    self.pidMatrix = np.load(self.pid_file)
                    # Check marix for squarness
                    if self.pidMatrix.shape[0] == self.pidMatrix.shape[1]:
                        # Check if elements in matrix is the same as sequences
                        if self.pidMatrix.shape[0] == len(sequences):
                            return self.pidMatrix

            self.pidMatrix = np.zeros((len(sequences), len(sequences)))
            for i in range(len(sequences)):
                pids = blast.calculatePIDs(sequences[i], sequences)
                self.pidMatrix[i] = pids

            # Check if pid_file was given
            if self.pid_file != None:
                # Save PIDs in pid_file
                np.save(self.pid_file, self.pidMatrix)

        return self.pidMatrix

    def clusterByPID(self, clustering_distance, verbose=False):
        """
        Calculate clusters based on PID distance matrix

        Parameters
        ----------
        clustering_distance : float
            PID value to cluster structures
        import pdb; pdb.set_trace()
        verbose : bool (False)
            Print clustering information to screen
        """

        if verbose:
            print('Calculating PID matrix')
        self.calculatePIDMatrix()
        if verbose:
            print('Calculating clusters with hierarchical clustering')
        self.pid_cluster = clustering.hierarchical.clusterByDistance(self.pidMatrix,
                           verbose=verbose)

    def plotDendrogram(self, dpi=100):
        """
        Plot dendogram containing the clustering information

        Parameters
        ----------
        dpi : int
            Resolution of the generated plot.
        """
        if self.pid_cluster != None:
            self.pid_cluster.plotDendrogram(dpi=dpi)
        else:
            raise ValueError('First you must calculate the PID clusters with clusterByPID()')

    def plotSequenceLengthDistribution(self, dpi=100, bins=None):
        """
        This method plots the sequence length distribution of the chains in the
        PDB file.

        Parameters
        ----------
        dpi : int
            Resolution of the generated plot.
        """

        sequence_lengths = []

        for pdb in self.sequence:
            for chain in self.sequence[pdb]:
                sequence_lengths.append(len(self.sequence[pdb][chain]))

        # Define bins as the single spaced sequence range
        if bins == None:
            bins = range(min(sequence_lengths), max(sequence_lengths))
        figure = plt.figure(dpi=dpi)
        hist = plt.hist(sequence_lengths, bins=bins)
        plt.xlabel('Sequence length')
        plt.ylabel('Nº of polypeptides')
        plt.xlim(min(sequence_lengths)-2, max(sequence_lengths)+2)
