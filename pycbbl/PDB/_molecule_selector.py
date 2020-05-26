import pycbbl
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
    sequence_lengths : dict
        Lenght of the sequences of each PDB protein chain ID.
    pidMatrix : np.array
        Matrix of PID values for all sequences.
    sequence_matrix_map : dict
        Dictionary mapping (pdb_code, chain_id) tuples to the position in the comparison
        square matrices.
    matrix_sequence_map : dict
        Dictionary mapping the positions in the comparison square matrices to the
        (pdb_code, chain_id) tuples.
    pid_clustering = pycbbl.clustering.hierarchical.clusterByDistance
        Object that contains pid clustering informations.
    pid_clusters = dict
        Dictionary containing the elements of each PID cluster or only the centroids.
    pid_file = str
        File to write or read PIDs from.
    selection : Dictionary containing selections

    Methods
    -------
    calculatePIDMatrix()
        Calculate pairwise PID values for all sequences in the PDBs
    clusterByPID()
        Cluster PDBs by PID values using hierarchical clustering.
    plotDendrogram()
        Plot a dendrogram based on the hierarchical clustering.
    selectBySequenceLengthDistribution()
        Apply selection based on the distribution of chain sequence lenghts
    selectByReferencePID()
        Apply selection based on the distribution of PID regarding a reference pdb
        chain.
    plotSequenceLengthVsAveragePID
        Plot a scatter plot to visualize the lenght of the pdb chain sequences vs
        the average PID of each row in the PID matrix.
    plotSequenceLengthVsReferencePID
        Plot a scatter plot to visualize the lenght of the pdb chain sequences vs
        the PID values to a reference PDB chain (a specific row in the PID matrix).
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
        self.sequence_lengths = {}
        self.sequence_matrix_map = {}
        self.matrix_sequence_map = {}
        self.excluded = []
        self.pid_file = pid_file
        self.selection = {}
        self.last_selection = 0

        # Define attributes holders (Initialized as None)
        self.pidMatrix = None
        self.pid_clustering = None



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
            self.sequence_lengths[pdb] = OrderedDict()
            for chain in self.structure[pdb].get_chains():
                # Store sequence and its length
                self.sequence[pdb][chain.id] = pycbbl.PDB.methods.getChainSequence(chain)
                self.sequence_lengths[pdb][chain.id] = len(self.sequence[pdb][chain.id])
                # Store mapping function to sort and find elements in comparison matrix
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

    def clusterByPID(self, clustering_distance, verbose=False, return_centroids=False):
        """
        Calculate clusters based on PID distance matrix. Please plot dendrogram
        to estimate the threshold values.

        Parameters
        ----------
        clustering_distance : float
            PID value to cluster structures.
        verbose : bool (False)
            Print clustering information to screen.
        return_centroids : bool (False)
            Return only the centroids of each cluster.
        """

        if verbose:
            print('Calculating PID matrix')
        # Calcualte PID matrix
        self.calculatePIDMatrix()
        if verbose:
            print('Calculating clusters at threshold distance %.2f' % clustering_distance)
        # Calculate clusters
        self.pid_clustering = clustering.hierarchical.clusterByDistance(self.pidMatrix,
                                                                       verbose=verbose)

        self.pid_clusters = self.pid_clustering.getClusters(clustering_distance,
                                              return_centroids=return_centroids)

        return self.pid_clusters

    def plotDendrogram(self, dpi=100):
        """
        Plot dendogram containing the clustering information

        Parameters
        ----------
        dpi : int
            Resolution of the generated plot.
        """
        if self.pid_clustering != None:
            self.pid_clustering.plotDendrogram(dpi=dpi)
        else:
            raise ValueError('First you must calculate the PID clusters with clusterByPID()')

    def selectBySequenceLengthDistribution(self, dpi=100, bins=None, vertical_line=None,
        select=0, new_selection=False, apply_to_selection=None):
        """
        This method allows to select the left (select=1) or the right side (select=2)
        of the sequences lengths distribution. Without keywords it generates a plot
        of the distribution. Adding the keyword vertical_line divides the distribution
        into left and right side. If, additionally, the keyword select is given,
        a new selection is created based on the side of the distribution selected.
        This selection is stored in the attribute "selection" of the "moleculeSelector"
        class and the previous selection is overwritten. To avoid overwritting the
        old selection the option new_selection must be set tot True, creating a
        new selection entry in the attribute "selection". The default value of select
        is 0 and means that no selection will be carried out. Additionally, if apply_to_selection
        option is used, only PDB chains pertaining to that particular selection
        will be considered.

        Parameters
        ----------
        dpi : int
            Resolution of the generated plot.
        bins : range
            Bins to use in the histogram of the distribution.
        vertical_line : float
            Sequence length value at which to divide the distribution.
        select : int
            Side of the distribution to select: 0 left and 1 for right side.
        apply_to_selection : int
            Use a previous selection of models stored in the selection attribute
            of this class.
        """

        if not isinstance(select, int):
            raise ValueError('select must be an integer. Please see documentation.')

        selection = set()

        # Divide values into two distributions by vertical_line
        if vertical_line != None:
            ls_sequence_lengths = []
            rs_sequence_lengths = []
            for pdb in self.sequence:
                for chain in self.sequence[pdb]:
                    i = self.sequence_matrix_map[(pdb,chain)]
                    if apply_to_selection != None:
                        if apply_to_selection in self.selection:
                            if i not in self.selection[apply_to_selection]:
                                continue
                        else:
                            raise ValueError('Selection given in "apply_to_selection" not found\
                            Please print moleculeSelector.selection attriubute to see available selections.')
                    sl = self.sequence_lengths[pdb][chain]
                    if sl <= vertical_line:
                        ls_sequence_lengths.append(sl)
                        if select == 1:
                            i = self.sequence_matrix_map[(pdb,chain)]
                            selection.add(i)
                    else:
                        rs_sequence_lengths.append(sl)
                        if select == 2:
                            i = self.sequence_matrix_map[(pdb,chain)]
                            selection.add(i)

            # Save selection
            if select > 0:
                if new_selection:
                    self.last_selection += 1
                self.selection[self.last_selection] = selection

        else:
            # Gather the values of sequence lengths
            sequence_lengths = []
            for pdb in self.sequence:
                for chain in self.sequence[pdb]:
                    i = self.sequence_matrix_map[(pdb,chain)]
                    if apply_to_selection != None:
                        if apply_to_selection in self.selection:
                            if i not in self.selection[apply_to_selection]:
                                continue
                        else:
                            raise ValueError('Selection given in "apply_to_selection" not found\
                            Please print moleculeSelector.selection attriubute to see available selections.')
                    sequence_lengths.append(self.sequence_lengths[pdb][chain])

        # Plot sequences length distribution
        figure = plt.figure(dpi=dpi)
        plt.xlabel('Sequence length')
        plt.ylabel('Nº of polypeptides')

        if vertical_line != None:
            # Define bins as the single spaced sequence range
            if bins == None:
                bins = range(min(ls_sequence_lengths), max(rs_sequence_lengths))
            hist = plt.hist(ls_sequence_lengths, bins=bins, color='r')
            hist = plt.hist(rs_sequence_lengths, bins=bins, color='b')
            plt.axvline(vertical_line, ls='--', c='k')
            plt.xlim(min(ls_sequence_lengths)-2, max(rs_sequence_lengths)+2)
            if select == 1:
                plt.axvspan(min(ls_sequence_lengths)-2, vertical_line, facecolor='g', alpha=0.25)
            elif select == 2:
                plt.axvspan(vertical_line, max(rs_sequence_lengths)+2, facecolor='g', alpha=0.25)
        else:
            if bins == None:
                bins = range(min(sequence_lengths), max(sequence_lengths))
            hist = plt.hist(sequence_lengths, bins=bins)
            plt.xlim(min(sequence_lengths)-2, max(sequence_lengths)+2)

    def selectByReferencePID(self, reference, dpi=100, vertical_line=None, select=0, new_selection=False, apply_to_selection=None):
        """
        This method allows to select the left (select=1) or the right side (select=2)
        of the PID distribution relative to a reference PDB chain. Without keywords
        it generates a plot of the distribution. Adding the keyword vertical_line
        divides the distribution into left and right side. If, additionally, the
        keyword select is given, a new selection is created based on the side of
        the distribution selected. This selection is stored in the attribute
        "selection" of the "moleculeSelector" class and the previous selection is
        overwritten. To avoid overwritting the old selection the option new_selection
        must be set tot True, creating a new selection entry in the attribute "selection".
        The default value of select is 0 and means that no selection will be carried
        out. Additionally, if apply_to_selection option is used, only PDB chains
        pertaining to that particular selection will be considered.

        Parameters
        ----------
        reference : tuples
            2 elements-tuple containing the pdb (str) and chain id (str) of the reference PDB chain.
        dpi : int
            Resolution of the generated plot.
        vertical_line : float
            Sequence length value at which to divide the distribution.
        select : int (0)
            Side of the distribution to select: 0 left and 1 for right side.
        apply_to_selection : int
            Use a previous selection of models stored in the selection attribute
            of this class.
        """

        if not isinstance(select, int):
            raise ValueError('select must be an integer. Please see documentation.')

        selection = set()

        # Define reference index
        reference_index = self.sequence_matrix_map[(reference[0], reference[1])]

        # Divide values into two distributions by vertical_line
        if vertical_line != None:
            ls_reference_PIDs = []
            rs_reference_PIDs = []

            # Gather reference PIDs values
            for i in range(self.pidMatrix.shape[0]):
                if apply_to_selection != None:
                    if apply_to_selection in self.selection:
                        if i not in self.selection[apply_to_selection]:
                            continue
                    else:
                        raise ValueError('Selection given in "apply_to_selection" not found\
                        Please print moleculeSelector.selection attriubute to see available selections.')
                rpid = self.pidMatrix[i][reference_index]
                if rpid <= vertical_line:
                    ls_reference_PIDs.append(rpid)
                    if select == 1:
                        selection.add(i)
                else:
                    rs_reference_PIDs.append(rpid)
                    if select == 2:
                        selection.add(i)

            if select > 0:
                if new_selection:
                    self.last_selection += 1
                self.selection[self.last_selection] = selection

        else:
            # Gather the values of sequence lengths
            reference_PIDs = []
            for i in range(self.pidMatrix.shape[0]):
                if apply_to_selection != None:
                    if apply_to_selection in self.selection:
                        if i not in self.selection[apply_to_selection]:
                            continue
                    else:
                        raise ValueError('Selection given in "apply_to_selection" not found\
                        Please print moleculeSelector.selection attriubute to see available selections.')
                rpid = self.pidMatrix[i][reference_index]
                reference_PIDs.append(rpid)

        # Plot sequences length distribution
        figure = plt.figure(dpi=dpi)
        plt.xlabel('PID to chain %s of pdb %s' % (reference[1], reference[0]))
        plt.ylabel('Nº of polypeptides')
        bins = [i/100 for i in range(0,100)]
        plt.xlim(0, 1)

        if vertical_line != None:
            hist = plt.hist(ls_reference_PIDs, bins=bins, color='r')
            hist = plt.hist(rs_reference_PIDs, bins=bins, color='b')
            plt.axvline(vertical_line, ls='--', c='k')
            if select == 1:
                plt.axvspan(0, vertical_line, facecolor='g', alpha=0.25)
            elif select == 2:
                plt.axvspan(vertical_line, 1, facecolor='g', alpha=0.25)
        else:
            hist = plt.hist(reference_PIDs, bins=bins)

    def selectByReferenceRMSD(self, reference, models_folder, dpi=100, vertical_line=None,
        select=0, new_selection=False, apply_to_selection=None, align=False):
        """
        This method allows to select the left (select=1) or the right side (select=2)
        of the RMSD distribution relative to a reference PDB chain. The method works
        with previously saved PDB chains, so it needs the path to the folder where
        those models were saved. Without keywords it generates a plot of the distribution.
        Adding the keyword vertical_line divides the distribution into left and
        right side. If, additionally, the keyword select is given, a new selection
        is created based on the side of the distribution selected. This selection
        is stored in the attribute "selection" of the "moleculeSelector" class and
        the previous selection is overwritten. To avoid overwritting the old selection
        the option new_selection must be set tot True, creating a new selection
        entry in the attribute "selection". The default value of select is 0 and
        means that no selection will be carried out. Additionally, if apply_to_selection
        option is used, only PDB chains pertaining to that particular selection
        will be considered.

        Parameters
        ----------
        reference : tuples
            2 elements-tuple containing the pdb (str) and chain id (str) of the reference PDB chain.
        models_folder : str
            Path to the folder where the PDB chains are stored.
        dpi : int
            Resolution of the generated plot.
        vertical_line : float
            Sequence length value at which to divide the distribution.
        select : int (0)
            Side of the distribution to select: 0 left and 1 for right side.
        apply_to_selection : int
            Use a previous selection of models stored in the selection attribute
            of this class.
        align : bool
            Whether to align models to the reference chain.
        """

        if not isinstance(select, int):
            raise ValueError('select must be an integer. Please see documentation.')

        selection = set()

        # Define reference index
        reference_index = self.sequence_matrix_map[(reference[0], reference[1])]

        # Read PDBs in given folder
        models = {}
        for m in os.listdir(models_folder):
            if m.endswith('.pdb'):
                pdb = m.split('.')[0].split('_')[0]
                chain = m.split('.')[0].split('_')[1]
                models[(pdb, chain)] = m

        # Check if reference model is in folder
        if not (reference[0], reference[1]) in models:
            raise ValueError('Reference model was not found in given folder.')

        # Gather the RMSD values
        pdbs = [models[(reference[0], reference[1])]]
        for i in range(self.pidMatrix.shape[0]):
            pdb, chain = self.matrix_sequence_map[i]
            if (pdb, chain) in models:
                if apply_to_selection != None:
                    if apply_to_selection in self.selection:
                        if i not in self.selection[apply_to_selection]:
                            continue
                    else:
                        raise ValueError('Selection given in "apply_to_selection" not found\
                        Please print moleculeSelector.selection attriubute to see available selections.')
                if (pdb, chain) in models:
                    pdbs.append(models[(pdb, chain)])

        alignment = pycbbl.PDB.pymol.alignPDBs(models_folder, pdb_list=pdbs, align=align)
        reference_RMSDs = [alignment[n]['RMS'] for n in pdbs[1:]]

        # Divide values into two distributions by vertical_line
        if vertical_line != None:
            ls_reference_RMSDs = []
            rs_reference_RMSDs = []
            for i in range(self.pidMatrix.shape[0]):
                pdb, chain = self.matrix_sequence_map[i]
                if (pdb, chain) in models:
                    if apply_to_selection != None:
                        if apply_to_selection in self.selection:
                            if i not in self.selection[apply_to_selection]:
                                continue
                        else:
                            raise ValueError('Selection given in "apply_to_selection" not found\
                            Please print moleculeSelector.selection attriubute to see available selections.')

                    rmsd = alignment[models[(pdb, chain)]]['RMS']
                    if rmsd <= vertical_line:
                        ls_reference_RMSDs.append(rmsd)
                        if select == 1:
                            selection.add(i)
                    else:
                        rs_reference_RMSDs.append(rmsd)
                        if select == 2:
                            selection.add(i)

                if select > 0:
                    if new_selection:
                        self.last_selection += 1
                    self.selection[self.last_selection] = selection

        # Plot sequences length distribution
        figure = plt.figure(dpi=dpi)
        plt.xlabel('RMSD to chain %s of pdb %s' % (reference[1], reference[0]))
        plt.ylabel('Nº of polypeptides')
        bins = [i/2 for i in range(0, int(max(reference_RMSDs))*2+3)]
        plt.xlim(0, max(reference_RMSDs)+1)

        if vertical_line != None:
            hist = plt.hist(ls_reference_RMSDs, bins=bins, color='r')
            hist = plt.hist(rs_reference_RMSDs, bins=bins, color='b')
            plt.axvline(vertical_line, ls='--', c='k')
            if select == 1:
                plt.axvspan(0, vertical_line, facecolor='g', alpha=0.25)
            elif select == 2:
                plt.axvspan(vertical_line, 1, facecolor='g', alpha=0.25)
        else:
            hist = plt.hist(reference_RMSDs, bins=bins)

    def getSelectionTuples(self, selection):
        """
        Get the tuples (pdb_code,chain) for a specific selection index.

        Parameters
        ----------
        selection : int
            Selection index to use
        """

        models = []

        for i in self.selection[selection]:
            models.append(self.matrix_sequence_map[i])

        return models

    def saveSelection(self, selection, output_folder):
        """
        Save separate PDB for each of the selected chains

        Parameters
        ----------
        selection : int
            Selection index to use
        output_folder : str
            Path to the folder to store models
        """

        if not os.path.exists(output_folder):
            os.mkdir(output_folder)

        for m in self.getSelectionTuples(selection):
            for chain in self.structure[m[0]].get_chains():
                if chain.id == m[1]:
                    structure = pycbbl.PDB.methods.chainAsStructure(chain)
                    output_path = output_folder+'/'+m[0]+'_'+m[1]+'.pdb'
                    pycbbl.PDB.methods.saveStructureToPDB(structure,  output_path)

    def plotSequenceLengthVsAveragePID(self, dpi=100, vertical_line=None, horizontal_line=None, **kwargs):
        """
        This method plots the sequence length distribution vs the average PIDs of
        the chains in the PDB files.

        Parameters
        ----------
        dpi : int
            Resolution of the generated plot.
        vertical_line : float
            Value to plot a vertical dashed line
        horizontal_line : float
            Value to plot a horizontal dashed line
        **kwargs
            Keywords for the matplotlib.pyplot.scatter method.
        """
        sequence_lengths = []
        average_PIDs = []

        # Gather sequence lengths and average PIDs values
        for i in range(self.pidMatrix.shape[0]):
            (pdb,chain) = self.matrix_sequence_map[i]
            sequence_lengths.append(self.sequence_lengths[pdb][chain])
            average_PIDs.append(np.average(self.pidMatrix[i]))

        # Plot sequence lengths vs average PIDs values
        figure = plt.figure(dpi=dpi)
        plt.scatter(sequence_lengths, average_PIDs, **kwargs)
        plt.xlabel('Sequence length')
        plt.ylabel('Average PID')
        if horizontal_line != None:
            plt.axhline(horizontal_line, ls='--')
        if vertical_line != None:
            plt.axvline(vertical_line, ls='--')
        plt.xlabel('Sequence length')
        plt.ylabel('Average PID')

    def plotSequenceLengthVsReferencePID(self, reference, dpi=100, vertical_line=None, horizontal_line=None, **kwargs):
        """
        This method plots the sequence length distribution vs PIDs to a reference of
        the chains in the PDB files

        Parameters
        ----------
        reference : tuples
            2 elements-tuple containing the pdb (str) and chain id (str) of the reference PDB chain.
        dpi : int
            Resolution of the generated plot.
        vertical_line : float
            Value to plot a vertical dashed line
        horizontal_line : float
            Value to plot a horizontal dashed line
        **kwargs
            Keywords for the matplotlib.pyplot.scatter method.
        """

        sequence_lengths = []
        reference_PIDs = []
        reference_index = self.sequence_matrix_map[(reference[0], reference[1])]

        # Gather sequence lengths and reference PIDs values
        for i in range(self.pidMatrix.shape[0]):
            (pdb,chain) = self.matrix_sequence_map[i]
            sequence_lengths.append(self.sequence_lengths[pdb][chain])
            reference_PIDs.append(self.pidMatrix[i][reference_index])

        # Plot sequence lengths vs reference PIDs values
        figure = plt.figure(dpi=dpi)
        plt.scatter(sequence_lengths, reference_PIDs, **kwargs)
        plt.xlabel('Sequence length')
        plt.ylabel('Average PID')
        if horizontal_line != None:
            plt.axhline(horizontal_line, ls='--')
        if vertical_line != None:
            plt.axvline(vertical_line, ls='--')
        plt.xlabel('Sequence length')
        plt.ylabel('PID to chain %s of pdb %s' % (reference[1], reference[0]))
