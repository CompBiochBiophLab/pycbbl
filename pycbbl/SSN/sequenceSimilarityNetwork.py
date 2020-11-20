import os
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from pycbbl import alignment

class sequenceSimilarityNetwork:
    """
    Class to perform sequence similarity networks.

    Attributes
    ----------
    sequences : dict
        Input dictionary mapping each entry code to its sequence.
    codes : list
        Reference list containing the protein codes in the order they are in the
        similarity matrix.
    similarity_matrix : np.ndarray
        Array containing the similarity values between each pair of sequences.
    node_attributes
        Dictionary holding attributes add to each node.

    Methods
    -------
    plotSimilarityMatrix()
        Generates a plot of the similarity matrix.
    createSifNetworkFile()
        Creates a network sif format file for a specific threshold.
    addNodeAttribute()
        Add attributes to the nodes in the network.
    createNodeAttributesFile()
        Generates a file containing all the attributes for the nodes in the network.
    """

    def __init__(self, sequences, similarity_matrix=None, similarity_matrix_file='pid.matrix.npy',
                 write_pid_file=True, overwrite=False, verbose=False):
        """
        Initializes the sequenceSimilarityNetwork class by giving a dictionary
        containing the protein codes as keys and sequences as values. The methods
        calculates and stores a PID matrix in a numpy format file. If the pid_file
        path is found then the matrix is read from this file, otherwise the PID
        matrix is calculated and the content is stored in the pid_file path. Alternatively
        a precalcualted similarity matrix can be given in numpy format.

        Parameters
        ----------
        sequences : dict
            Dictionary mapping each entry protein code to its sequence.
        similarity_matrix : numpy.ndarray
            The similarity matrix for the input sequences as a numpy array.
        similarity_matrix_file : str
            Path to the similarity matrix file in numpy format.
        write_pid_file : bool
            Write the PID matrix to a file?
        overwrite : bool
            Forces the recalcualation of the PID matrix and overwrites the pid_file
            file.
        """

        # Create class variables
        self.sequences = sequences
        if not similarity_matrix_file.endswith('.npy'):
            similarity_matrix_file = similarity_matrix_file+'.npy'
        self.similarity_matrix_file = similarity_matrix_file
        self.node_attributes = {}
        # Store sequences and order of codes
        self.codes = [code for code in self.sequences]

        # Check input parameters
        if not isinstance(sequences, dict):
            raise ValueError('Input sequences should be a dictionary.')

        if similarity_matrix != None:
            self.similarity_matrix = similarity_matrix
            overwrite = True # Write the similarity matrix into a file
        else:
            # Check if the similarity matrix file exists
            if os.path.exists(self.similarity_matrix_file) and not overwrite:
                print('Similarity matrix file found! Reading the similarity matrix from '+
                self.similarity_matrix_file)
                self.similarity_matrix = np.load(self.similarity_matrix_file)
            else:
                # Calcualte the similarity matrix
                self.similarity_matrix = np.zeros((len(self.codes), len(self.codes)))
                # Calculate similarity matrix using blast
                sequences = [self.sequences[code] for code in self.codes]
                for i in range(len(self.codes)):
                    if verbose:
                        print('Calculating PID values for : '+self.codes[i])
                    self.similarity_matrix[i] = alignment.blast.calculatePIDs(sequences[i],
                                                sequences)
                overwrite = True # Write the similarity matrix into a file

        # Check problems with the similarity matrix
        if self.similarity_matrix.shape[0] != self.similarity_matrix.shape[1]:
            raise ValueError('The matrix in the file is not a square matrix')
        if self.similarity_matrix.shape[0] != len(self.codes):
            raise ValueError('The first matrix dimmension do not mathch the number of sequences')
        if self.similarity_matrix.shape[1] != len(self.codes):
            raise ValueError('The seconde matrix dimmension do not mathch the number of sequences')

        if overwrite:
            np.save(self.similarity_matrix_file, self.similarity_matrix)

    def createSifNetworkFile(self, threshold, output_file='network.sif', overwrite=None):
        """
        Generates a network file based on the given threshold.

        Parameters
        ----------
        threshold : float
            Similarity value to create the network edges.
        output_file : str
            Name of the network output file.
        overwrite : bool
            Whether to overwrite an existing network file.
        """

        if not output_file.endswith('.sif'):
            output_file = output_file+'.sif'

        if output_file == 'network.sif':
            output_file = output_file.replace('.sif', str(threshold)+'.sif')

        if os.path.exists(output_file) and not overwrite:
            raise ValueError('Network file exists. Give overwrite=True to replace it.')

        added = []
        with open(output_file, 'w') as nf:

            # Write edges above threshold
            for i in range(self.similarity_matrix.shape[0]):
                for j in range(self.similarity_matrix.shape[1]):
                    if j > i:
                        if self.similarity_matrix[i][j] >= threshold:
                            nf.write('%s 1.0 %s\n' % (self.codes[i], self.codes[j]))
                            if self.codes[i] not in added:
                                added.append(self.codes[i])
                            if self.codes[j] not in added:
                                added.append(self.codes[j])

            # Write remaining codes
            for code in self.codes:
                if code not in added:
                    nf.write('%s\n' % code)

    def addNodeAttribute(self, attribute_name, attribute_values, overwrite=False):
        """
        Add an attribute category to the proteins in the network. The input is a
        dictionary mapping the attribute values to the protein codes. All missing
        codes will be given a 'none' value. If the attribute name or any attribute
        value have a space " " character it will be replaced by an underscore "_"
        character.

        Parameters
        ----------
        attribute_name : str
            Name of the node attribute to be added.
        attribute_values : dict
            Dictionary mapping the protein codes to the attribute values.
        overwrite : bool
            Whether to overwrite the previous attribute values with the same attributes
            name.
        """

        if attribute_name not in self.node_attributes:
            self.node_attributes[attribute_name] = {}
        elif attribute_name in self.node_attributes and not overwrite:
            raise ValueError('Attribute name already found. Give overwrite=True to rewrite the values.')

        # Remove empty spaces from the attribute name
        if ' ' in attribute_name:
            print('Space " " found in attribute %s' % attribute_name)
            attribute_name = attribute_name.replace(' ', '_')
            print('The new attribute name will be %s' % attribute_name)

        # Check for nodes not present in the network.
        for code in attribute_values:
            if code not in self.codes:
                print('Code %s is not present in the network!' % code)

        # Iterate for codes
        for code in self.codes:
            if code not in attribute_values:
                self.node_attributes[attribute_name][code] = 'none'
            else:
                if  ' ' in attribute_values[code]:
                    print('Space " " found in attribute value %s' % attribute_values[code])
                    attribute_values[code] = attribute_values[code].replace(' ', '_')
                    print('The new attribute value will be %s' % attribute_values[code])
                self.node_attributes[attribute_name][code] = attribute_values[code]

    def createNodeAttributesFile(self, output_file, overwrite=None):
        """
        Create a node property file using the properties given to the sequenceSimilarityNetwork.

        Parameters
        ----------
        output_file : str
            Name of the node property output file.
        overwrite : bool
            Whether to overwrite an existing node attributes file.
        """

        if os.path.exists(output_file) and not overwrite:
            raise ValueError('Attribute file exists. Give overwrite=True to replace it.')

        with open(output_file, 'w') as af:

            # Write header line
            af.write('code\t'+'\t'.join(self.node_attributes)+'\n')

            for code in self.codes:
                attribute_line = ''
                for attribute in self.node_attributes:
                    attribute_line += self.node_attributes[attribute][code]+'\t'
                af.write('%s\t%s\n' % (code, attribute_line))

    def colorNodeFillByAttribute(self, attribute_name, overwrite=False):
        """
        Create a node attribute to color nodes based on the value of an attribute.
        The colors are limited to 10 different attribute names, they are colored from
        the most numerous to the less numerous. Further than 10 attribute names will
        not be assigned a color.

        Parameters
        ----------
        attribute_name : str
            Name of the attribute to use to assign colors
        """

        color_palette = [matplotlib.colors.rgb2hex(matplotlib.pyplot.cm.tab10(i)) for i in range(10)]
        names = [self.node_attributes[attribute_name][n] for n in self.node_attributes[attribute_name]]
        names_set = set(names)
        names_count = sorted([ (name,names.count(name)) for name in names_set], key=lambda x:x[1], reverse=True)
        colors = {}
        for i,name in enumerate(names_count):
            if i < 10:
                for code in self.codes:
                    if self.node_attributes[attribute_name][code] == name[0]:
                        colors[code] = color_palette[i]

        self.addNodeAttribute('fill_color', colors, overwrite=overwrite)

    def plotSimilarityMatrix(self):
        """
        Plots the similarity matrix of the given sequences
        """

        plt.matshow(self.similarity_matrix, vmin=0.0, vmax=1.0)
        cbar = plt.colorbar()
        cbar.set_label('Sequence identity')
        plt.xlabel('Sequence index')
        plt.ylabel('Sequence index')
