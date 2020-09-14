import os
import numpy as np
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

    """

    def __init__(self, sequences):
        """
        Initializes the sequenceSimilarityNetwork class by giving a dictionary
        containing the protein codes as keys and sequences as values. The methods
        calculates and stores the PID matrix.

        Parameters
        ----------
        sequences : dict
            Dictionary mapping each entry protein code to its sequence.
        """

        # Check input parameters
        if not isinstance(sequences, dict):
            raise ValueError('Input sequences should be a dictionary.')

        # Store sequences and order of codes
        self.sequences = sequences
        self.codes = [code for code in self.sequences]

        # Calculate similarity matrix
        sequences = [sequences[c] for c in sequences ]
        self.similarity_matrix = np.zeros((len(self.codes), len(self.codes)))
        for i in range(len(self.codes)):
            self.similarity_matrix[i] = alignment.blast.calculatePIDs(sequences[i],
                                        sequences)

        # Create attributes
        self.node_attributes = {}

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
            print('The new attribute name will be %' % attribute_name)

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
                    print('The new attribute value will be %' % attribute_values[code])
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
            af.write('code '+' '.join(self.node_attributes)+'\n')

            for code in self.codes:
                attribute_line = ''
                for attribute in self.node_attributes:
                    attribute_line += self.node_attributes[attribute][code]+' '
                af.write('%s %s\n' % (code, attribute_line))

    def plotSimilarityMatrix(self):
        """
        Plots the similarity matrix of the given sequences
        """

        plt.matshow(self.similarity_matrix, vmin=0.0, vmax=1.0)
        cbar = plt.colorbar()
        cbar.set_label('Sequence identity')
        plt.xlabel('Sequence index')
        plt.ylabel('Sequence index')
