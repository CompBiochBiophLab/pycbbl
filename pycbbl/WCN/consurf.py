import os
import numpy as np

class consurf:
    """
    Class for reading and holding evolutionary sequence conservation from the server
    Consurf.

    Attributes
    ----------
    consurf_file : str
        Path to the file of residue conservation data.

    Methods
    -------


    """
    def __init__(self, consurf_file):
        """
        Initializes consurf class. The input requires the file that contains the
        aminoacid conservation scores derived from consurf.

        Parameters
        ----------
        consurf_file : str
            Path to folder containing the files of a consurf analysis results.
        """
        self.consurf_file = consurf_file
        self.data = {}
        cond = False

        # Define score variables to store
        self.data['POS'] = []
        self.data['SEQ'] = []
        self.data['3LATOM'] = []
        self.data['SCORE'] = []
        self.data['COLOR'] = []

        # Read consurf file.
        with open(self.consurf_file, 'r') as cf:
            lines = cf.readlines()
            for l in lines[1:]:
                if '(normalized)' in l:
                    cond = True
                elif l.split() == []:
                    cond = False
                elif cond == True:
                    self.data['POS'].append(int(l.split()[0]))
                    self.data['SEQ'].append(l.split()[1])
                    self.data['3LATOM'].append(l.split()[2])
                    self.data['SCORE'].append(float(l.split()[3]))
                    self.data['COLOR'].append(int(l.split()[4].replace('*','')))

        # Convert all data list into numpy arrays
        for d in self.data:
            self.data[d] = np.array(self.data[d])
