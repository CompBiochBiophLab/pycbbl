import os
from Bio import AlignIO
import shutil

class mafft:
    """
    Class to hold methods to work with mafft executable.

    Methods
    -------
    multipleSequenceAlignment()
        Execute a multiple sequence alignment of the input sequences
    """

    def multipleSequenceAlignment(sequences, output=None, anysymbol=False):
        """
        Use the mafft executable to perform a multiple sequence alignment.

        Parameters
        ----------
        sequences : dict
            Dictionary containing as values the strings representing the sequences
            of the proteins to align and their identifiers as keys.
        output : str
            File name to write the fasta formatted alignment output.
        anysymbol : bool
            Use unusual symbols in the alignment

        Returns
        -------
        alignment : Bio.AlignIO
            Multiple sequence alignment in Biopython format.
        """

        # Write input file containing the sequences
        with open('sequences.fasta.tmp', 'w') as iff:
            for name in sequences:
                iff.write('>'+name+'\n')
                iff.write(sequences[name]+'\n')

        # Calculate alignment
        command = 'mafft --auto'
        if anysymbol:
            command += ' --anysymbol'
        command += ' sequences.fasta.tmp > sequences.aligned.fasta.tmp'
        os.system(command)

        # Read aligned file
        alignment = AlignIO.read("sequences.aligned.fasta.tmp", "fasta")

        # Remove temporary file
        os.remove('sequences.fasta.tmp')
        if output != None:
            shutil.copyfile('sequences.aligned.fasta.tmp', output)
        os.remove('sequences.aligned.fasta.tmp')

        return alignment

    def readSequenceFastaFile(fasta_file):
        """
        Function to read the sequences in a fasta file into a dictionary.

        Parameters
        ----------
        fasta_file : str
            Path to the input fasta file

        Returns
        -------

        sequences : dict
            Dictionary containing as values the sequences and their identifiers
            as keys.
        """

        sequences = {}
        sequence = ''
        with open(fasta_file) as ff:
            for l in ff:
                if l.startswith('>'):
                    if sequence != '':
                        sequences[identifier] = sequence
                    identifier = l.strip().replace('>','')
                    sequence = ''
                else:
                    sequence += l.strip()
            if sequence != '':
                sequences[identifier] = sequence

        return sequences
