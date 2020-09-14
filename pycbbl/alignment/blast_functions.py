import os
import numpy as np
from collections import OrderedDict

class blast:
    """
    Class to hold methods to work with blast executable.

    Methods
    -------
    calculatePIDs()
        Fast method to calculate the PID of a sequence against many.
    _getPIDsFromBlastpOutput()
        Method to parse the ouput of blast and retrieve the PIDs.
    """

    def calculatePIDs(target_sequence, comparison_sequences):
        """
        Calculate the percentage of identity (PID) of a target sequence against a group
        of other sequences.

        Parameters
        ----------
        target_sequence : str
            Target sequence.
        comparison_sequences : list of str
            List of sequences to compare the target sequence against.

        Returns
        -------
        pids : numpy.array
            Array containg the PIDs of the target_sequence to all the comparison_sequences.
        """

        if isinstance(comparison_sequences, str):
            comparison_sequences = [comparison_sequences]
        elif not isinstance(comparison_sequences, list):
            raise ValueError('Comparison sequences must be given a single string or \
            as a list of strings.')

        # Write sequences into a temporary file
        with open('seq1.fasta.tmp', 'w') as sf1:
            sf1.write('>seq1\n')
            sf1.write(str(target_sequence))
        with open('seq2.fasta.tmp', 'w') as sf2:
            for i,s in enumerate(comparison_sequences):
                sf2.write('>seq'+str(i)+'\n')
                sf2.write(str(s)+'\n')
        # Execute blastp calculation
        try:
            os.system('blastp -query seq1.fasta.tmp -subject seq2.fasta.tmp -out ssblast.out.tmp -max_target_seqs '+str(len(comparison_sequences)))
        except:
            raise ValueError('blastp executable failed!')

        # Parse blastp output to extract PIDs
        pids = blast._getPIDsFromBlastpOutput('ssblast.out.tmp', len(comparison_sequences))
        pids = np.array(list(pids.values()))

        # Remove temporary files
        os.remove('seq1.fasta.tmp')
        os.remove('ssblast.out.tmp')
        os.remove('seq2.fasta.tmp')

        return pids

    def _getPIDsFromBlastpOutput(blastp_outfile, n_sequences):
        """
        Parse output file from a blastp comparison to extract pids

        Parameters
        ----------
        blastp_outfile : str
            Path to the blastp outputfile.
        n_sequences : str
            number of sequences in the comparison file.

        Returns
        -------
        values : OrderedDict
            Dictionary containing the PID values.
        """

        # Create dictionary integer entries
        values = OrderedDict()
        for i in range(n_sequences):
            values[i] = 0

        # Read PID from blastp output file
        with open(blastp_outfile) as bf:
            for l in bf:
                if l.startswith('>'):
                    seq = int(l.split()[1].replace('seq',''))
                elif 'Identities' in l:
                    pid = eval(l.split()[2])
                    values[seq] = pid

        return values

    def blastPDB(target_sequence, path_to_pdb_fasta, max_target_seqs=500):
        """
        Blast a specific sequence against the PDB database. An updated file of the
        PDB sequences in fasta format must be supplied.
        """
        max_target_seqs = str(max_target_seqs)

        # Write target sequence into a temporary file
        with open('seq1.fasta.tmp', 'w') as sf1:
            sf1.write('>seq1\n')
            sf1.write(str(target_sequence))

        # Execute blastp calculation
        try:
            os.system('blastp -query seq1.fasta.tmp -subject '+path_to_pdb_fasta+' -out pdbblast.out.tmp -max_target_seqs '+max_target_seqs)
        except:
            raise ValueError('blastp executable failed!')

        with open('pdbblast.out.tmp') as bo:
            for l in bo:
                print(l)

        # Remove temporary files
        os.remove('seq1.fasta.tmp')
        os.remove('pdbblast.out.tmp')
