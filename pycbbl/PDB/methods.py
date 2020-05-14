from Bio import PDB
import os
import shutil
import numpy as np
from collections import OrderedDict

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

    Parameters
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

def chainAsStructure(chains):
    """
    This method creates a new Structure object containing only the given chains.

    Parameters
    ----------
    chains : list or Bio.PDB.Chain
        Chain or chains to be added to the new structure object.

    Returns
    -------
    structure : Bio.PDB.Structure
    """

    if not isinstance(chains, list):
        chains = [chains]

    structure = PDB.Structure.Structure(0)
    model = PDB.Model.Model(0)
    for chain in chains:
        model.add(chain)
    structure.add(model)

    return structure

def saveStructureToPDB(structure, output):
    """
    Saves a structure into a PDB file

    Parameters
    ----------
    structure : list or Bio.PDB.Structure
        Structure to save
    """

    io = PDB.PDBIO()
    io.set_structure(structure)
    io.save(output)

class blast:
    """
    Class to hold methods to work with blast executable.

    Methods
    -------

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
