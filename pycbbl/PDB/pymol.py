import __main__
__main__.pymol_argv = ['pymol','-qc']
import pymol
from pymol import cmd, stored
import sys
stdout = sys.stdout
stderr = sys.stderr
pymol.finish_launching(['pymol', '-qc'])
sys.stdout = stdout
sys.stderr = stderr

import os

def alignPDBs(folder_name, pdb_list=None, alignment_index=0, align=True, verbose=True):
    """
    Performs alpha-carbon based alignment of PDB models inside a specific folder.
    The PDBs to be aligned can be narrowed down by giving a list of PDB file names.
    An index can be used to use a specific PDB in the list as the reference for
    the full alignment. The method returns a dictionary containing the RMSD values
    and the number of atoms used into a specific alignment.

    The method can be used only to calculate the RMSD values without saving the
    aligned models by giving the option align=False.

    Parameters
    ----------
    folder_name : str
        Path to the folder where PDB models are found.
    pdb_list : list
        Optional list with the PDBs to consider in the alignment.
    alignment_index : int
        Index of the file in the list to be used as the reference in the alignment.
    align : bool
        Whether to align or not the PDBs in the folder (set to False to only calculate
        RMSDs).
    verbose : bool
        Verbose mode

    Return
    ------

    alignment_details : dict
        Dictionary containing the RMSD values and the number of alpha carbon atoms
        used in the alignment.
    """

    cwd = os.getcwd()
    #Align all structures inside this folder
    os.chdir(folder_name)

    #Read pdbs in folder or read given pdbs
    if pdb_list == None:
        pdbs = []
        for f in os.listdir('.'):
            if f.endswith('.pdb'):
                pdbs.append(f)
        if pdbs == []:
            os.chdir(cwd)
            raise ValueError('There is no files with .pdb extension in the input folder. Please check that your input folder is correct.')
    else:
        pdbs = pdb_list

    # Load pdbs and align them
    for i in range(len(pdbs)):
        if i == alignment_index:
            cmd.load(pdbs[i], 'reference')
        else:
            cmd.load(pdbs[i], pdbs[i])
    if verbose:
        print('Aligning to model: '+pdbs[alignment_index])
    sys.stdout = open('pymolout.tmp', 'w')
    cmd.extra_fit( 'name CA', 'reference', 'super', object='aln_super')
    sys.stdout = stdout

    alignment_details = {}
    with open('pymolout.tmp') as pof:
        for l in pof:
            f = l.split()[0]
            rms = float(l.split()[3])
            atoms = int(l.split()[4].replace('(',''))
            alignment_details[f] = {}
            alignment_details[f]['RMS'] = rms
            alignment_details[f]['atoms'] = atoms
    if align:
        for i in range(len(pdbs)):
            if i != alignment_index:
                cmd.save(pdbs[i], pdbs[i])
    cmd.delete('all')
    os.remove('pymolout.tmp')
    os.chdir(cwd)

    return alignment_details
