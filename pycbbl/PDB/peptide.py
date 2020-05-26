from pyrosetta import *
from Bio import PDB
import os
init()

def buildAlaninePeptide(input_pdb, output_pdb, chain_id='A', beta_carbon=True):
    """
    Function to build a poly alanine peptide using the backbone and, optionally,
    beta carbon coordinates.

    Parameters
    ----------
    input_pdb : str
        Path to the PDB file containing only the backbone coordinates of the peptide.
    chain_id : str
        Optional Chain ID to give to peptide
    beta_carbon : str
        Wther to use beta carbons to build alanine sidechain positions.
    output_pdb : str
        Path to the PDB output file.
    """

    # Read input pdb
    parser = PDB.PDBParser()
    io = PDB.PDBIO()
    input_structure = parser.get_structure(input_pdb, input_pdb)
    tmp_pdb = input_pdb.replace('.pdb','')+'.tmp.pdb'

    # Create new structure, model and chain objects
    structure = PDB.Structure.Structure(0)
    model = PDB.Model.Model(0)
    chain = PDB.Chain.Chain(chain_id)

    # Define atom names to keep
    bb_atoms = ['N','C','CA','O']
    if beta_carbon:
        bb_atoms.append('CB')
    # Create new alanine residues and add them to the chain
    for i,residue in enumerate(input_structure.get_residues()):
        new_residue = PDB.Residue.Residue((' ', i+1, ' '), 'ALA', residue.segid)
        for atom in residue.get_atoms():
            if atom.name in bb_atoms:
                new_residue.add(atom)
        chain.add(new_residue)

    # Build structure object to write new temporary PDB file
    model.add(chain)
    structure.add(model)
    io.set_structure(structure)
    io.save(tmp_pdb)

    # Read temporary PDB with rosetta to add missing atoms
    temp_pose = Pose()
    pose_from_file(temp_pose, tmp_pdb)
    temp_pose.dump_pdb(output_pdb)

    # Remove temporary PDb file
    os.remove(tmp_pdb)
