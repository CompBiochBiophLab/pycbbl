#!/usr/bin/env python
# coding: utf-8

import shutil
import re
from collections import OrderedDict

from Bio import PDB

from pyrosetta import *
init('-ignore_zero_occupancy false')

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

class rosetta:

    def readPoseFromPdb(pdb_filename, zero_occupancy=True):
        """
        Create a pose from the desired PDB file

        Args:
            pdb_filename: path to pdb file

        Returns:
            Returns the pose object
        """

        if zero_occupancy == False:
            init('-ignore_zero_occupancy false')

        temp_pose = Pose()
        pose_from_file(temp_pose, pdb_filename)

        return temp_pose

    def getResidues(pose):
        """
        Get all residues from a pose.

        Args:
            pose: pyrosetta pose object; pyrosetta.rosetta.core.pose()

        Returns:
            Returns the list of residues in the pose
        """
        #Build a list of residues
        residues = range( 1 , pose.total_residue()  + 1 )

        return residues

    def getChains(pose):

        chains = []
        residues = rosetta.getResidues(pose)
        for r in residues:
            c = pose.pdb_info().chain(r)
            if c not in chains:
                chains.append(c)

        return chains

    def getChainResidues(pose,chain):
        """
        Get residues from a pose that belongs to a single chain

        Args:
            pose: pyrosetta pose object; pyrosetta.rosetta.core.pose()
            chain: chain or chain_id; type: integer or str, respectively

        Returns:
            Returns list of residues in chain
        """
        #Evaluate if chain is given as string or as an integer.
        if isinstance(chain, int):
            #Get first and last residue of chain
            lrb = pose.conformation().chain_begin(chain)
            lre = pose.conformation().chain_end(chain)

        elif isinstance(chain, str):
            #Get first and last residue of chain
            chain = pyrosetta.rosetta.core.pose.get_chain_id_from_chain(chain, pose)
            lrb = pose.conformation().chain_begin(chain)
            lre = pose.conformation().chain_end(chain)

        #Return list of residues as list
        residues = range(lrb, lre+1)

        #If no residues has been selected return None
        if residues == []:
            residues = None

        return residues

class removeHydrogensAndDisordered(PDB.Select):

    def accept_atom(self, atom):

        _hydrogen = re.compile("[123 ]*H.*")

        if atom.get_altloc() == 'A':
            atom.set_altloc(' ')

        if _hydrogen.match(atom.get_id()):
            return 0
        elif atom.get_altloc() == 'B':
            return 0
        else:
            return 1

class smog2:
    """
    Class for autamating the generation of SBM contact maps using the SMOG2 executable.

    The Executable can be obtained at:

    http://smog-server.org/smog2/

    The program needs to be installed and the executables avaiable through the linux PATH variable.

    Methods
    -------
    generateSBMTopology(pdb_file, kwarg**)
        Generates an all-heavy atom (AA) or alpha-carbon (CA) SBM topology and a protein contact file
        using SMOG2.
    """

    def correctPDBformat(pdb_file, output_file):

        pose = rosetta.readPoseFromPdb(pdb_file)
        chains = {c:i for i,c in enumerate(rosetta.getChains(pose))}
        pose.dump_pdb(output_file)

        parser = PDB.PDBParser()
        structure = parser.get_structure(pdb_file.replace('.pdb',''), output_file)
        io = PDB.PDBIO()
        io.set_structure(structure)
        io.save(output_file, removeHydrogensAndDisordered())

        with open(output_file, 'r') as pdbf:
            lines = pdbf.readlines()
            if lines[-2].startswith('TER'):
                del lines[-2]
        with open(output_file, 'w') as pdbf:
            for line in lines:
                if line.startswith('ATOM') or line.startswith('TER'):
                    pdbf.write(line)
            pdbf.write('END')

        return output_file


    def generateSBMTopology(pdb_file, model_type='CA', contact_file=None, adjust_pdb=False, output_path=None, change_terminal_names=False):
        """
        Generates an AA or CA SBM topology and a protein contact file using SMOG2.

        Parameters
        ----------
        pdb_file : string
            Path to the input PDB file.
        model_type : string ('CA')
            Whether to generate an all atom 'AA' or a alpha-carbon 'CA' SBM model.
        contact_file : string
            Path to the output contact file.
        adjust_pdb : boolean (False)
            Wheter to adjust the format of the pdb file to be read by smog2. If given True this will
            overwrite the input file.
        output_path : str
            Path to save the adjusted file
        Returns
        -------
        None
        """

        # Executes the adjust pdb command and replaces the input pdb file
        if adjust_pdb:
            adjust_command = 'smog_adjustPDB -default -i '+pdb_file
            os.system(adjust_command)
            if output_path == None:
                output_path = pdb_file.split('/')[-1]
            else:
                output_path = output_path+'/'+pdb_file.split('/')[-1]
            shutil.copyfile('adjusted.pdb', output_path)
            os.remove('adjusted.pdb')
            pdb_file = output_path

        # Executes the SBM topology command for smog2
        if model_type == 'AA':
            command = 'smog2 -i '+pdb_file+' -AA -dname '+pdb_file.replace('.pdb','')
        elif model_type == 'CA':
            command = 'smog2 -i '+pdb_file+' -CA -dname '+pdb_file.replace('.pdb','')
        if contact_file != None:
            command.replace('smog2','smog2 -c '+contact_file)
        os.system(command)

        if change_terminal_names:
            with open(pdb_file) as pf:
                lines = pf.readlines()
            with open(pdb_file, 'w') as of:
                for l in lines:
                    if l.startswith('ATOM'):
                        if len(l.split()) < 9:
                            nname = l.split()[3][:3]+' '+l.split()[3][-1]
                            of.write(l.replace(l.split()[3], nname))
                        else:
                            of.write(l)
                    else:
                        of.write(l)

    def renameTerminalResidues(pdb):
        lines = []
        with open(pdb) as sf:
            for l in sf:
                if l.startswith('ATOM'):
                    if  len(l.split()) == 8:
                        l = l[:20]+' '+l[21:]
                lines.append(l)
        with open(pdb, 'w') as sf:
            for l in lines:
                sf.write(l)

