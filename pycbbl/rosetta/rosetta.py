import os
import shutil
import pandas as pd
import numpy as np
from pycbbl.PDB import PDBsToTrajectory

def readScoreFile(score_file, skip_failed=False):
    """
    Function to read score files from Rosetta outputs as a panda DataFrame. The
    structure tags are the indexes of the dataframe.

    Parameters
    ----------
    score_file : string
        Path to the rosetta score file.
    skip_failed : bool
        Whether to skip incomplete score lines.
    Returns
    -------
    scores : pandas.DataFrame
    """

    with open(score_file) as sf:
        # Define lines in file
        lines = sf.readlines()

        # Gather score term names
        score_terms = None
        score_count = None
        for c, line in enumerate(lines):
            if line.startswith('SCORE:'):
                score_terms = line.split()
                score_count = c
                break

        if score_count == None:
            raise InputError('Score file is empty!')
        # Remove score terms names from lines
        lines.pop(score_count)

        # Store score values in dictionary
        scores = {}
        for c, line in enumerate(lines):
            if line.startswith('SCORE:'):
                if len(line.split()) != len(score_terms):
                    if skip_failed:
                        print('Warning! Skipped line %s in file %s' % (c+2, score_file))
                        continue
                    else:
                        raise Warning('File %s has incorrect format. Please review the\
 number of items in line %s.' % (score_file, c+2))
                for i, score in enumerate(score_terms):
                    if score not in scores:
                        scores[score] = []
                    try:
                        scores[score].append(float(line.split()[i]))
                    except:
                        scores[score].append(line.split()[i])

    scores.pop('SCORE:')
    for s in scores:
        scores[s] = np.array(scores[s])

    scores = pd.DataFrame(scores)
    scores = scores.set_index('description')
    scores = scores.sort_index()

    return scores

def extractPDBsFromSilent(silent_file, output_folder):
    """
    Extract pdb files from a silent file into a specified folder.

    Parameters
    ----------
    silent_file : str
        Path to silent file
    output_folder : str
        Path to output_folder

    Returns
    -------
    pdbs : list
        Return a list of the PDBs inside the output folder. Independently if they
        were extracted from the silent file.
    """
    # Get current working directory
    cwd = os.getcwd()

    if output_folder.endswith('/'):
        back_path = '../'*(len(output_folder.split('/'))-1)
    else:
        back_path = '../'*len(output_folder.split('/'))

    os.chdir(output_folder)
    os.system('extract_pdbs.linuxgccrelease -silent '+back_path+silent_file)
    os.chdir(back_path)

    # return PDBs inside the folder
    pdbs = sorted([x for x in os.listdir(output_folder)])
    return pdbs

def silentToDCDTrajectory(silent_file, superpose=True, tmp_dir='TMP_SILENT_PDBs',
                          return_tags=False, verbose=True):
    """
    Convert a silent file into a dcd trajectory. The silent file must contain PDBs
    with the same topology (not design mode). PDBs are written momentarily into a
    temporary folder (tmp_dir) and then converted into a dcd file. They are sorted
    by filename before reading them.The temporary folder is deleted at the end of
    the conversion.

    Parameters
    ----------
    silent_file : str
        Path to silent file.
    superpose : bool
        Whether to align the trajectory before returning it.
    return_tags : bool
        If True the function will also return a list containg the tags of the PDBs.

    Returns
    -------
    trajectory : md.Trajectory
        A MD traj trajectory object containing the coordinates of the loaded PDBs

    if return_tags == True
        (trajectory, tags): tuple
            Tuple containing the trajectory and a list with the PDB tags.
    """

    if not os.path.exists(tmp_dir):
        os.mkdir(tmp_dir)
    if verbose:
        print('Extracting PDBs from silent file...')
    extractPDBsFromSilent(silent_file, tmp_dir)
    if verbose:
        print('Done.')
    if verbose:
        print('Converting PDBs into a MDtraj trajectory object...')
    if return_tags:
        trajectory, tags = PDBsToTrajectory(tmp_dir, return_filenames=True)
    else:
        trajectory = PDBsToTrajectory(tmp_dir)
    if verbose:
        print('Done.')
    trajectory.superpose(trajectory[0])
    shutil.rmtree(tmp_dir)

    if return_tags:
        return (trajectory, tags)
    else:
        return trajectory

class scoring:

    """
    Class to hold methods for scoring analysis
    """

    def boltzmannAveragedScore(total_score, score, KT=1):
        """
        Calculate the Boltzmann averaged value of a particular score. Boltzmann
        probabilities are calculated for each member of the esemmble of poses based
        on their total scores. These probabilities are used to calculated the expectation
        value of a particular "score" of interest by summing all the score terms
        pondered by their respective Boltzmann probability.

        Parameters
        ----------
        total_score : np.ndarray
            Array containing the total score values
        score : np.ndarray
            Score values of interest.
        KT : float
            The energy partition constant

        Returns
        -------
        ba_score : float
            The Boltzmann-averaged score.
        """

        relative_energy = total_score - np.amin(total_score)    #Relative total_score to the minimum energy pose
        partition_coefficient = np.sum(np.exp(-relative_energy/KT))
        boltzmann_probabilities = np.exp(-relative_energy/KT)/partition_coefficient
        ba_score = np.sum([score[i]*boltzmann_probabilities[i] for i in range(len(score))])
        return ba_score
