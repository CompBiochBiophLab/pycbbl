import os
import pandas as pd
import numpy as np

def readScoreFile(score_file):
    """
    Function to read score files from Rosetta outputs as a panda DataFrame.

    Parameters
    ----------
    score_file : string
        Path to the rosetta score file.

    Returns
    -------
    scores : pandas.DataFrame
    """

    with open(score_file) as sf:
        lines = [x for x in sf.readlines() if x.startswith('SCORE:')]
        score_terms = lines[0].split()
        scores = {}
        for line in lines[1:]:
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

    return scores

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
