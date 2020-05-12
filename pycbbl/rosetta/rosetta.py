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
