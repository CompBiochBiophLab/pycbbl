import numpy as np
from scipy.spatial.distance import squareform
from scipy.cluster.hierarchy import fcluster, dendrogram, linkage
from matplotlib import pyplot as plt

class clusterByDistance:
    """
    Class to hold methods for generating clusters based on a distance matrix given as input.

    The method uses the cluster hierarchical package from scipy and methods within.
    """

    def __init__(self, distance_matrix, verbose=False):
        """
        At initialization the method calls the linkage method to generate an initialization
        set of clusters.
        """

        self.distances = squareform(distance_matrix, checks=False)
        self.verbose = verbose

        #Calculate clusters
        self.distance_clusters = linkage(reduced_distances, method='average')

        if verbose:
            # Ignore diagonal values for max and min printing
            mask = np.ones(D.shape, dtype=bool)
            np.fill_diagonal(mask, 0)
            max_value = D[mask].max()
            min_value = D[mask].min()

            print('Max pairwise distance: %.3f [diagonal values ignored]' % max_value)
            print('Min pairwise distance: %.3f [diagonal values ignored]' % min_value)
