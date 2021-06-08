distance_matrix = pairwise_distances(blobs)
clusterer = hdbscan.HDBSCAN(metric='precomputed',min_cluster_size=3)
clusterer.fit(distance_matrix)
clusterer.labels_

clusterer = hdbscan.HDBSCAN(min_cluster_size=10)
cluster_labels = clusterer.fit_predict(data)

import math
import numpy as np
from event_model import event_model as em

class HDBCluster:

    def __init__(self, event):
        self.hits = sort_by_phi(event.hits) # sort hits by phi
        # self.hits = sort_by_phi_projected(event.hits) # sort by projected phi into the plane
    
    # method that checks previous tracks
    def solve(self):
        distance_matrix = computeDistances(self.hits)
        clusterer = hdbscan.HDBSCAN(metric='precomputed',min_cluster_size=3)
        clusterer.fit(distance_matrix)
        clusterer.labels_

def computeDistances(hits):
    distanceMatrix = np.array((len(hits),len(hits)))
    for i in range(len(hits)):
        for j in range(i,len(hits)):

def distance(hit1,hit2):
    if(hit1.module_number ==  hit2.module_number):
        return 100000
    elif((hit1.module_number % 2) == (hit2.module_number%2)):
        return (abs(math.atan2(hit1.y, hit1.x) - math.atan2(hit2.y, hit2.x)/abs(hit1.module_number - hit2.module_number))
    else:
        return abs(math.atan2(hit1.y, hit1.x) - math.atan2(hit2.y, hit2.x)

def sort_by_phi(hits):
    phis = []
    for h in hits:
        # h.update_polar()
        phis.append(math.atan2(h.y, h.x))
    sorted_index = np.argsort(phis)
    x = np.array(hits)[sorted_index]
    return x