import math

import numpy as np
from hdbscan import HDBSCAN
from sklearn.cluster import AgglomerativeClustering, AffinityPropagation, MeanShift
from event_model import event_model as em

# Parameters
distance_threshold = 0.009
z_shift = -0


class Clustering:
    def __init__(self, event, precomputed_affinity=False):
        self.hits_a = [hit_angles(h) for h in event.hits]
        self.hits = event.hits
        self.precomputed_affinity = precomputed_affinity

    def solve_agglomerative_clustering(self):
        if self.precomputed_affinity:
            distance_matrix = compute_distances(self.hits)
            clusterer = AgglomerativeClustering(n_clusters=None, affinity='precomputed', linkage='complete',
                                                distance_threshold=distance_threshold).fit(distance_matrix)
        else:
            clusterer = AgglomerativeClustering(n_clusters=None, affinity='euclidean',
                                                distance_threshold=distance_threshold).fit(self.hits_a)

        return self.tracks_from_clusterer(clusterer)

    def solve_mean_shift(self):
        clusterer = MeanShift(bandwidth=0.009, bin_seeding=False,
                              min_bin_freq=3, cluster_all=False, max_iter=1000).fit(self.hits_a)
        return self.tracks_from_clusterer(clusterer)

    def solve_hdbscan(self):
        if self.precomputed_affinity:
            distance_matrix = compute_distances(self.hits)
            clusterer = HDBSCAN(metric='precomputed', min_cluster_size=3).fit(distance_matrix)

        else:
            clusterer = HDBSCAN(metric='euclidean', min_cluster_size=3).fit(self.hits_a)

        return self.tracks_from_clusterer(clusterer)

    def solve_affinity_propagation(self):
        if self.precomputed_affinity:
            distance_matrix = compute_distances(self.hits)
            clusterer = AffinityPropagation(affinity='precomputed', random_state=None).fit(distance_matrix)

        else:
            clusterer = AffinityPropagation(affinity='euclidean').fit(self.hits_a)

        return self.tracks_from_clusterer(clusterer)

    def tracks_from_clusterer(self, clusterer):
        tracks = []
        for t in range(clusterer.labels_.max() + 1):
            currentTrack = []
            for k in range(len(self.hits)):
                if clusterer.labels_[k] == t:
                    currentTrack.append(self.hits[k])
            tracks.append(em.track(currentTrack))
        return tracks


def hit_angles(hit):
    return [math.atan2(hit.y, hit.x),
            math.atan2(hit.z+z_shift, hit.y),
            math.atan2(hit.x, hit.z+z_shift)]


def compute_distances(hits):
    distanceMatrix = np.empty((len(hits), len(hits)))
    for i in range(len(hits)):
        for j in range(i, len(hits)):
            distanceMatrix[i][j] = distance(hits[i], hits[j])
            distanceMatrix[j][i] = distanceMatrix[i][j]
    return distanceMatrix


def distance(hit1, hit2):
    if hit1.module_number == hit2.module_number:
        return 10
    elif (hit1.module_number % 2) == (hit2.module_number % 2) and hit1.module_number != hit2.module_number:
        return abs(math.atan2(hit1.y, hit1.x) - math.atan2(hit2.y, hit2.x))
    else:
        return 1
        # abs(math.atan2(hit1.y, hit1.x) - math.atan2(hit2.y, hit2.x))*abs(hit1.module_number - hit2.module_number)


def distance2(hit1, hit2):
    if hit1.module_number == hit2.module_number:
        return 10
    elif (hit1.module_number % 2) == (hit2.module_number % 2) and hit1.module_number != hit2.module_number:
        return (abs(math.atan2(hit1.y, hit1.x) - math.atan2(hit2.y, hit2.x))) / math.sqrt(hit1.x ** 2 + hit2.x ** 2)
    else:
        return 1
        # abs(math.ata


def distance3(hit1, hit2):
    if hit1.module_number == hit2.module_number:
        return 10
    elif (hit1.module_number % 2) == (hit2.module_number % 2) and hit1.module_number != hit2.module_number:
        return (abs(math.atan2(hit1.y, hit1.x) - math.atan2(hit2.y, hit2.x))) / math.sqrt(hit1.z ** 2 + hit2.z ** 2)
    else:
        return 1


def distance4(hit1,hit2):
    if hit1.module_number == hit2.module_number:
        return 100000
    elif((hit1.module_number % 2) == (hit2.module_number % 2) and hit1.module_number !=  hit2.module_number):
        return 10000*(abs(math.atan2(hit1.y, hit1.x) - math.atan2(hit2.y, hit2.x)))+ abs(hit1.module_number - hit2.module_number)
    else:
        return 50000*(abs(math.atan2(hit1.y, hit1.x) - math.atan2(hit2.y, hit2.x)))+ abs(hit1.module_number - hit2.module_number)


def distance5(hit1,hit2):
    if math.atan2(hit2.y, hit2.x) - math.atan2(hit1.y, hit1.x)<0.01:

        if((hit1.module_number % 2) == (hit2.module_number % 2) and hit1.module_number !=  hit2.module_number and abs(hit1.module_number - hit2.module_number) < 10):
            return abs(hit1.module_number - hit2.module_number)
        else:
            return 500
    else:
        return 500


def distance6(hit1,hit2):
    if hit1.module_number ==  hit2.module_number:
        return 100000
    elif((hit1.module_number % 2) == (hit2.module_number % 2) and hit1.module_number !=  hit2.module_number):
        return (abs(math.atan2(hit1.y, hit1.x) - math.atan2(hit2.y, hit2.x)))/(2*math.pi/20000) + abs(hit1.module_number - hit2.module_number)
    else:
        return 10000


def distance7(hit1,hit2):
    if hit1.module_number ==  hit2.module_number:
        return 100000
    elif (hit1.module_number % 2) == (hit2.module_number % 2) and hit1.module_number !=  hit2.module_number:
        polar_distance_hit1 = math.sqrt(hit1.x ** 2 + hit1.y ** 2)
        polar_distance_hit2 = math.sqrt(hit2.x ** 2 + hit2.y ** 2)
        if(hit1.module_number >  hit2.module_number):
            return 500*(abs(math.atan2(hit1.y, hit1.x) - math.atan2(hit2.y, hit2.x)))/(polar_distance_hit1-polar_distance_hit2+1e-7)
        else:
            return 500*(abs(math.atan2(hit1.y, hit1.x) - math.atan2(hit2.y, hit2.x)))/(polar_distance_hit2-polar_distance_hit1+1e-7)
    else:
        return 10000


def distance8(hit1,hit2):
    if hit1.module_number ==  hit2.module_number:
        return 100000
    elif (hit1.module_number % 2) == (hit2.module_number % 2) and hit1.module_number !=  hit2.module_number:
        polar_distance_hit1 = math.sqrt(hit1.x ** 2 + hit1.y ** 2)
        polar_distance_hit2 = math.sqrt(hit2.x ** 2 + hit2.y ** 2)
        if(hit1.module_number >  hit2.module_number):
            return 500*(abs(math.atan2(hit1.y, hit1.x) - math.atan2(hit2.y, hit2.x)))/(polar_distance_hit1-polar_distance_hit2+1e-7)
        else:
            return 500*(abs(math.atan2(hit1.y, hit1.x) - math.atan2(hit2.y, hit2.x)))/(polar_distance_hit2-polar_distance_hit1+1e-7)
    else:
        return 10000
