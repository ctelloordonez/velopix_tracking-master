from event_model import event_model as em
from algorithms.AffinityPropgation import computeDistances
from sklearn.cluster import AgglomerativeClustering


class AgglomerativeClusteringWrapper:
    def __init__(self, event):
        self.hits = event.hits

    def solve(self):
        distance_matrix = computeDistances(self.hits)
        clusterer = AgglomerativeClustering(n_clusters=None, affinity='precomputed', linkage='complete',
                                            distance_threshold=0.01)
        clusterer.fit_predict(distance_matrix)
        clusterer.labels_.max()

        tracks = []
        for t in range(clusterer.labels_.max()+1):
            currentTrack = []
            for k in range(len(self.hits)):
                if clusterer.labels_[k] == t:
                    currentTrack.append(self.hits[k])
            tracks.append(em.track(currentTrack))

        return tracks
