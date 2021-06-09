
import math
import numpy as np
from event_model import event_model as em
import sklearn

from sklearn import cluster
from sklearn.cluster import AffinityPropagation

class AffinityPropagationWrapper:

    def __init__(self, event):
        self.hits = event.hits 
    
    # method that checks previous tracks
    def solve(self):

        distance_matrix = computeDistances(self.hits)
        clusterer = AffinityPropagation(affinity='precomputed', random_state=None) # affinity='precomputed'
        clusterer.fit_predict(distance_matrix)
        clusterer.labels_
        clusterer.labels_.max()

        tracks = []
        for t in range(clusterer.labels_.max()+1):
            currentTrack = []
            for k in range(len(self.hits)):
                if(clusterer.labels_[k]== t):
                    currentTrack.append(self.hits[k])
            tracks.append(em.track(currentTrack))

        return tracks

    

        

def computeDistances(hits):
    distanceMatrix = np.empty((len(hits),len(hits)))
    for i in range(len(hits)):
        for j in range(i,len(hits)):
            distanceMatrix[i][j] = distance3(hits[i],hits[j])
            distanceMatrix[j][i] = distanceMatrix[i][j]
    return distanceMatrix

def distance(hit1,hit2):
    if(hit1.module_number ==  hit2.module_number):
        return 10
    elif((hit1.module_number % 2) == (hit2.module_number % 2) and hit1.module_number !=  hit2.module_number):
        return ((abs(math.atan2(hit1.y, hit1.x) - math.atan2(hit2.y, hit2.x))))
    else:
        return  1
        # abs(math.atan2(hit1.y, hit1.x) - math.atan2(hit2.y, hit2.x))*abs(hit1.module_number - hit2.module_number)
    
def distance2(hit1,hit2):
    if(hit1.module_number ==  hit2.module_number):
        return 10
    elif((hit1.module_number % 2) == (hit2.module_number % 2) and hit1.module_number !=  hit2.module_number):
        return (abs(math.atan2(hit1.y, hit1.x) - math.atan2(hit2.y, hit2.x)))/math.sqrt(hit1.x**2 + hit2.x**2)
    else:
        return  1
        # abs(math.ata

def distance3(hit1,hit2):
    if(hit1.module_number ==  hit2.module_number):
        return 10
    elif((hit1.module_number % 2) == (hit2.module_number % 2) and hit1.module_number !=  hit2.module_number):
        return (abs(math.atan2(hit1.y, hit1.x) - math.atan2(hit2.y, hit2.x)))/math.sqrt(hit1.z**2 + hit2.z**2)
    else:
        return  1
        # abs(math.ata

def sort_by_phi(hits):
    phis = []
    for h in hits:
        # h.update_polar()
        phis.append(math.atan2(h.y, h.x))
    sorted_index = np.argsort(phis)
    x = np.array(hits)[sorted_index]
    return x