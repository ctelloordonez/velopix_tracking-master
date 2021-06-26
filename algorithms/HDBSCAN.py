
import math
import numpy as np
from event_model import event_model as em
import sklearn

from sklearn import cluster

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
            distanceMatrix[i][j] = distance7(hits[i],hits[j])
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
      

def distance3(hit1,hit2):
    if(hit1.module_number ==  hit2.module_number):
        return 10
    elif((hit1.module_number % 2) == (hit2.module_number % 2) and hit1.module_number !=  hit2.module_number):
        return (abs(math.atan2(hit1.y, hit1.x) - math.atan2(hit2.y, hit2.x)))/math.sqrt(hit1.z**2 + hit2.z**2)
    else:
        return  1
 

def distance4(hit1,hit2):
    if(hit1.module_number ==  hit2.module_number):
        return 100000
    elif((hit1.module_number % 2) == (hit2.module_number % 2) and hit1.module_number !=  hit2.module_number):
        return 10000*(abs(math.atan2(hit1.y, hit1.x) - math.atan2(hit2.y, hit2.x)))+ abs(hit1.module_number - hit2.module_number)
    else:
        return 50000*(abs(math.atan2(hit1.y, hit1.x) - math.atan2(hit2.y, hit2.x)))+ abs(hit1.module_number - hit2.module_number)

def distance5(hit1,hit2):
    if(math.atan2(hit2.y, hit2.x) - math.atan2(hit1.y, hit1.x)<0.01):

        if((hit1.module_number % 2) == (hit2.module_number % 2) and hit1.module_number !=  hit2.module_number and abs(hit1.module_number - hit2.module_number) < 10):
            return abs(hit1.module_number - hit2.module_number)
        else:
            return 500
    else:
        return 500

def distance6(hit1,hit2):
    if(hit1.module_number ==  hit2.module_number):
        return 100000
    elif((hit1.module_number % 2) == (hit2.module_number % 2) and hit1.module_number !=  hit2.module_number):
        return (abs(math.atan2(hit1.y, hit1.x) - math.atan2(hit2.y, hit2.x)))/(2*math.pi/20000) + abs(hit1.module_number - hit2.module_number)
    else:
        return 10000

def distance7(hit1,hit2):
    if(hit1.module_number ==  hit2.module_number):
        return 100000
    elif((hit1.module_number % 2) == (hit2.module_number % 2) and hit1.module_number !=  hit2.module_number):
        if(hit1.module_number >  hit2.module_number):
            return 500*(abs(math.atan2(hit1.y, hit1.x) - math.atan2(hit2.y, hit2.x)))/(polarDistance(hit1)-polarDistance(hit2)+1e-7)
        else:
            return 500*(abs(math.atan2(hit1.y, hit1.x) - math.atan2(hit2.y, hit2.x)))/(polarDistance(hit2)-polarDistance(hit1)+1e-7)
    else:
        return 10000

def distance8(hit1,hit2):
    if(hit1.module_number ==  hit2.module_number):
        return 100000
    elif((hit1.module_number % 2) == (hit2.module_number % 2) and hit1.module_number !=  hit2.module_number):
        if(hit1.module_number >  hit2.module_number):
            return 500*(abs(math.atan2(hit1.y, hit1.x) - math.atan2(hit2.y, hit2.x)))/(polarDistance(hit1)-polarDistance(hit2)+1e-7)
        else:
            return  500*(abs(math.atan2(hit1.y, hit1.x) - math.atan2(hit2.y, hit2.x)))/(polarDistance(hit2)-polarDistance(hit1)+1e-7)
    else:
        return 10000
def polarDistance(hit1):
    return math.sqrt(hit1.x**2 + hit1.y**2)

def sort_by_phi(hits):
    phis = []
    for h in hits:
        # h.update_polar()
        phis.append(math.atan2(h.y, h.x))
    sorted_index = np.argsort(phis)
    x = np.array(hits)[sorted_index]
    return x