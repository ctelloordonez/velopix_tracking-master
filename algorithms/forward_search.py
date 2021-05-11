import math
import numpy as np
from event_model import event_model as em

class ForwardSearch:

    def __init__(self, event):
        self.hits = sort_by_phi(event.hits) # sort hits by phi
        # self.hits = sort_by_phi_projected(event.hits) # sort by projected phi into the plane
    
    # method that checks previous tracks
    def solve(self):

        '***********Parameters ***************'
        NumberOfPreviousTracks = 200     # Parameter for number of previous tracks to check
        XminY = 0.01 # accepted deviation between the x and y ratio values of the track and a hit h
        YminZ = 0.01 # accepted deviation between the z and y ratio values of the track and a hit h
        XminZ = 0.01 # accepted deviation between the x and z ratio values of the track and a hit h
        angleDifference = 0.005 # accepted deviation between the polar angle in x,y values of the current track and a hit h
        trackGreaterThan = 1 # we only make the current_track a real track if its length is greater than this value
        beginOrEnd = -1 # set to 0 to compare hit polar angle to polar angle of the first hit of the track, or -1 to compare to last hit of the track
        '*************************************'
        usedHits = np.zeros(len(self.hits))
        tracks = [] # list of tracks found
        current_track = [] # list of hits representing the track currently under consideration
        current_track.append(self.hits[0]) # initialize by considering the first hit
        self.hits = self.hits[1:] # Take the first hit out of consideration for forming a track with itself
        skipped = 0 # variable to keep track of the number of modules skipped
        
        for i in range(len(self.hits)-1) : # loop over the hits in the event
            if(self.hits[i].module_number != self.hits[i+1].module_number and checkEven(self.hits[i]) == checkEven(self.hits[i+1]) ):
               current_track.append(self.hits[i+1])
               current_track.append(self.hits[i])

               
        return tracks

def sort_by_phi(hits):
    phis = []
    for h in hits:
        phis.append(math.atan2(h.y, h.x))
    sorted_index = np.argsort(phis)
    x = np.array(hits)[sorted_index]
    return x

def sort_by_phi_projected(hits):
    phis = []
    for h in hits:
        phis.append(math.atan2(h.y/abs(h.z-1), h.x/abs(h.z-1)))
    sorted_index = np.argsort(phis)
    x = np.array(hits)[sorted_index]
    return x

# function that returns module number, needed for sort_by_moduleNumber
def moduleNumber(hit):
    return hit.module_number

# this sorts a set of hits by module number in ascending order
def sort_by_moduleNumber(hits):
    return hits.sort(key = moduleNumber)

def calculateMonotoneApproximation(hit1,hit2):
    if(hit1.module_number > hit2.module_number):
        x = (polarDistance(hit1) - polarDistance(hit2))/(hit1.z-hit2.z)
        return x
    if(hit1.module_number < hit2.module_number):
        x = (polarDistance(hit2) - polarDistance(hit1))/(hit2.z-hit1.z)
        return x
    else:
        return 1000000.0

# Function that checks whether 3 hits are monotone w.r.t polar distance r as a function of module number
# Put in in order of module number!
def checkMonotonicity(hit1,hit2,hit3):
    if(hit1.pol_r < hit2.pol_r):
        if(hit2.pol_r <= hit3.pol_r):
            return True
        else:
            return False

    if(hit1.pol_r > hit2.pol_r):
        if(hit2.pol_r >= hit3.pol_r):
            return True
        else:
            return False
    else:
        return True

# check whether a hit is in even or odd module 
def checkEven(hit):
    if (hit.module_number % 2 == 0):
        return True
    else:
        return False

# calculate the direction vector of line by subtracting x,y and z for two hits assumed to be on this line
def calculateDirectionVector(hit1,hit2):
    if(hit1.module_number < hit2.module_number):
        xDirection = hit1.x - hit2.x
        yDirection = hit1.y - hit2.y
        zDirection = hit1.z - hit2.z
    else:
        xDirection = hit2.x - hit1.x
        yDirection = hit2.y - hit1.y
        zDirection = hit2.z - hit1.z
    
    if(xDirection ==0):
        xDirection = hit1.x_
    elif(yDirection == 0):
        yDirection = hit1.y
    elif(zDirection ==0):
        zDirection == hit1.z
    
    
    return xDirection,yDirection,zDirection

def polarDistance(hit1):
    return math.sqrt(hit1.x**2 + hit1.y**2)

def removeTracks(tracks, minLength):
    i = 0
    length = len(tracks)
    while(i < length):
        if(len(tracks[i].hits) < minLength):
            del tracks[i]
            length = length-1
        else:
            i=i+1
    return tracks
    
def matchTracks(tracks,numberOfTracks,MonotoneDifference):
    for i in range(len(tracks)):
        for j in range(numberOfTracks):
             if(abs(calculateMonotoneApproximation(tracks[i].hits[0],tracks[i+1+j].hits[0]) - calculateMonotoneApproximation(tracks[i].hits[-1],tracks[i+1+j].hits[-1])) < MonotoneDifference ):
                for z in range(len(tracks[i+1+j].hits)):
                    tracks[i].hits.append(tracks[i+1+j].hits[0]) # add it at the end