import math
import numpy as np
from event_model import event_model as em

class TemplateMatching:

    def __init__(self, event):
        self.hits = event.hits # sort hits by phi
        # self.hits = sort_by_phi_projected(event.hits) # sort by projected phi into the plane
    
    # method that checks previous tracks
    def solve(self):
        '********parameters***********'
        # numberOfTemplates =  the number of parts we split the polar angle into
        numberOfTemplates = 900 # 900 does best on small data set
        templateArray = np.zeros((numberOfTemplates,52)) # Array to store hit indexes in for template matching
        '*****************************'
        # In this loop, we give each hit its place in the template matching array
        for i in range(len(self.hits)):
           polarAngle = math.atan2(self.hits[i].y, self.hits[i].x)
           indexInArray = math.floor(((polarAngle+math.pi)/(2*math.pi/numberOfTemplates)))
           # were we want to store the hit is decided by 2pi/numberOfTemplates(i.e. how big do we want the angle difference to be) 
           # and then checking to see how often this angle difference fits into the polar angle + pi
           
           if(indexInArray > numberOfTemplates-1): # floats can go just over the last index in the array
               indexInArray = numberOfTemplates-1

            # we add the index +1 in the array at the spot consistent with its polar angle an module number
           templateArray[indexInArray][self.hits[i].module_number-1] = i+1 #h.id +1 because we start with an array of zero's
        
        # Methods that thakes the template array, number of templates and the hits and returns the tracks
        tracks = checkTemplate(templateArray,numberOfTemplates,self.hits)
    
        return tracks

# This function goes over the template array to find tracks
def checkTemplate(templateArray,numberOfTemplates,hits):
    tracks = [] # initialize list to store tracks in.

    # numberOfTemplates =  the number of parts we split the polar angle into.
    # Loop over the different angles
    for t in range(numberOfTemplates):
        
        # Counters for number of hits in sequence and skipped modules for both odd and even modules
        skippedEven = 0
        skippedOdd = 0
        countEven = 0
        countOdd = 0

        hitIDsEven = [] # list to keep potential track hits in even modules
        hitIDsOdd = [] # list to keep potential track hits in odd modules

        # loop over the module numbers
        for k in range(52):
    
            if(k %2 ==0): # checking even
                if(templateArray[t,k]!= 0): # this means there is a hit in it
                    countEven += 1
                    index = int(templateArray[t,k]-1)
                    hitIDsEven.append(hits[index]) # append to even
                
                elif(countEven > 0 and templateArray[t,k]== 0 and skippedEven == 0): 
                    skippedEven +=1 # allow to skip this one module
                    
                elif( templateArray[t,k]== 0 and skippedEven > 0): # When we cannot skip anymore
                    if(countEven<3):   # The track is to short
                        countEven = 0
                        skippedEven = 0
                    else: # add track and reset counters
                        countEven = 0
                        skippedEven = 0
                        tracks.append(em.track(hitIDsEven))
                        hitIDsEven = []

            else: # odd
                if(templateArray[t,k]!= 0): # This means there is a hit in this element of the array
                    countOdd += 1 # update count
                    index = int(templateArray[t,k]-1)
                    hitIDsOdd.append(hits[index]) # append hit
    
                elif(countOdd > 0 and templateArray[t,k]== 0 and skippedOdd == 0): # we allow it to skip once
                    skippedOdd +=1
                elif(countOdd> 0 and templateArray[t,k] == 0 and skippedOdd > 0): # We cannot skip anymore
                    if(countOdd<3): # The track is to short
                        countOdd = 0
                        skippedOdd = 0
                    else: # track is 3 or more hits, so add to tracks and reset counters and list
                        countOdd = 0
                        skippedOdd = 0
                        tracks.append(em.track(hitIDsOdd))
                        hitIDsOdd = []
                   
        # If the last track is longer than 2, it can still be added
        if(countEven > 2): 
            tracks.append(em.track(hitIDsEven))
        if(countOdd > 2): 
            tracks.append(em.track(hitIDsOdd))
           
        
    return tracks

def sort_by_phi(hits):
    phis = []
    for h in hits:
        phis.append(math.atan2(h.y, h.x))
    sorted_index = np.argsort(phis)
    x = np.array(hits)[sorted_index]
    return x

def convertToTrack(listOfHitIDs):
    currentTrack = []
    for i in listOfHitIDs:
       tempHit = em.hit( x = 0, y = 0, z = 0, hit_id = i, module=-1, t=0)
       currentTrack.append(tempHit)
    track = em.track(currentTrack)
    return track

def convertToTrack2(listOfHitIDs,hits):
    currentTrack = []
    for h in hits:
        for i in listOfHitIDs:
            if(h.id == i):
                currentTrack.append(h)
    track = em.track(currentTrack)
    return track

# This function might be reworked into a template checking method
def checkTemplate2(templateArray,numberOfTemplates,hits):
    tracks = []
    for t in range(numberOfTemplates):
        
        skippedEven = 0
        skippedOdd = 0
        countEven = 0
        countOdd = 0

        hitIDsEven = []
        hitIDsOdd = []

        for k in range(52):
    
            if(k %2 ==0): # checking even
                if(templateArray[t,k]!= 0 ):
                    countEven += 1
                    index = int(templateArray[t,k]-1)
                    hitIDsEven.append(hits[index])
                
                elif(countEven > 0 and templateArray[t,k]== 0 and skippedEven == 0):
                    skippedEven +=1
                    
                elif( templateArray[t,k]== 0 and skippedEven > 0):
                    if(countEven<3):
                        countEven = 0
                        skippedEven = 0
                    else:
                        countEven = 0
                        skippedEven = 0
                        # newTrack=convertToTrack2(hitIDsEven,hits)
                        tracks.append(em.track(hitIDsEven))
                        hitIDsEven = []

            else:
                if(templateArray[t,k]!= 0):
                    countOdd += 1
                    index = int(templateArray[t,k]-1)
                    hitIDsOdd.append(hits[index])
    
                elif(countOdd > 0 and templateArray[t,k]== 0 and skippedOdd == 0):
                    skippedOdd +=1
                elif(countOdd> 0 and templateArray[t,k] == 0 and skippedOdd > 0):
                    if(countOdd<3):
                        countOdd = 0
                        skippedOdd = 0
                    else:
                        countOdd = 0
                        skippedOdd = 0
                        # newTrack = convertToTrack2(hitIDsOdd,hits)
                        tracks.append(em.track(hitIDsOdd))
                        hitIDsOdd = []
                   
    
        if(countEven > 2):
            # newTrack=convertToTrack2(hitIDsEven,hits)
            tracks.append(em.track(hitIDsEven))
        if(countOdd > 2):
            # newTrack = convertToTrack2(hitIDsOdd,hits)
            tracks.append(em.track(hitIDsOdd))
           
        
    return tracks