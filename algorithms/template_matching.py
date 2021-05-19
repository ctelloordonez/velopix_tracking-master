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
        numberOfTemplates = 900 # 900 does best on small data set
        templateArray = np.zeros((numberOfTemplates,52))
        '*****************************'
        for h in self.hits:
           polarAngle = math.atan2(h.y, h.x)
           indexInArray = math.floor(((polarAngle+math.pi)/(2*math.pi/numberOfTemplates)))
           # were we want to store the hit is decided by 2pi/numberOfTemplates(i.e. how big do we want the angle difference to be) 
           # and then checking to see how often this angle difference fits into the polar angle + pi
           templateArray[indexInArray][h.module_number-1] = h.id
        
        tracks = checkTemplate(templateArray,numberOfTemplates,self.hits)
        # print(len(tracks))
        return tracks

# This function might be reworked into a template checking method
def checkTemplate(templateArray,numberOfTemplates,hits):
    tracks = []
    for t in range(numberOfTemplates):
        
        skipped = 0
        countEven = 0
        countOdd = 0

        hitIDsEven = []
        hitIDsOdd = []

        temporaryListEven = []
        temporaryListOdd = []

        for k in range(52):
    
            if(k %2 ==0):
                if(templateArray[t,k]!= 0 ):
                    countEven += 1
                    hitIDsEven.append(templateArray[t,k])
                
                elif(countEven > 0 and templateArray[t,k]== 0 and skipped == 0):
                    skipped +=1
                    
                elif( templateArray[t,k]== 0 and skipped > 0):
                    if(countEven<3):
                        countEven = 0
                        skipped = 0
                    else:
                        countEven = 0
                        skipped = 0
                        newTrack=convertToTrack2(hitIDsEven,hits)
                        tracks.append(newTrack)
                        hitIDsEven = []

            else:
                if(templateArray[t,k]!= 0):
                    countOdd += 1
                    hitIDsOdd.append(templateArray[t,k])
    
                elif(countOdd > 0 and templateArray[t,k]== 0 and skipped == 0):
                    skipped +=1
                elif(countOdd> 0 and templateArray[t,k] == 0 and skipped > 0):
                    if(countOdd<3):
                        countOdd = 0
                        skipped = 0
                    else:
                        countOdd = 0
                        skipped = 0
                        newTrack = convertToTrack2(hitIDsOdd,hits)
                        tracks.append(newTrack)
                        hitIDsOdd = []
                   
    
        if(countEven > 2):
            newTrack=convertToTrack2(hitIDsEven,hits)
            tracks.append(newTrack)
        if(countOdd > 2):
            newTrack = convertToTrack2(hitIDsOdd,hits)
            tracks.append(newTrack)
        
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