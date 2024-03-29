import math
import numpy as np
from event_model import event_model as em

class TemplateMatching:

    def __init__(self, event):
        self.hits = event.hits # sort hits by phi
        
    def solve(self):
        '********parameters***********'
        # numberOfTemplates =  the number of parts we split the 2 pi range of the polar angle into.
        numberOfTemplates = 700 # 900 does best on small data set

        #For method 1,2,3: use this template array and comment next out
        #templateArray = np.zeros((numberOfTemplates,52,1)) # Array to store hit indexes in for template matching

        #For method 4,5 use this template array and comment the above out
        templateArray = np.zeros((numberOfTemplates,52,1)) # Array to store hit indexes in for template matching


        '*****************************'
        # In this loop, we give each hit its place in the template matching array
        for i in range(len(self.hits)):
           polarAngle = math.atan2(self.hits[i].y, self.hits[i].x)
           indexInArray = math.floor(((polarAngle+math.pi)/(2*math.pi/numberOfTemplates)))
           # were we want to store the hit is decided by 2pi/numberOfTemplates(i.e. how big do we want the angle difference to be) 
           # and then checking to see how often this angle difference fits into the polar angle + pi
           
           if(indexInArray > numberOfTemplates-1): # floats can go just over the last index in the array
               indexInArray = numberOfTemplates-1



            # we add the hit index +1 in the array at the spot consistent with its polar angle and module number
           if templateArray[indexInArray][self.hits[i].module_number-1][0]==0:
            templateArray[indexInArray][self.hits[i].module_number-1][0] = i+1 #h.id +1 because we start with an array of zero's

           #comment the elif below out for method 1,2,3 (this is for filling in the next instances)
           elif templateArray[indexInArray][self.hits[i].module_number-1][1]==0:
            templateArray[indexInArray][
                   self.hits[i].module_number - 1][1] = i + 1  # h.id +1 because we start with an array of zero's

          # else:
          #  templateArray[indexInArray][self.hits[i].module_number - 1][2] = i + 1

        # Methods that takes the template array, number of templates and the hits and returns the tracks
        tracks = checkTemplate3(templateArray,numberOfTemplates,self.hits)
    
        return tracks

# This function goes over the template array to find tracks
def checkTemplate(templateArray,numberOfTemplates,hits):
    tracks = [] # initialize list to store tracks in.

    # numberOfTemplates =  the number of parts we split the 2 pi range of the polar angle into.
    # Loop over the different angle
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


# Template matching with monotonicity check
def checkTemplate2(templateArray,numberOfTemplates,hits):
    tracks = []
    for t in range(numberOfTemplates):
        
        countEven = 0

        
        countOdd = 0

        hitIDsEven = []
        hitIDsOdd = []

        for k in range(52):
    
            if(k %2 ==0): # checking even
                if(templateArray[t,k]!= 0):
                    if(countEven <= 1):
                        index = int(templateArray[t,k]-1) 
                        countEven += 1
                        hitIDsEven.append(hits[index])
                    else:
                        xTemp = (hitIDsEven[0].x-hits[index].x)*(hitIDsEven[0].x-hitIDsEven[1].x) 
                        yTemp = (hitIDsEven[0].y-hits[index].y)*(hitIDsEven[0].y-hitIDsEven[1].y) 
                        if(xTemp > 0 and yTemp > 0 ):
                            countEven += 1
                            index = int(templateArray[t,k]-1)
                            hitIDsEven.append(hits[index])
                        elif(countEven>2):
                            countEven = 0
                            tracks.append(em.track(hitIDsEven))
                            hitIDsEven = []
                        else:
                            countEven = 1
                            hitIDsEven = []
                            index = int(templateArray[t,k]-1)
                            hitIDsEven.append(hits[index])

            
                elif( templateArray[t,k]== 0 and countEven > 0):
                    if(countEven<3):
                        countEven = 0
                        hitIDsEven = []
                    else:
                        countEven = 0
                        tracks.append(em.track(hitIDsEven))
                        hitIDsEven = []
                else:
                    hitIDsEven = []
                


            else:
                if(templateArray[t,k]!= 0):
                    if(countOdd <= 1):
                        countOdd += 1
                        index = int(templateArray[t,k]-1)
                        hitIDsOdd.append(hits[index])
                    else:
                        xTemp = (hitIDsOdd[0].x-hits[index].x)*(hitIDsOdd[0].x-hitIDsOdd[1].x) 
                        yTemp = (hitIDsOdd[0].y-hits[index].y)*(hitIDsOdd[0].y-hitIDsOdd[1].y) 
                        if(xTemp > 0 and yTemp > 0 ):
                            countOdd += 1
                            index = int(templateArray[t,k]-1)
                            hitIDsOdd.append(hits[index])
                        elif(countOdd>2):
                            countOdd = 0
                            tracks.append(em.track(hitIDsEven))
                            hitIDsOdd = []
                        else:
                            countOdd = 1
                            hitIDsOdd = []
                            index = int(templateArray[t,k]-1)
                            hitIDsOdd.append(hits[index])
    
                
                elif(countOdd> 0 and templateArray[t,k] == 0 ):
                    if(countOdd<3):
                        countOdd = 0
                        hitIDsOdd = []
                    else:
                        countOdd = 0
                        tracks.append(em.track(hitIDsOdd))
                        hitIDsOdd = []
                else:
                    hitIDsOdd = []
    
        if(countEven > 2):
            # newTrack=convertToTrack2(hitIDsEven,hits)
            tracks.append(em.track(hitIDsEven))
        if(countOdd > 2):
            # newTrack = convertToTrack2(hitIDsOdd,hits)
            tracks.append(em.track(hitIDsOdd))
        
    return tracks


#This monotonicity check function is used only for method 3
def MonotonicityCheckEven1Instance(countEven,index,hitIDsOdd,hits,templateArray,tracks,hitIDsEven,t,k):
    if (countEven <= 1):
        index = int(templateArray[t, k] - 1)
        countEven += 1
        hitIDsEven.append(hits[index])
        templateArray[t, k] = 0

    else:
        xTemp = (hitIDsEven[0].x - hits[index].x) * (hitIDsEven[0].x - hitIDsEven[1].x)
        yTemp = (hitIDsEven[0].y - hits[index].y) * (hitIDsEven[0].y - hitIDsEven[1].y)
        if (xTemp > 0 and yTemp > 0):
            countEven += 1
            index = int(templateArray[t, k] - 1)
            hitIDsEven.append(hits[index])
            templateArray[t, k] = 0

        elif (countEven > 2):
            countEven = 0
            tracks.append(em.track(hitIDsEven))
            hitIDsEven = []
        else:
            countEven = 1
            hitIDsEven = []
            index = int(templateArray[t, k] - 1)
            hitIDsEven.append(hits[index])
            templateArray[t, k] = 0

    return (countEven, hitIDsOdd,templateArray,tracks,hitIDsEven)
#This monotonicity check function is used only for method 3
def MonotonicityCheckOdd1Instance(countOdd,index,hitIDsOdd,hits,templateArray,tracks,hitIDsEven,t,k):
    if (countOdd <= 1):
        countOdd += 1
        index = int(templateArray[t, k] - 1)
        hitIDsOdd.append(hits[index])
        templateArray[t, k] = 0
    else:
        xTemp = (hitIDsOdd[0].x - hits[index].x) * (hitIDsOdd[0].x - hitIDsOdd[1].x)
        yTemp = (hitIDsOdd[0].y - hits[index].y) * (hitIDsOdd[0].y - hitIDsOdd[1].y)
        if (xTemp > 0 and yTemp > 0):
            countOdd += 1
            index = int(templateArray[t, k] - 1)
            hitIDsOdd.append(hits[index])
            templateArray[t, k] = 0

        elif (countOdd > 2):
            countOdd = 0
            tracks.append(em.track(hitIDsEven))
            hitIDsOdd = []

        else:
            countOdd = 1
            hitIDsOdd = []
            index = int(templateArray[t, k] - 1)
            hitIDsOdd.append(hits[index])
            templateArray[t, k] = 0

    return((countOdd, hitIDsOdd,templateArray,tracks,hitIDsEven))


# Template matching with monotonicity check and multiple angle sections to check
def checkTemplate3(templateArray, numberOfTemplates, hits):
    tracks = []
    for t in range(numberOfTemplates):

        countEven = 0
        countOdd = 0

        hitIDsEven = []
        hitIDsOdd = []

        for k in range(52):

            if (k % 2 == 0):  # checking even
                if (templateArray[t, k] != 0):
                    index = int(templateArray[t, k] - 1)
                    countEven,hitIDsOdd,templateArray,tracks,hitIDsEven=MonotonicityCheckEven1Instance(countEven,index,hitIDsOdd,hits,templateArray,tracks,hitIDsEven,t,k)

                elif (t < numberOfTemplates - 1 and templateArray[t + 1, k] != 0):  # if no hit was found in t,k  look for hit in t+1,k
                    index = int(templateArray[t + 1, k] - 1)
                    countEven, hitIDsOdd,templateArray,tracks,hitIDsEven=MonotonicityCheckEven1Instance(countEven, index,hitIDsOdd,hits,templateArray,tracks,hitIDsEven,t+1,k)

                elif (t > 0 and templateArray[t - 1, k] != 0):  # if no hit was found in t,k or t+1,k, look for hit in t-1,k
                    index = int(templateArray[t - 1, k] - 1)
                    countEven, hitIDsOdd,templateArray,tracks,hitIDsEven=MonotonicityCheckEven1Instance(countEven,index,hitIDsOdd,hits,templateArray,tracks,hitIDsEven,t-1,k)

                elif (templateArray[t, k] == 0 and countEven > 0):
                    if (countEven < 3):
                        countEven = 0
                        hitIDsEven = []
                    else:
                        countEven = 0
                        tracks.append(em.track(hitIDsEven))
                        hitIDsEven = []
                else:
                    hitIDsEven = []

            else:
                if (templateArray[t, k] != 0):
                    index = int(templateArray[t, k] - 1)
                    countOdd, hitIDsOdd,templateArray,tracks,hitIDsEven=MonotonicityCheckOdd1Instance(countOdd, index,hitIDsOdd,hits,templateArray,tracks,hitIDsEven,t,k)

                elif (t < numberOfTemplates - 1 and templateArray[t + 1, k] != 0):  # if no hit was found in t,k  look for hit in t+1,k
                    index = int(templateArray[t + 1, k] - 1)
                    countOdd, hitIDsOdd,templateArray,tracks,hitIDsEven=MonotonicityCheckOdd1Instance(countOdd, index,hitIDsOdd,hits,templateArray,tracks,hitIDsEven,t+1,k)

                elif (t > 0 and templateArray[t - 1, k] != 0):  # if no hit was found in t,k or t+1,k, look for hit in t-1,k
                    index = int(templateArray[t - 1, k] - 1)
                    countOdd, hitIDsOdd,templateArray,tracks,hitIDsEven=MonotonicityCheckOdd1Instance(countOdd, index,hitIDsOdd,hits,templateArray,tracks,hitIDsEven,t-1,k)

                elif (countOdd > 0 and templateArray[t, k] == 0):
                    if (countOdd < 3):
                        countOdd = 0
                        hitIDsOdd = []
                    else:
                        countOdd = 0
                        tracks.append(em.track(hitIDsOdd))
                        hitIDsOdd = []
                else:
                    hitIDsOdd = []

        if (countEven > 2):
            # newTrack=convertToTrack2(hitIDsEven,hits)
            tracks.append(em.track(hitIDsEven))
        if (countOdd > 2):
            # newTrack = convertToTrack2(hitIDsOdd,hits)
            tracks.append(em.track(hitIDsOdd))

    return tracks

#This monotonicity check function is used only for method 4 and 5
def MonotonicityCheckEven(countEven,index,hitIDsOdd,hits,templateArray,tracks,hitIDsEven,t,k):
    if (countEven <= 1):
        countEven += 1
        hitIDsEven.append(hits[index])
        templateArray[t, k][0] = templateArray[t, k][1]
        # templateArray[t, k][1] = templateArray[t, k][2]
        # templateArray[t, k][2] = 0
        templateArray[t, k][1] = 0


    else:
        xTemp = (hitIDsEven[0].x - hits[index].x) * (hitIDsEven[0].x - hitIDsEven[1].x)
        yTemp = (hitIDsEven[0].y - hits[index].y) * (hitIDsEven[0].y - hitIDsEven[1].y)
        if (xTemp > 0 and yTemp > 0):
            countEven += 1
            hitIDsEven.append(hits[index])
            templateArray[t, k][0] = templateArray[t, k][1]
            # templateArray[t, k][1] = templateArray[t, k][2]
            # templateArray[t, k][2] = 0
            templateArray[t, k][1] = 0

        # from here to the next '#' the monotonicity of the next hit in line is checked, this can be left out
        #which for now still gives better performance
        elif templateArray[t,k][1]!=0:
            index = int(templateArray[t, k][1] - 1)
            xTemp = (hitIDsEven[0].x - hits[index].x) * (hitIDsEven[0].x - hitIDsEven[1].x)
            yTemp = (hitIDsEven[0].y - hits[index].y) * (hitIDsEven[0].y - hitIDsEven[1].y)
            if (xTemp > 0 and yTemp > 0):
                countEven += 1
                hitIDsEven.append(hits[index])
                #templateArray[t, k][1] = templateArray[t, k][2]
                #templateArray[t, k][2] = 0
                templateArray[t, k][1] = 0
        # end of monotonicity check for next hit



        elif (countEven > 2):
            countEven = 0
            tracks.append(em.track(hitIDsEven))
            hitIDsEven = []
        else:
            countEven = 1
            hitIDsEven = []
            hitIDsEven.append(hits[index])
            templateArray[t, k][0] = templateArray[t, k][1]
            # templateArray[t, k][1] = templateArray[t, k][2]
            # templateArray[t, k][2] = 0
            templateArray[t, k][1] = 0
    return (countEven, hitIDsOdd,templateArray,tracks,hitIDsEven)


#This monotonicity check function is used only for method 4 and 5
def MonotonicityCheckOdd(countOdd,index,hitIDsOdd,hits,templateArray,tracks,hitIDsEven,t,k):
    if (countOdd <= 1):
        countOdd += 1
        hitIDsOdd.append(hits[index])
        templateArray[t, k][0] = templateArray[t, k][1]
        # templateArray[t, k][1] = templateArray[t, k][2]
        # templateArray[t, k][2] = 0
        templateArray[t, k][1] = 0
    else:
        xTemp = (hitIDsOdd[0].x - hits[index].x) * (hitIDsOdd[0].x - hitIDsOdd[1].x)
        yTemp = (hitIDsOdd[0].y - hits[index].y) * (hitIDsOdd[0].y - hitIDsOdd[1].y)
        if (xTemp > 0 and yTemp > 0):
            countOdd += 1
            hitIDsOdd.append(hits[index])
            templateArray[t, k][0] = templateArray[t, k][1]
            # templateArray[t, k][1] = templateArray[t, k][2]
            # templateArray[t, k][2] = 0
            templateArray[t, k][1] = 0

        # from here to the next '#' the monotonicity of the next hit in line is checked, this can be left out
        #which for now still gives better performance
        elif templateArray[t,k][1]!=0:
            index = int(templateArray[t, k][1] - 1)
            xTemp = (hitIDsOdd[0].x - hits[index].x) * (hitIDsOdd[0].x - hitIDsOdd[1].x)
            yTemp = (hitIDsOdd[0].y - hits[index].y) * (hitIDsOdd[0].y - hitIDsOdd[1].y)
            if (xTemp > 0 and yTemp > 0):
                countOdd += 1
                hitIDsOdd.append(hits[index])
                #templateArray[t, k][1] = templateArray[t, k][2]
                #templateArray[t, k][2] = 0
                templateArray[t, k][1] = 0
        #end point of monotonicity check for second in line

        elif countOdd > 2:
            countOdd = 0
            tracks.append(em.track(hitIDsEven))
            hitIDsOdd = []

        else:
            countOdd = 1
            hitIDsOdd = []
            hitIDsOdd.append(hits[index])
            templateArray[t, k][0] = templateArray[t, k][1]
            # templateArray[t, k][1] = templateArray[t, k][2]
            # templateArray[t, k][2] = 0
            templateArray[t, k][1] = 0


    return((countOdd, hitIDsOdd,templateArray,tracks,hitIDsEven))

# Template matching with monotonicity check and multiple angle sections to check
def checkTemplate4(templateArray, numberOfTemplates, hits):
    tracks = []
    for t in range(numberOfTemplates):

        countEven = 0
        countOdd = 0

        hitIDsEven = []
        hitIDsOdd = []

        for k in range(52):

            if (k % 2 == 0):  # checking even
                if (templateArray[t, k][0] != 0):
                    index = int(templateArray[t, k][0] - 1)
                    countEven,hitIDsOdd,templateArray,tracks,hitIDsEven=MonotonicityCheckEven(countEven,index,hitIDsOdd,hits,templateArray,tracks,hitIDsEven,t,k)

                elif (t < numberOfTemplates - 1 and templateArray[t + 1, k][0] != 0):  # if no hit was found in t,k  look for hit in t+1,k
                    index = int(templateArray[t + 1, k][0] - 1)
                    countEven, hitIDsOdd,templateArray,tracks,hitIDsEven=MonotonicityCheckEven(countEven, index,hitIDsOdd,hits,templateArray,tracks,hitIDsEven,t+1,k)

                elif (t > 0 and templateArray[t - 1, k][0] != 0):  # if no hit was found in t,k or t+1,k, look for hit in t-1,k
                    index = int(templateArray[t - 1, k][0] - 1)
                    countEven, hitIDsOdd,templateArray,tracks,hitIDsEven=MonotonicityCheckEven(countEven,index,hitIDsOdd,hits,templateArray,tracks,hitIDsEven,t-1,k)

                elif (templateArray[t, k][0] == 0 and countEven > 0):
                    if (countEven < 3):
                        countEven = 0
                        hitIDsEven = []
                    else:
                        countEven = 0
                        tracks.append(em.track(hitIDsEven))
                        hitIDsEven = []
                else:
                    hitIDsEven = []

            else:
                if (templateArray[t, k][0] != 0):
                    index = int(templateArray[t, k][0] - 1)
                    countOdd, hitIDsOdd,templateArray,tracks,hitIDsEven=MonotonicityCheckOdd(countOdd, index,hitIDsOdd,hits,templateArray,tracks,hitIDsEven,t,k)

                elif (t < numberOfTemplates - 1 and templateArray[t + 1, k][0] != 0):  # if no hit was found in t,k  look for hit in t+1,k
                    index = int(templateArray[t + 1, k][0] - 1)
                    countOdd, hitIDsOdd,templateArray,tracks,hitIDsEven=MonotonicityCheckOdd(countOdd, index,hitIDsOdd,hits,templateArray,tracks,hitIDsEven,t+1,k)

                elif (t > 0 and templateArray[t - 1, k][0] != 0):  # if no hit was found in t,k or t+1,k, look for hit in t-1,k
                    index = int(templateArray[t - 1, k][0] - 1)
                    countOdd, hitIDsOdd,templateArray,tracks,hitIDsEven=MonotonicityCheckOdd(countOdd, index,hitIDsOdd,hits,templateArray,tracks,hitIDsEven,t-1,k)

                elif (countOdd > 0 and templateArray[t, k][0] == 0):
                    if (countOdd < 3):
                        countOdd = 0
                        hitIDsOdd = []
                    else:
                        countOdd = 0
                        tracks.append(em.track(hitIDsOdd))
                        hitIDsOdd = []
                else:
                    hitIDsOdd = []

        if (countEven > 2):
            # newTrack=convertToTrack2(hitIDsEven,hits)
            tracks.append(em.track(hitIDsEven))
        if (countOdd > 2):
            # newTrack = convertToTrack2(hitIDsOdd,hits)
            tracks.append(em.track(hitIDsOdd))

    return tracks

#this method checks 4 neighbours and 2 instances
def checkTemplate5(templateArray, numberOfTemplates, hits):
    tracks = []
    for t in range(numberOfTemplates):

        countEven = 0
        countOdd = 0

        hitIDsEven = []
        hitIDsOdd = []

        for k in range(52):

            if (k % 2 == 0):  # checking even
                if (templateArray[t, k][0] != 0):
                    index = int(templateArray[t, k][0] - 1)
                    countEven,hitIDsOdd,templateArray,tracks,hitIDsEven=MonotonicityCheckEven(countEven,index,hitIDsOdd,hits,templateArray,tracks,hitIDsEven,t,k)

                elif (t < numberOfTemplates - 1 and templateArray[t + 1, k][0] != 0):  # if no hit was found in t,k  look for hit in t+1,k
                    index = int(templateArray[t + 1, k][0] - 1)
                    countEven, hitIDsOdd,templateArray,tracks,hitIDsEven=MonotonicityCheckEven(countEven, index,hitIDsOdd,hits,templateArray,tracks,hitIDsEven,t+1,k)

                elif (t > 0 and templateArray[t - 1, k][0] != 0):  # if no hit was found in t,k or t+1,k, look for hit in t-1,k
                    index = int(templateArray[t - 1, k][0] - 1)
                    countEven, hitIDsOdd,templateArray,tracks,hitIDsEven=MonotonicityCheckEven(countEven,index,hitIDsOdd,hits,templateArray,tracks,hitIDsEven,t-1,k)

                elif (t < numberOfTemplates - 2 and templateArray[t + 2, k][0] != 0):  # if no hit was found in t,k  look for hit in t+1,k
                    index = int(templateArray[t + 2, k][0] - 1)
                    countEven, hitIDsOdd,templateArray,tracks,hitIDsEven=MonotonicityCheckEven(countEven, index,hitIDsOdd,hits,templateArray,tracks,hitIDsEven,t+1,k)


                elif (t > 1 and templateArray[t - 2, k][0] != 0):  # if no hit was found in t,k or t+1,k, look for hit in t-1,k
                    index = int(templateArray[t - 2, k][0] - 2)
                    countEven, hitIDsOdd,templateArray,tracks,hitIDsEven=MonotonicityCheckEven(countEven,index,hitIDsOdd,hits,templateArray,tracks,hitIDsEven,t-1,k)

                elif (templateArray[t, k][0] == 0 and countEven > 0):
                    if (countEven < 3):
                        countEven = 0
                        hitIDsEven = []
                    else:
                        countEven = 0
                        tracks.append(em.track(hitIDsEven))
                        hitIDsEven = []
                else:
                    hitIDsEven = []

            else:
                if (templateArray[t, k][0] != 0):
                    index = int(templateArray[t, k][0] - 1)
                    countOdd, hitIDsOdd,templateArray,tracks,hitIDsEven=MonotonicityCheckOdd(countOdd, index,hitIDsOdd,hits,templateArray,tracks,hitIDsEven,t,k)

                elif (t < numberOfTemplates - 1 and templateArray[t + 1, k][0] != 0):  # if no hit was found in t,k  look for hit in t+1,k
                    index = int(templateArray[t + 1, k][0] - 1)
                    countOdd, hitIDsOdd,templateArray,tracks,hitIDsEven=MonotonicityCheckOdd(countOdd, index,hitIDsOdd,hits,templateArray,tracks,hitIDsEven,t+1,k)

                elif (t > 0 and templateArray[t - 1, k][0] != 0):  # if no hit was found in t,k or t+1,k, look for hit in t-1,k
                    index = int(templateArray[t - 1, k][0] - 1)
                    countOdd, hitIDsOdd,templateArray,tracks,hitIDsEven=MonotonicityCheckOdd(countOdd, index,hitIDsOdd,hits,templateArray,tracks,hitIDsEven,t-1,k)

                elif (countOdd > 0 and templateArray[t, k][0] == 0):
                    if (countOdd < 3):
                        countOdd = 0
                        hitIDsOdd = []
                    else:
                        countOdd = 0
                        tracks.append(em.track(hitIDsOdd))
                        hitIDsOdd = []
                else:
                    hitIDsOdd = []

        if (countEven > 2):
            # newTrack=convertToTrack2(hitIDsEven,hits)
            tracks.append(em.track(hitIDsEven))
        if (countOdd > 2):
            # newTrack = convertToTrack2(hitIDsOdd,hits)
            tracks.append(em.track(hitIDsOdd))

    return tracks