import math
import numpy as np
from event_model import event_model as em


class ForwardSearch:

    def __init__(self, event):
        self.hits = sort_by_phi(event.hits) # sort hits by phi
        # self.hits = sort_by_phi_projected(event.hits) # sort by projected phi into the plane
    
    # Main method, to go over a list of hits sorted by phi.
    # It then uses the direction vector of two consequtive hits to see if there are hits around these two hits that fall on this vector.
    # If so, these are added to the a track together with the 2 initial hits, and then these hits are flagged as being used.
    def solve(self):

        '*********** Parameters **************'
        # NumberOfPreviousTracks = 200     # Parameter for number of previous tracks to check
        XminY = 0.13 # accepted deviation between the x and y ratio values of the track and a hit h
        YminZ = 0.13# accepted deviation between the z and y ratio values of the track and a hit h
        XminZ = 0.13 # accepted deviation between the x and z ratio values of the track and a hit h
        moduleDifferenceAllowed = 6 # How far apart can hits be to allow them to form a initial direction vector
        tolerance = 0.0013 # monotone deviation allowed when combining tracks to reduce ghosts
        lookAround = 6 # number of hits to look at around the direction vector hits # 10 standard
        '*************************************'
        usedHits = np.zeros(len(self.hits)) # array to track used hits, if hit used in track set index of hit to 1.
        tracks = [] # list of tracks found
        
        
        for i in range(len(self.hits)-1) : # loop over the hits in the event
            if(self.hits[i].module_number != self.hits[i+1].module_number and checkEven(self.hits[i]) == checkEven(self.hits[i+1]) ):
                if(abs(self.hits[i].module_number - self.hits[i+1].module_number) < moduleDifferenceAllowed):
                    if(usedHits[i] == 0 and usedHits[i+1] == 0):
                        xDirection,yDirection,zDirection = calculateDirectionVector(self.hits[i],self.hits[i+1])
                        checkForTrack3(i,self.hits,xDirection,yDirection,zDirection,tracks,usedHits,XminY,YminZ,XminZ,lookAround)
        
        # combineTracks reduces clone tracks
        # tracks = combineTracks2(tracks,5,tolerance)
        # tracks = combineTracks2(tracks,5,tolerance)
        # tracks = combineTracks2(tracks,3,tolerance)  


        # tracks = combineTracks3(tracks,5,tolerance)
        # tracks = combineTracks3(tracks,5,tolerance)
        # tracks = combineTracks3(tracks,3,tolerance) 

        # #check for gaps reduces ghost tracks
        tracks = checkForGaps(tracks)
        return tracks
    
    # Similar to the first solve, now with the ability to use a hit if its neighbour has already been used.
    def solve2(self):

        '***********Parameters ***************'
        XminY = 0.09 # accepted deviation between the x and y ratio values of the track and a hit h #0.05 works well
        YminZ = 0.09 # accepted deviation between the z and y ratio values of the track and a hit h
        XminZ = 0.09# accepted deviation between the x and z ratio values of the track and a hit h
        moduleDifferenceAllowed = 6
        tolerance = 0.004
        lookAround = 10
        '*************************************'
        usedHits = np.zeros(len(self.hits)) # array to track used hits, if hit used in track set index of hit to 1.
        tracks = [] # list of tracks found
        skipped = 0 # variable to keep track of the number of modules skipped
        
        for i in range(len(self.hits)-1) : # loop over the hits in the event
            if(self.hits[i].module_number != self.hits[i+1].module_number and checkEven(self.hits[i]) == checkEven(self.hits[i+1]) ):
                if(abs(self.hits[i].module_number - self.hits[i+1].module_number) <moduleDifferenceAllowed):
                    if(usedHits[i] == 0 and usedHits[i+1] == 0):
                        xDirection,yDirection,zDirection = calculateDirectionVector(self.hits[i],self.hits[i+1])
                        tracks,usedHits = checkForTrack2(i,self.hits,xDirection,yDirection,zDirection,tracks,usedHits,XminY,YminZ,XminZ,lookAround)
                        
                    elif(usedHits[i]==0):
                        temp = i+1
                        while(usedHits[temp]==1 and temp < len(self.hits)-1 ):
                            temp = temp+1
                        if(usedHits[temp] == 0):
                            Direction,yDirection,zDirection = calculateDirectionVector(self.hits[i],self.hits[temp])
                            tracks,Usedhits = checkForTrack2(i,self.hits,xDirection,yDirection,zDirection,tracks,usedHits,XminY,YminZ,XminZ,lookAround)
                            
                    elif(usedHits[i+1]==0):
                        temp = i-1
                        while(usedHits[temp]==1 and temp >0 ):
                            temp = temp-1
                        if(usedHits[temp] == 0):
                            Direction,yDirection,zDirection = calculateDirectionVector(self.hits[temp],self.hits[i+1])
                            tracks,Usedhits = checkForTrack2(i,self.hits,xDirection,yDirection,zDirection,tracks,usedHits,XminY,YminZ,XminZ,lookAround)

        # tracks = combineTracks2(tracks,4,tolerance)
        # tracks = combineTracks2(tracks,3,tolerance)
        # tracks = combineTracks2(tracks,3,tolerance)
        

        # track = checkForGaps(tracks)

        # tracks = combineTracks2(tracks,3,tolerance)
        # tracks = combineTracks2(tracks,4,tolerance)

        tracks = combineTracks3(tracks,5,tolerance)
        tracks = combineTracks3(tracks,5,tolerance)
        tracks = combineTracks3(tracks,3,tolerance) 

        tracks = checkForGaps(tracks)
        return tracks

    # This is just a method that blindly adds tracks of 2 that have a similar angle. Just to see how many we can find.
    def solve3(self):

        tracks = []
        currentTrack = []
        for i in range(len(self.hits)-1):
            # contributed = 0
            for j in range(1,150): # 4 does well
                if(i <len(self.hits)-j):
                    if(0<abs(self.hits[i].module_number -self.hits[i+j].module_number) <6 ):
                        currentTrack.append(self.hits[i])
                        currentTrack.append(self.hits[i+j])
                        tracks.append(em.track(currentTrack))
                        currentTrack = []
                        # contributed = 1


            
           
        return tracks

def checkForTrack(i,hits,xDir,yDir,zDir,tracks,usedHits,XminY,YminZ,XminZ):

    indexList = []
    indexList.append(i)
    indexList.append(i+1)

    for j in range(1,8):  # 30 good # 8 gave 72% #6 provided 70 on large data set.
        if( i-j >= 0):
            if(xDir != 0 and yDir != 0 and zDir != 0 and usedHits[i-j]==0):
                xRelation = (hits[i-j].x - hits[i].x)/xDir # difference in x values between track point and h,divided by x-direction of track
                yRelation = (hits[i-j].y - hits[i].y)/yDir # difference in y values between track point and h,divided by y-direction of track
                zRelation = (hits[i-j].z - hits[i].z)/zDir # difference in z values between track point and h,divided by z-direction of track

            # for straight lines, xRelation, yRelation and zRelation should have a the same value
                if(abs(xRelation-yRelation) <=XminY and abs(yRelation-zRelation) <=YminZ and (xRelation-zRelation) <=XminZ ):
                    indexList.insert(0,i-j)
                    usedHits[i-j]=1

        if( i+1+j <= len(hits)-1):
            if(xDir != 0 and yDir != 0 and zDir != 0 and usedHits[i+1+j]==0 ):
                xRelation = (hits[i+1+j].x - hits[i+1].x)/xDir # difference in x values between track point and h,divided by x-direction of track
                yRelation = (hits[i+1+j].y - hits[i+1].y)/yDir # difference in y values between track point and h,divided by y-direction of track
                zRelation = (hits[i+1+j].z - hits[i+1].z)/zDir # difference in z values between track point and h,divided by z-direction of track

                # for straight lines, xRelation, yRelation and zRelation should have a the same value
                if(abs(xRelation-yRelation) <=XminY and abs(yRelation-zRelation) <=YminZ and (xRelation-zRelation) <=XminZ ):
                    indexList.append(i+1+j)
                    usedHits[i-j]=1
    
    if(len(indexList)> 2):
        usedHits[i]=1
        usedHits[i+1]=1
        tempTrack = []
        followingModuleNumbers = np.zeros(52)
        for t in range(len(indexList)):
            tempTrack.append(hits[indexList[t]])
            # usedHits[t] = 1
        tracks.append(em.track(tempTrack))
                
def combineTracks(tracks,lookAhead):
    length = len(tracks)
    for t in range(len(tracks)-1):
        for k in range(1,lookAhead): # 5 works good
            if(t+k < length):
                MonotoneTrackT = calculateMonotoneApproximation(tracks[t].hits[0],tracks[t].hits[-1])
                tempMonotone = calculateMonotoneApproximation(tracks[t+k].hits[0],tracks[t+k].hits[-1])

            
                # firstHitSame = (tracks[t].hits[0] != tracks[t+k].hits[0])
                # lastHitSame = (tracks[t].hits[-1] != tracks[t+k].hits[-1])
                # if(abs(MonotoneTrackT-tempMonotone) < 0.0015 and firstHitSame and lastHitSame):
                # Increase value: more ghost, fewer clones , decrease value: fewer ghost, more clones
                if(abs(MonotoneTrackT-tempMonotone) < 0.0007): #0.0011 and 0.0015 wok well. Increase value, more ghost, fewer clones
                    for z in range(len(tracks[t+k].hits)):
                        tracks[t].hits.append(tracks[t+k].hits[z]) # add it at the end
                    del tracks[t+k]
                    length = length - 1
    return tracks

def combineTracks2(tracks,lookAhead,tolerance):
    length = len(tracks)
    for t in range(len(tracks)-1):
        temp = 1
        if(t+temp< length):
            MonotoneTrackT = calculateMonotoneApproximation(tracks[t].hits[0],tracks[t].hits[-1])
            tempMonotone = calculateMonotoneApproximation(tracks[t+temp].hits[0],tracks[t+temp].hits[-1])

            if(abs(MonotoneTrackT-tempMonotone) < tolerance):
                for z in range(len(tracks[t+temp].hits)):
                    tracks[t].hits.append(tracks[t+temp].hits[z]) # add it at the end
                del tracks[t+temp]
                
                length = length - 1
    return tracks

def combineTracks3(tracks,lookAhead,tolerance):
    length = len(tracks)

    for t in range(len(tracks)-1):
        temp = 1
        if(t+temp< length):
            MonotoneTrackT = calculateMonotoneApproximation(tracks[t].hits[0],tracks[t].hits[-1])
            tempMonotone = calculateMonotoneApproximation(tracks[t+temp].hits[0],tracks[t+temp].hits[-1])

            firstHitSame = (tracks[t].hits[-1].module_number != tracks[t+temp].hits[-1].module_number)
            lastHitSame = (tracks[t].hits[-1].module_number!= tracks[t+temp].hits[-1].module_number)
            if(abs(MonotoneTrackT-tempMonotone) < tolerance and firstHitSame and lastHitSame):
                for z in range(len(tracks[t+temp].hits)):
                    tracks[t].hits.append(tracks[t+temp].hits[z]) # add it at the end
                del tracks[t+temp]
                
                length = length - 1
    return tracks

def combineTracks4(tracks,lookAhead):
    length = len(tracks)
  
    for t in range(len(tracks)-1):
        toBeRemoved = []
        for k in range(1,lookAhead): # 5 works good
            if(t+k < length):
                MonotoneTrackT = calculateMonotoneApproximation(tracks[t].hits[0],tracks[t].hits[-1])
                tempMonotone = calculateMonotoneApproximation(tracks[t+k].hits[0],tracks[t+k].hits[-1])

            
                # firstHitSame = (tracks[t].hits[0] != tracks[t+k].hits[0])
                # lastHitSame = (tracks[t].hits[-1] != tracks[t+k].hits[-1])
                # if(abs(MonotoneTrackT-tempMonotone) < 0.0015 and firstHitSame and lastHitSame):
                if(abs(MonotoneTrackT-tempMonotone) < 0.0015): #0.0011 
                    for z in range(len(tracks[t+k].hits)):
                        tracks[t].hits.append(tracks[t+k].hits[z]) # add it at the end
                        toBeRemoved.append(t+k)
        for i in toBeRemoved:
            del i
            # print(i)
            length = length -1
    return tracks

# def combineTracksViaDVector(tracks,NumberOfTracksToCheck,tolerance):
#     length = len(tracks)
#     for t in range(len(tracks)):
#         if(t < length):
#             xChange,yChange,zChange = calculateDirectionVector(tracks[t].hits[0],tracks[t].hits[1])
#             for k in range(1,NumberOfTracksToCheck):
#                 tempX,tempY,TempZ = calculateDirectionVector(tracks[t+k].hits[0],tracks[t+k].hits[0])
#                 xSquaredDif = (xChange-tempX)**2
#                 ySquaredDif = (yChange-tempY)**2
#                 zSquaredDif = (zChange-tempZ)**2
#                 if(math.sqrt(xSquaredDif+ySquaredDif+zSquaredDif) < tolerance):

#     return tracks


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


def checkForTrack2(i,hits,xDir,yDir,zDir,tracks,usedHits,XminY,YminZ,XminZ,lookAround):

    indexList = []
    indexList.append(i)
    indexList.append(i+1)

    for j in range(1,lookAround):  # 30 good # 8 gave 72s
        if( i-j >= 0):
            if(xDir != 0 and yDir != 0 and zDir != 0 and usedHits[i-j]==0):
                xRelation = (hits[i-j].x - hits[i].x)/xDir # difference in x values between track point and h,divided by x-direction of track
                yRelation = (hits[i-j].y - hits[i].y)/yDir # difference in y values between track point and h,divided by y-direction of track
                zRelation = (hits[i-j].z - hits[i].z)/zDir # difference in z values between track point and h,divided by z-direction of track

            # for straight lines, xRelation, yRelation and zRelation should have a the same value
                if(abs(xRelation-yRelation) <=XminY and abs(yRelation-zRelation) <=YminZ and (xRelation-zRelation) <=XminZ ):
                    indexList.insert(0,i-j)
                    # usedHits[i-j]=1

        if( i+1+j <= len(hits)-1):
            if(xDir != 0 and yDir != 0 and zDir != 0 and usedHits[i+1+j]==0 ):
                xRelation = (hits[i+1+j].x - hits[i+1].x)/xDir # difference in x values between track point and h,divided by x-direction of track
                yRelation = (hits[i+1+j].y - hits[i+1].y)/yDir # difference in y values between track point and h,divided by y-direction of track
                zRelation = (hits[i+1+j].z - hits[i+1].z)/zDir # difference in z values between track point and h,divided by z-direction of track

                # for straight lines, xRelation, yRelation and zRelation should have a the same value
                if(abs(xRelation-yRelation) <=XminY and abs(yRelation-zRelation) <=YminZ and (xRelation-zRelation) <=XminZ ):
                    indexList.append(i+1+j)
                    # usedHits[i-j]=1
    
    if(len(indexList)> 2):
        # usedHits[i]=1
        # usedHits[i+1]=1
        tempTrack = []
        for t in range(len(indexList)):
            tempTrack.append(hits[indexList[t]])
            usedHits[indexList[t]] = 1
        tracks.append(em.track(tempTrack))
    return tracks, usedHits

def checkForTrack3(i,hits,xDir,yDir,zDir,tracks,usedHits,XminY,YminZ,XminZ,lookAround):

    indexList = []
    indexList.append(i)
    indexList.append(i+1)

    for j in range(1,lookAround):  # 30 good # 8 gave 72s
        if( i-j >= 0):
            if(xDir != 0 and yDir != 0 and zDir != 0 and usedHits[i-j]==0):
                xRelation = (hits[i-j].x - hits[i].x)/xDir # difference in x values between track point and h,divided by x-direction of track
                yRelation = (hits[i-j].y - hits[i].y)/yDir # difference in y values between track point and h,divided by y-direction of track
                zRelation = (hits[i-j].z - hits[i].z)/zDir # difference in z values between track point and h,divided by z-direction of track

            # for straight lines, xRelation, yRelation and zRelation should have a the same value
                if(abs(xRelation-yRelation) <=XminY and abs(yRelation-zRelation) <=YminZ and (xRelation-zRelation) <=XminZ ):
                    indexList.insert(0,i-j)
                    # usedHits[i-j]=1

        if( i+1+j <= len(hits)-1):
            if(xDir != 0 and yDir != 0 and zDir != 0 and usedHits[i+1+j]==0 ):
                xRelation = (hits[i+1+j].x - hits[i].x)/xDir # difference in x values between track point and h,divided by x-direction of track
                yRelation = (hits[i+1+j].y - hits[i].y)/yDir # difference in y values between track point and h,divided by y-direction of track
                zRelation = (hits[i+1+j].z - hits[i].z)/zDir # difference in z values between track point and h,divided by z-direction of track

                # for straight lines, xRelation, yRelation and zRelation should have a the same value
                if(abs(xRelation-yRelation) <=XminY and abs(yRelation-zRelation) <=YminZ and (xRelation-zRelation) <=XminZ ):
                    indexList.append(i+1+j)
                    # usedHits[i-j]=1

    if(len(indexList)> 2):
        # usedHits[i]=1
        # usedHits[i+1]=1
        tempTrack = []
        for t in range(len(indexList)):
            tempTrack.append(hits[indexList[t]])
            usedHits[indexList[t]] = 1
        tracks.append(em.track(tempTrack))
    return tracks, usedHits

# function that returns module number, needed for sort_by_moduleNumber
def moduleNumber(hit):
    return hit.module_number

# this sorts a set of hits by module number in ascending order
def sort_by_moduleNumber(hits):
    return hits.sort(key = moduleNumber)

def calculateMonotoneApproximation(hit1,hit2):
    if(hit1.module_number > hit2.module_number):
        x = (polarDistance(hit1) - polarDistance(hit2))/abs(hit1.z-hit2.z)
        return x
    if(hit1.module_number < hit2.module_number):
        x = (polarDistance(hit2) - polarDistance(hit1))/abs(hit2.z-hit1.z)
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
        xDirection = hit1.x
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
    
def takeOutNonConsecutiveTracks(tracks):
    
    for t in tracks:
        modulesHit = np.zeros(52)
        for i in range(len(t.hits)):
            modulesHit[t.hits[i].module_number-1] = 1
        skipped = 0
        countEven = 0
        countOdd = 0
        for k in range(len(modulesHit)):
            if(countEven <  3 and countOdd < 3):
                if(k %2 ==0):
                    if(modulesHit[k]==1 and countEven < 3):
                        countEven += 1
                    elif(countEven == 3):
                        break
                    elif(countEven > 0 and modulesHit[k]==0 and skipped == 0):
                        skipped +=1
                    elif(countEven> 0 and modulesHit[k]==0 and skipped > 0):
                        countEven = 0
                        skipped = 0
                else:
                    if(modulesHit[k]==1 and countOdd < 3):
                        countOdd += 1
                    elif(countOdd == 3):
                        break
                    elif(countOdd > 0 and modulesHit[k]==0 and skipped == 0):
                        skipped +=1
                    elif(countOdd> 0 and modulesHit[k]==0 and skipped > 0):
                        countOdd = 0
                        skipped = 0
        
        if(countEven <  3 and countOdd < 3):
            del t
    return tracks

def takeOutDuplicateHits(tracks):
    for t in tracks:
        for i in range(len(t.hits)-1):
            length = len(t.hits)
            for j in range(i,len(t.hits)):
                if(j < length and t.hits[i].id == t.hits[j].id):
                    del t.hits[j]
                    length = length -1
            
    return tracks
def matchTracks(tracks,numberOfTracks,MonotoneDifference):
    for i in range(len(tracks)):
        for j in range(numberOfTracks):
             if(abs(calculateMonotoneApproximation(tracks[i].hits[0],tracks[i+1+j].hits[0]) - calculateMonotoneApproximation(tracks[i].hits[-1],tracks[i+1+j].hits[-1])) < MonotoneDifference ):
                for z in range(len(tracks[i+1+j].hits)):
                    tracks[i].hits.append(tracks[i+1+j].hits[0]) # add it at the end

def checkForGaps(tracks):
    
    toBeRemoved = []
    for i in range(len(tracks)):
        temporary = np.zeros(52)
        for h in tracks[i].hits:
            if(temporary[h.module_number] != 0):
                toBeRemoved.append(i)
            else:
                temporary[h.module_number] = 1
        
        validTrack = False
        for i in range(52):
            if(i < 49):
                sum = temporary[i]+temporary[i+1]+temporary[i+2] +temporary[i+3]
                if(sum>=3):
                    validTrack = True
        if (validTrack == False):
            toBeRemoved.append(i)

    copy = []
    
    for x in range(len(tracks)):
        notToBeAdded = False
        for j in range(len(toBeRemoved)):
            if(x==toBeRemoved[j]):
                notToBeAdded =True
        if(notToBeAdded == False):
            copy.append(tracks[x])
    return copy

