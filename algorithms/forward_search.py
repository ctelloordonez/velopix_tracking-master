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
        # NumberOfPreviousTracks = 200     # Parameter for number of previous tracks to check
        XminY = 0.04 # accepted deviation between the x and y ratio values of the track and a hit h
        YminZ = 0.04 # accepted deviation between the z and y ratio values of the track and a hit h
        XminZ = 0.04 # accepted deviation between the x and z ratio values of the track and a hit h
        moduleDifferenceAllowed = 6
        '*************************************'
        usedHits = np.zeros(len(self.hits))
        tracks = [] # list of tracks found
        
        
        for i in range(len(self.hits)-1) : # loop over the hits in the event
            if(self.hits[i].module_number != self.hits[i+1].module_number and checkEven(self.hits[i]) == checkEven(self.hits[i+1]) ):
                if(abs(self.hits[i].module_number - self.hits[i+1].module_number) < moduleDifferenceAllowed):
                    if(usedHits[i] == 0 and usedHits[i+1] == 0):
                        xDirection,yDirection,zDirection = calculateDirectionVector(self.hits[i],self.hits[i+1])
                        checkForTrack2(i,self.hits,xDirection,yDirection,zDirection,tracks,usedHits,XminY,YminZ,XminZ)
                        # checkForTrack2(i,self.hits,xDirection,yDirection,zDirection,tracks,usedHits,XminY,YminZ,XminZ)
        
        tracks = combineTracks(tracks,5)
        tracks = combineTracks(tracks,5)
        tracks = combineTracks(tracks,3)
        # tracks = combineTracks(tracks,2)

        # tracks = combineTracks2(tracks,3, 0.0011)
        # tracks = combineTracks2(tracks,2,0.0009)
        # tracks = combineTracks2(tracks,1,0.0007)
        # tracks = combineTracks2(tracks,3,0.0011)
        # tracks = combineTracks2(tracks,2,0.0009)
        # tracks = combineTracks2(tracks,1,0.0007)
        
        tracks = removeTracks(tracks,3)
        # tracks = takeOutNonConsecutiveTracks(tracks)
        # tracks = takeOutDuplicateHits(tracks)
        # for i in tracks[1].hits:
        #     print(i.id)
        
        return tracks

    def solve2(self):

        '***********Parameters ***************'
        XminY = 0.05 # accepted deviation between the x and y ratio values of the track and a hit h #0.05 works well
        YminZ = 0.05 # accepted deviation between the z and y ratio values of the track and a hit h
        XminZ = 0.05 # accepted deviation between the x and z ratio values of the track and a hit h
        moduleDifferenceAllowed = 6
        '*************************************'
        usedHits = np.zeros(len(self.hits))
        tracks = [] # list of tracks found
        # current_track = [] # list of hits representing the track currently under consideration
        # current_track.append(self.hits[0]) # initialize by considering the first hit
        # self.hits = self.hits[1:] # Take the first hit out of consideration for forming a track with itself
        skipped = 0 # variable to keep track of the number of modules skipped
        
        for i in range(len(self.hits)-1) : # loop over the hits in the event
            if(self.hits[i].module_number != self.hits[i+1].module_number and checkEven(self.hits[i]) == checkEven(self.hits[i+1]) ):
                if(abs(self.hits[i].module_number - self.hits[i+1].module_number) <moduleDifferenceAllowed):
                    if(usedHits[i] == 0 and usedHits[i+1] == 0):
                        xDirection,yDirection,zDirection = calculateDirectionVector(self.hits[i],self.hits[i+1])
                        # checkForTrack(i,self.hits,xDirection,yDirection,zDirection,tracks,usedHits,XminY,YminZ,XminZ)
                        tracks,usedHits = checkForTrack2(i,self.hits,xDirection,yDirection,zDirection,tracks,usedHits,XminY,YminZ,XminZ)
                        
                    elif(usedHits[i]==0):
                        temp = i+1
                        while(usedHits[temp]==1 and temp < len(self.hits)-1 ):
                            temp = temp+1
                        if(usedHits[temp] == 0):
                            Direction,yDirection,zDirection = calculateDirectionVector(self.hits[i],self.hits[temp])
                            # checkForTrack(i,self.hits,xDirection,yDirection,zDirection,tracks,usedHits,XminY,YminZ,XminZ)
                            tracks,Usedhits = checkForTrack2(i,self.hits,xDirection,yDirection,zDirection,tracks,usedHits,XminY,YminZ,XminZ)
                            
                    elif(usedHits[i+1]==0):
                        temp = i-1
                        while(usedHits[temp]==1 and temp >0 ):
                            temp = temp-1
                        if(usedHits[temp] == 0):
                            Direction,yDirection,zDirection = calculateDirectionVector(self.hits[temp],self.hits[i+1])
                            # checkForTrack(i,self.hits,xDirection,yDirection,zDirection,tracks,usedHits,XminY,YminZ,XminZ)
                            tracks,Usedhits = checkForTrack2(i,self.hits,xDirection,yDirection,zDirection,tracks,usedHits,XminY,YminZ,XminZ)
            # if(i==len(self.hits)-2):
            #     print(usedHits)

        # tracks = takeOutNonConsecutiveTracks(tracks)

        # tracks = combineTracks(tracks,4)
        # tracks = combineTracks(tracks,3)
        # tracks = combineTracks(tracks,3)
        
        tracks = combineTracks(tracks,3)
        tracks = combineTracks(tracks,4)


        # tracks = combineTracks2(tracks,3)
        # tracks = combineTracks2(tracks,2)
        # tracks = combineTracks2(tracks,1)
        # tracks = combineTracks2(tracks,3)
        # tracks = combineTracks2(tracks,2)
        # tracks = combineTracks2(tracks,1)
        
        # tracks = takeOutNonConsecutiveTracks(tracks)
        # tracks= takeOutDuplicateHits(tracks)
        tracks = removeTracks(tracks,3)
        # print(usedHits)
        return tracks

    def solve3(self):
        tracks = []
        currentTrack = []
        for i in range(len(self.hits)-1):
            contributed = 0
            for j in range(1,2): # 4 does well
                if(i <len(self.hits)-j and contributed == 0):
                    if(0<abs(self.hits[i].module_number -self.hits[i+j].module_number) <4 ):
                        currentTrack.append(self.hits[i])
                        currentTrack.append(self.hits[i+j])
                        tracks.append(em.track(currentTrack))
                        currentTrack = []
                        contributed = 1

            
           
        return tracks

def checkForTrack(i,hits,xDir,yDir,zDir,tracks,usedHits,XminY,YminZ,XminZ):

    indexList = []
    indexList.append(i)
    indexList.append(i+1)

    for j in range(1,6):  # 30 good # 8 gave 72% #6 provided 70 on large data set.
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

        # skipped = False
        # countEven = 0
        # countOdd = 0
   
        # for k in range(len(followingModuleNumbers)):
        #     if(followingModuleNumbers[k] == 1):
        #         temp = k+2
        #         while(followingModuleNumbers[temp] == 1 and count <3 and temp <= len(followingModuleNumbers)-3):
        #             count = count + 1
        #             temp= temp+2
        #             if(skipped == False and followingModuleNumbers[temp] == 0 ):
        #                 skipped = True
        #                 count = count + 1
        #             elif(count >= 3):
        #                 tracks.append(em.track(tempTrack))
                
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
                if(abs(MonotoneTrackT-tempMonotone) < 0.0011): #0.0011 and 0.0015 wok well. Increase value, more ghost, fewer clones
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

            firstHitSame = (tracks[t].hits[0] != tracks[t+temp].hits[0])
            lastHitSame = (tracks[t].hits[-1] != tracks[t+temp].hits[-1])

            if(abs(MonotoneTrackT-tempMonotone) < tolerance and firstHitSame and lastHitSame):
                for z in range(len(tracks[t+temp].hits)):
                    tracks[t].hits.append(tracks[t+temp].hits[z]) # add it at the end
                del tracks[t+temp]
                
                length = length - 1
    return tracks

def combineTracks3(tracks,lookAhead):
    length = len(tracks)
    for t in range(len(tracks)-1):
        temp = 1
        if(t+temp< length):
            MonotoneTrackT = calculateMonotoneApproximation(tracks[t].hits[0],tracks[t].hits[-1])
            tempMonotone = calculateMonotoneApproximation(tracks[t+temp].hits[0],tracks[t+temp].hits[-1])

            firstHitSame = (tracks[t].hits[-1] != tracks[t+temp].hits[-1])
            lastHitSame = (tracks[t].hits[-1] != tracks[t+temp].hits[-1])
            if(abs(MonotoneTrackT-tempMonotone) < 0.0013 and firstHitSame and lastHitSame):
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

# def sort_by_direction(tracks):
#     phis = []

#     for t in tracks:
#         if(tracks[t].hits[0].module_number > tracks[t].hits[-1].module_number):
#             phis.append((tracks[t].hits[0].x-tracks[t].hits[-1].x)
#         else:
#             phis.append((tracks[t].hits[-1].x-tracks[t].hits[0].x)

#     sorted_index = np.argsort(phis)
#     x = np.array(tracks)[sorted_index]
#     return x

def checkForTrack2(i,hits,xDir,yDir,zDir,tracks,usedHits,XminY,YminZ,XminZ):

    indexList = []
    indexList.append(i)
    indexList.append(i+1)

    for j in range(1,8):  # 30 good # 8 gave 72s
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
    
def matchTracks(tracks,numberOfTracks,MonotoneDifference):
    for i in range(len(tracks)):
        for j in range(numberOfTracks):
             if(abs(calculateMonotoneApproximation(tracks[i].hits[0],tracks[i+1+j].hits[0]) - calculateMonotoneApproximation(tracks[i].hits[-1],tracks[i+1+j].hits[-1])) < MonotoneDifference ):
                for z in range(len(tracks[i+1+j].hits)):
                    tracks[i].hits.append(tracks[i+1+j].hits[0]) # add it at the end