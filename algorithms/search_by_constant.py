import math
import numpy as np
from event_model import event_model as em


class SearchByConstant:

    def __init__(self, event):
        self.hits = sort_by_phi(event.hits) # sort hits by phi

    # method that checks previous tracks
    def solve(self):

        '***********Parameters ***************'
        NumberOfPreviousTracks = 1      # Parameter for number of previous tracks to check
        XminY = 0.01 # accepted deviation between the x and y ratio values of the track and a hit h
        YminZ = 0.01 # accepted deviation between the z and y ratio values of the track and a hit h
        XminZ = 0.01 # accepted deviation between the x and z ratio values of the track and a hit h
        angleDifference = 0.01 # accepted deviation between the polar angle in x,y values of the current track and a hit h
        trackGreaterThan = 1 # we only make the current_track a real track if its length is greater than this value
        beginOrEnd = 0 # set to 0 to compare hit polar angle to polar angle of the first hit of the track, or -1 to compare to last hit of the track
        '*************************************'
        
        tracks = [] # list of tracks found
        current_track = [] # list of hits representing the track currently under consideration
        current_track.append(self.hits[0]) # initialize by considering the first hit
        self.hits = self.hits[1:] # Take the first hit out of consideration for forming a track with itself
        skipped = 0 # variable to keep track of the number of modules skipped
        # KeepingTrackOfModules = np.zeros(NumberOfPreviousTracks,52)

     
        for h in self.hits: # loop over the hits in the event
            # if the hits polar angle is withing range of the track
            if abs(math.atan2(h.y, h.x) - math.atan2(current_track[beginOrEnd].y, current_track[beginOrEnd].x)) < angleDifference: 
                # We skipped 1 module in the track already, hit is in module after the first hit of the track
                if skipped == 1 and h.module_number == current_track[0].module_number + 2: 
                    current_track.insert(0, h) # add the hit after the first track hit
                # We skipped 1 module in the track already, hit is in module before the first hit of the track
                elif skipped == 1 and h.module_number == current_track[-1].module_number - 2:
                    current_track.append(h) # add the hit at the end of the track
                # We did not yet skip a module in the current track,  hit is in module before the first hit of the track
                elif h.module_number == current_track[0].module_number - 2:
                    current_track.insert(0, h) # add the hit after the first track hit
                # We did not yet skip a module in the current track,  hit is in module after the last hit of the track
                elif h.module_number == current_track[-1].module_number + 2:
                    current_track.append(h) # add the hit at the end of the track
                # We did not yet skip a module in the current track,  hit is 2 modules before the first hit of the track
                elif skipped <= 0 and h.module_number == current_track[0].module_number - 4:
                    current_track.insert(0, h) # add the hit 
                    skipped += 1 # increase skipped module counter
                # We did not yet skip a module in the current track,  hit is 2 modules after the last hit of the track
                elif skipped <= 0 and h.module_number == current_track[-1].module_number + 4:
                    current_track.append(h) # add the hit, skip next module counter
                    skipped += 1 # increase skipped module counter
                else:
                    hitHhasBeenAdded = False # Boolean to see whether the hit gets added to a track
                    for i in range(NumberOfPreviousTracks):
                        if(i < len(tracks)):
                            xDir,yDir,zDir = calculateDirectionVector(tracks[-1-i].hits[0],tracks[-1-i].hits[-1]) # direction values of the track

                            xRelation = h.x-tracks[-1-i].hits[0].x/xDir # difference in x values between track point and h,divided by x-direction of track
                            yRelation = h.y-tracks[-1-i].hits[0].y/yDir # difference in y values between track point and h,divided by y-direction of track
                            zRelation = h.z-tracks[-1-i].hits[0].z/zDir # difference in z values between track point and h,divided by z-direction of track

                            # for straight lines, xRelation, yRelation and zRelation should have a the same value
                            if(abs(xRelation-yRelation) <=0.01 and abs(yRelation-zRelation) <=0.01 and (xRelation-zRelation) <=0.01 ):
                                if(h.module_number > tracks[-1-i].hits[-1].module_number ): # hit h is after the current end of the track
                                    tracks[-1-i].hits.append(h) # add it at the end
                                    hitHhasBeenAdded = True
                                    break 
                                if(h.module_number < tracks[-1-i].hits[0].module_number ): # hit h is before the current start of the track
                                    tracks[-1-i].hits.insert(0,h) # add it at the beginning
                                    hitHhasBeenAdded = True
                                    break
                                if(tracks[-1-i].hits[0].module_number < h.module_number < tracks[-1-i].hits[-1].module_number ): # hit h is in between the current track ends
                                    tracks[-1-i].hits.insert(1,h) # add it somewhere in between
                                    hitHhasBeenAdded = True
                                    break

                    # When a hit has not been added to any track, use it as a starting point        
                    if(hitHhasBeenAdded == False):
                        skipped = 0 # reset skipped counter
                        if len(current_track) > trackGreaterThan: # If the current track has a minimal number of hits
                            tracks.append(em.track(current_track)) # add it to the list of tracks found
                        current_track = [] # Reset current track
                        current_track.append(h) # Add the hit as the initial hit of the track

            else: # if hit h is not sufficiently close to the current track
                hitHhasBeenAdded = False # Boolean to see whether the hit gets added to a track
                for i in range(NumberOfPreviousTracks):
                    if(i < len(tracks)):
                        xDir,yDir,zDir = calculateDirectionVector(tracks[-1-i].hits[0],tracks[-1-i].hits[-1]) # direction values of the track

                        xRelation = h.x-tracks[-1-i].hits[0].x/xDir # difference in x values between track point and h,divided by x-direction of track
                        yRelation = h.y-tracks[-1-i].hits[0].y/yDir # difference in y values between track point and h,divided by y-direction of track
                        zRelation = h.z-tracks[-1-i].hits[0].z/zDir # difference in z values between track point and h,divided by z-direction of track

                        # for straight lines, xRelation, yRelation and zRelation should have a the same value
                        if(abs(xRelation-yRelation) <= XminY and abs(yRelation-zRelation) <= YminZ and (xRelation-zRelation) <= XminZ ):
                            if(h.module_number > tracks[-1-i].hits[-1].module_number ): # hit h is after the current end of the track
                                tracks[-1-i].hits.append(h) # add it at the end
                                hitHhasBeenAdded = True
                                break 
                            if(h.module_number < tracks[-1-i].hits[0].module_number ): # hit h is before the current start of the track
                                tracks[-1-i].hits.insert(0,h) # add it at the beginning
                                hitHhasBeenAdded = True
                                break
                            if(tracks[-1-i].hits[0].module_number < h.module_number < tracks[-1-i].hits[-1].module_number ): # hit h is in between the current track ends
                                tracks[-1-i].hits.insert(1,h) # add it somewhere in between
                                hitHhasBeenAdded = True
                                break

                # When a hit has not been added to any track, use it as a starting point        
                if(hitHhasBeenAdded == False):
                    skipped = 0 # reset skipped counter
                    if len(current_track) > trackGreaterThan: # If the current track has a minimal number of hits
                        tracks.append(em.track(current_track)) # add it to the list of tracks found
                    current_track = [] # Reset current track
                    current_track.append(h) # Add the hit as the initial hit of the track

        if len(current_track) > trackGreaterThan:
            tracks.append(em.track(current_track)) # append final track if longer than 1

        return tracks

    def solve2(self):
        tracks = [] # list of tracks found
        current_track = [] # list of hits representing the track currently under consideration
        current_track.append(self.hits[0]) # initialize by considering the first hit
        self.hits = self.hits[1:] # Take the first hit out of consideration for forming a track with itself
        skipped = 0 # variable to keep track of the number of modules skipped

        for h in self.hits: # loop over the hits in the event
            # if the hits polar angle is withing range of the last hit in the track ( rather than first in first implementation)
            if abs(math.atan2(h.y, h.x) - math.atan2(current_track[-1].y, current_track[-1].x)) < 0.01: 
                # put the hits apart if track hits are too far apart
                if 0 <abs(h.module_number - current_track[0].module_number) < 13 : 
                    if(len(current_track) >= 2):
                        if( current_track[0].module_number < h.module_number < current_track[-1].module_number):
                            if(checkMonotonicity(current_track[0], h, current_track[-1]) and abs(calculateConstant(current_track[0] , h ) - calculateConstant( h,current_track[-1])) < 0.01):
                                current_track.insert(-2,h)
                            # else:
                            #     # check other track, add to new track
                        elif(current_track[0].module_number > h.module_number):
                            if(len(current_track) >= 2):
                                if(checkMonotonicity(h,current_track[0] , current_track[-1]) and abs(calculateConstant(current_track[0] , h ) - calculateConstant( h,current_track[-1])) < 0.01):
                                    current_track.insert(0,h)
                                # else:
                                #     # check other track, add to new track
                        elif(h.module_number > current_track[-1].module_number):
                            if(len(current_track) >= 2):
                                if(checkMonotonicity(current_track[0], current_track[-1],h ) and abs(calculateConstant(current_track[0] , h ) - calculateConstant( h,current_track[-1])) < 0.01):
                                   current_track.append(h)
                                # else:
                    elif(h.module_number == current_track[0].module_number + 2):
                        current_track.append(h)
                    elif(h.module_number == current_track[0].module_number - 2):
                        current_track.insert(0,h)
                    elif(h.module_number == current_track[0].module_number + 4):
                        current_track.append(h)
                    elif(h.module_number == current_track[0].module_number - 4):
                        current_track.insert(0,h)                    
                elif(0 == abs(h.module_number - current_track[0].module_number)):
                    # do not add hit to track
        # if len(current_track) > 1:
        #     tracks.append(em.track(current_track)) # append final track if longer than 1     
        return tracks

# sort the hits by polar angle
def sort_by_phi(hits):
    phis = []
    for h in hits:
        phis.append(math.atan2(h.y, h.x))
    sorted_index = np.argsort(phis)
    x = np.array(hits)[sorted_index]
    return x

# function that returns module number, needed for sort_by_moduleNumber
def moduleNumber(hit):
    return hit.module_number

# this sorts a set of hits by module number in ascending order
def sort_by_moduleNumber(hits):
    return hits.sort(key = moduleNumber)

# Calculate a proposed constant
def calculateConstant(hit1,hit2):
    return math.sqrt((hit1.x-hit2.x)^2+(hit1.y+hit2.y)^2)/math.sqrt((hit1.z-hit2.z)^2)

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
    return xDirection,yDirection,zDirection

# returns sum of absolute difference between two direction vectors between 2 pairs of hits 
def calculateDifferenceOfDirectionVector(hit1,hit2,hit3,hit4):
    dirX1, dirY1, dirZ1 = calculateDirectionVector(hit1, hit2)
    dirX2, dirY2, dirZ2 = calculateDirectionVector(hit3, hit4)

    difference = abs(dirX1-dirX2)+abs(dirY1-dirY2)+abs(dirZ1-dirZ2)

    return difference

# function to look ahead at multiple hits, perhaps something for later
def compareHitsByDirectionVector(hits) :
    temp = np.array(len(hits)-1)
    xDir12, yDir12,zDir12 = calculateDirectionVector(hits[0],hits[1])
    xDir13, yDir13,zDir13 = calculateDirectionVector(hits[1],hits[2])
    
    track1 = []
    track2 = []
    
    track1.append(hits[0])
    track1.append(hits[1])

    track2.append(hits[0])
    track2.append(hits[2])

    for i in range (2,len(hits)):
        xRelation12 = hits[i].x-hits[0].x/xDir12
        yRelation12 = hits[i].y-hits[0].y/yDir12
        zRelation12 = hits[i].z-hits[0].z/zDir12

        if(abs(xRelation12-yRelation12)<=0.01 and abs(yRelation12-zRelation12) <=0.01 and (xRelation12-zRelation12) <=0.01 ):
            track1.append(hits[i])

    for i in range (3,len(hits)):

        xRelation13 = hits[i].x-hits[0].x/xDir13
        yRelation13 = hits[i].y-hits[0].y/yDir13
        zRelation13 = hits[i].z-hits[0].z/zDir13  

        if(abs(xRelation13-yRelation13)<=0.01 and abs(yRelation13-zRelation13) <=0.01 and (xRelation13-zRelation13) <=0.01 ):
                track1.append(hits[i])   
    
    if(len(track1) < len(track2)):
        if(len(track2) >= 3):
            tracks.append(em.track(track2))
    else:
        if(len(track1 >= 3)):
            tracks.append(em.track(track1))
