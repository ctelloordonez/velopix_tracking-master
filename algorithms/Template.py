import math
import copy

from event_model import *
from event_model.event_model import track

# Idea to be created: if there exist three or more hits at a line, the template should true, there is a (potential track here)
def template(panels,slopeX,slopeY,startingPanel,trackLength,error,distancePermitted):
    hitSum = 0 #count the number of panels in which a valid hit occurs
    numberOfConsecutivePanelsWithoutAHit = 0
    allowedConsecutiveNonHits = 3
    # Checking consecutive hits and number of hits
    for i in range(trackLength):
        if(numberOfConsecutivePanelsWithoutAHit < allowedConsecutiveNonHits):
            xValue = 0.0
            yValue = 0.0
            if(checkPanelRegion(panels[startingPanel+i],xValue,yValue,distancePermitted
            ) == True):
                hitSum = hitSum + 1
            else:
                numberOfConsecutivePanelsWithoutAHit += 1
        else:
            return False
    if(numberOfConsecutivePanelsWithoutAHit < allowedConsecutiveNonHits):
        if(hitSum >= trackLength - error and hitsum >= 3):
            return True
        else:
            return False
    else:
        return False

def calculateXvalue():
    return
def calculateYvalue():
    return


 Function to check whether hits have monotone distance as a function of z
# Let hits be an array of hits that have been identified as potentially forming a track
def checkMonotoneDistance(hits):
    i=0
    j=1

    distanceI = get_polar_distance(hits[i].x,hits[i].y)
    distanceJ = get_polar_distance(hits[j].x,hits[j].y)

    while(distanceI  == distanceJ and j != length(hits)-1):
        i = i+1
        j = j+1

        distanceI = distanceJ
        distanceJ = get_polar_distance(hits[j].x,hits[j].y)

    if(j == length(hits)-1 and distanceI  ==distanceJ):
        return True

    elif(distanceI y > distanceJ):
        for k in range(i,length(hits)-1):
            if(get_polar_distance(hits[k].x,hits[k].y) < get_polar_distance(hits[k+1].x,hits[k+1].y)):
                return False

    elif(distanceI  < distanceJ):
        for k in range(i,length(hits)-1):
            if(get_polar_distance(hits[k].x,hits[k].y) > get_polar_distance(hits[k+1].x,hits[k+1].y)):
                return False
    else:
        return True

# Function to check whether hits have monotone distance as a function of z
# Let hits be an array of hits that have been identified as potentially forming a track
def checkMonotoneTime(hits):
    i=0
    j=1

    while(hits[i].t == hits[j].t and j != length(hits)-1):
        i = i+1
        j = j+1

    if(j == length(hits)-1 and hits[i].t== hits[j].t):
        return True

    elif(hits[i].t > hits[j].t):
        for k in range(i,length(hits)-1):
            if(hits[k].te < hits[k+1].t):
                return False

    elif(hits[i].t < hits[j].t):
        for k in range(i,length(hits)-1):
            if(hits[k].t > hits[k+1].t):
                return False
    else:
        return True

# Would be better if we use binary search
def checkPanelRegion(panel,pointX,pointY, distance):
    numberOfHits = length(panel) # number of hits in the panel
    for j in (numberOfHits):
        difference = math.sqrt((panel[j].x-pointX)^2+(panel[j].y-pointY)^2) # calculate how far the pixel hit is from the predicited point
        if(difference < distance):
            return True
        elif(j = numberOfHits-1):
            return False

# Function to determine the amount of feasible hints in the x,y region of the panel.
# if this is more than 2, presumably the template is not specific enough yet
def checkNumberOfFeasibleHitsInPanel(panel,xUpper,xLower,yUpper,yLower):
    feasibleHits = 0
    # Binary search the panel for feasible points