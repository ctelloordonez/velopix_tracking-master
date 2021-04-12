import math


class template:
    def __init__(self, K):
        self.bins_count = K
        self.bins = [[] for i in range(self.bins_count)]

    def get_bins(self, event):
        for hit in event.hits:
            hit_phi = get_polar_angle(hit.x, hit.y)
            for b, bin in enumerate(self.bins):
                if b * (2 * math.pi / self.bins_count) < abs(hit_phi) <= (b+1) * (2 * math.pi / self.bins_count):
                    bin.append(hit)
        return self.bins


def get_polar_angle(x, y):
    phi = math.atan2(x, y)
    if (phi < 0):
        phi = math.pi - phi
    return phi

def get_polar_distance(x, y):
    r = math.sqrt(x ** 2 + y ** 2)
    return r

# Create function to sort each bin by z (or module)
def zSort():
    return

# Create template function
# Question: can we use module number of the hits?
def template():
    return

    
# Function to check whether hits have monotone distance as a function of z
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