# get the relevant dependencies 
import json
import os
import sys 
import numpy as np
import inspect, os.path
from math import pi, atan

# fileimport that should always work when executing this file
filename = inspect.getframeinfo(inspect.currentframe()).filename
file_path = os.path.dirname(os.path.abspath(filename))
project_root = os.path.dirname(file_path)

if project_root not in sys.path:
    sys.path.append(project_root)
from event_model import event_model as em


#### Interesting Implementations
# https://github.com/andreasfelix/hopfieldnetwork/tree/9e06c43c71c858afe4d8fc74f4c11c63ef83c22c

class smallHoppfiedNetwork():
    # takes 3 modules
    def __init__(self, m1, m2, m3) -> None:
        # m for module
        self.m1_hits = m1.hits()
        self.m2_hits = m2.hits()
        self.m3_hits = m3.hits()

        # c for count
        self.c1 = len(self.m1_hits)
        self.c2 = len(self.m2_hits)
        self.c3 = len(self.m3_hits)

        self.setup_neurons()
        self.init_weights()
        self.energy()

    def setup_neurons(self):
        # how are these orderd:
        # i propose in the order all hit_1_0 -> all hit_2_i 
        # TODO agreable?
        self.N1 = np.random.uniform(size = (self.c1 * self.c2)) # all connections between m1 & m2 -> dtype np.array size |m1| * |m2| 
        self.N2 = np.random.uniform(size = (self.c2 * self.c3)) # all connections between m2 & m3 -> dtype np.array size |m1| * |m2|
        
        # this is where we describe the angles in 3d space
        # I suppose lets go with the angles on xz and yz plane to describe the ordering of the direction in space (alt store a vector form one point to the other)
        # angles undirected so 0 - 180deg or 0 to pi 
        self.N1_info = np.zeros(shape = (self.c1 * self.c2, 2))
        self.N2_info = np.zeros(shape = (self.c2 * self.c3, 2))

        for i, hit1 in enumerate(self.m1_hits):
            for j, hit2 in enumerate(self.m2_hits):
                n_idx = i * self.c2 + j
                # this calculates angles of the line projections of the line segment defined by the hitpoints on zx, zy planes  
                angle_xz = atan((hit2.x - hit1.x)/(hit2.z - hit1.z))
                angle_yz = atan((hit2.y - hit1.y)/(hit2.z - hit1.z))

                self.N1_info[n_idx, 0] = angle_xz
                self.N1_info[n_idx, 1] = angle_yz
        
        for i, hit1 in enumerate(self.m2_hits):
            for j, hit2 in enumerate(self.m3_hits):
                n_idx = i * self.c3 + j
                # this calculates angles of the line projections of the line segment defined by the hitpoints on zx, zy planes  
                angle_xz = atan((hit2.x - hit1.x)/(hit2.z - hit1.z))
                angle_yz = atan((hit2.y - hit1.y)/(hit2.z - hit1.z))

                self.N2_info[n_idx, 0] = angle_xz
                self.N2_info[n_idx, 1] = angle_yz
                

    # fill the weight matrix based on the geometric properties of the hits / segments
    def init_weights(self):
        # its a bit like looking at small summatrices that need weights: |m1| segments of N1 far connected to |m2| segemtns in N2 
        self.weights = np.zeros(shape = (self.N1.shape[0], self.N2.shape[0])) 
        # i will init the weights as the angels between the lines -> just to see make the logic which ones to 
        # important here t handle the weights for segments that dont share a hit accordingly

        # con representing the idx of the hit in m2 that connects the segements in N1 and N2
        # intuitively i need to calc weights for all segments that connect via a point in m2
        for con in range(self.c2):
            # go over all segments (coming from hits in m1) that connect to the current con 
            for i in range(self.c1):
                n1_idx = i * self.c2 + con 
                # go over all segments coming from 'con' that connect to all hits in m3
                for j in range(self.c3): 
                    n2_idx = con * self.c3 + j
                    # this is the angle difference between the projected line segments in zx
                    self.weights[n1_idx, n2_idx] = abs(self.N1_info[n1_idx,0] - self.N1_info[n1_idx,1])
        pass
    
    def energy(self):
        # test calculation for the energy using the matrices 
        E = self.N1.T @ self.weights @ self.N2.T
        print(E)
        

    def update(self):
        pass

    def wibble_wabble(self):
        pass

    def tracks(self):
        pass


if __name__ == "__main__":
    # loading some event
    event_path = os.path.join(project_root, "events/bsphiphi/velo_event_12.json")
    f = open(event_path)
    json_data = json.loads(f.read())
    event = em.event(json_data)

    # taking a subsection of the event ( 3 modules to feed into the small hoppfield net)
    # lets only take the first 3 even modules
    m1 = event.modules[0]
    m2 = event.modules[2]
    m3 = event.modules[4]
    my_hoppfield = smallHoppfiedNetwork(m1=m1, m2=m2, m3=m3)

exit()

class HopfieldNetwork():
    def __init__(self, NEURONS, weights, ):
        self.NEURONS = NEURONS
        self.weights = weights

        self.n = NEURONS.shape[0] # Number of neurons
        self.steps = 0 # completed updated steps (one step is updating all neurons)

    def energy(self):
        return (-1) * self.goodness()

    def goodness(self):
        g = (self.weights.dot(self.NEURONS)).dot(self.NEURONS)
        return g/2
    

    # Stimpfl-Abele Paper - simple update rule (1) 
    # The update rule we want should be (16)
    def update_simple(self, steps=1):

        for i in range(steps):

            # For each step we go through all neurons (in random order) and update them
            for i in np.random.permutation(self.n): # is this really random? [TODO]
                
                # Local Update Rule 
                # Si = sign(sum_over_j(Tij * Sj))
                T_i = self.weights[i,:]

                # The Neurons Activation  - sum_over_j(Tij * Sj)
                neurons_activation = sum(T_i * self.NEURONS)

                # new updated value (either 1 or -1)
                S_i =  np.sign(neurons_activation)

                self.NEURONS[i] = S_i
                
            

            self.steps += 1

    def update(self, steps=1):

        for i in range(steps):
            pass



##############
# 1. Construct NEURONS
##############

# 1.1 Construct Track Segments
#       - "identified by a set of consecutive neurons connecting a trail of aligned hits" (Passaleva)
#       - "The angle among two neurons is defined by their scalar product;" (Passaleva)
#           + Can we even do this? With Segments close to the Vertex Collision Region?

# 1.2 Filter Track Segments

# 1.3 Initialize NEURONS 
#       - randomly?

neurons = np.array([1,1,1])


##############
# 2. Construct WEIGHTS
##############

weights = np.array([[0,3,2],
                    [3,0,-3],
                    [2,-3,0]])

##############
# find minumum
##############

hop_net = HopfieldNetwork(neurons, weights)

print('[Energy] ' + str(hop_net.energy()))


hop_net.update_simple(steps=3)
print('..updating')

print('[Energy] ' + str(hop_net.energy()))

print(hop_net.NEURONS)

