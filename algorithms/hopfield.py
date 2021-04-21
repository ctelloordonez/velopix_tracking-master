import numpy as np


import json
import os
import sys 


module_path = os.path.abspath(os.path.join('..'))
if module_path not in sys.path:
    sys.path.append(module_path)

import event_model.event_model as em

#### Interesting Implementations
# https://github.com/andreasfelix/hopfieldnetwork/tree/9e06c43c71c858afe4d8fc74f4c11c63ef83c22c


##############
# 0. UTILS
##############

# Functions to calculate Polar angles
def cart2pol(x, y):
    rho = np.sqrt(x**2 + y**2)
    phi = np.arctan2(y, x)
    return(rho, phi)

def cart2polAngle(x, y):
    phi = np.arctan2(y, x)
    return phi

def pol2cart(rho, phi):
    x = rho * np.cos(phi)
    y = rho * np.sin(phi)
    return(x, y)




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

# -------------------------------------------------------------------------------------
# NOTES:
#
# 1.1 Construct Track Segments
#       - "identified by a set of consecutive neurons connecting a trail of aligned hits" (Passaleva)
#       - "The angle among two neurons is defined by their scalar product;" (Passaleva)
#           + Can we even do this? With Segments close to the Vertex Collision Region?
#
# 1.2 Filter Track Segments
#
# 1.3 Initialize NEURONS 
#       - randomly?
#
# -------------------------------------------------------------------------------------

# Class Representing a Neuron in the Hopfield Net
# Not sure whether we'll need this
# But maybe we can use it for Weight Matrix creation?

class Neuron:

    # Maybe add some overrides later for update calculation. 
    # For example override MatrixMultiplication / Vector Multiplication / Addition etc.
    # Maybe we have to form an over class: NeuronCollection which we cans use as array.

    def __init__(self, head, tail):
        self.head = head
        self.tail = tail

        self.state = -1 # TODO random initialization?

        #Angle
        #self.angle = head * tail ## Is this correct? [Passaleva] TODO

        #Polar Angles
        # Of the Coordinates or the Lines which are formed by the segment?
        #theta_xz = np.arctan2(y, x) # TODO
        #phi_yz = np.arctan2(y, x) # TODO


## Load an Event
f = open("../events/bsphiphi/velo_event_12.json")
json_data = json.loads(f.read())
event = em.event(json_data)
print(len(event.hits))

## Get the modules
modules = event.modules ## For testing reduce number of used modules

## Go through the modules and create Neurons for all consecutive modules
neurons = []
for i in range(len(modules)-1):
    if i + 2 > len(modules) - 1:
        continue
    print(str(modules[i].module_number) + "->" + str(modules[i+2].module_number))

    combis = 0
    for hit in modules[i]:
        # Go through all hits in modules
        for hit_fl in modules[i+1]:
            # And combine with all hits in following layer to create Neurons
            # TODO Filtering
            neurons.append(Neuron(hit, hit_fl))
            
# Array of Neuron States (We might pass this to the net or the array of Neuron Objects)
neurons_states = [n.state for n in neurons]

print(len(neurons_states))
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

