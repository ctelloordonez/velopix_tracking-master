# get the relevant dependencies 
import json
import os
import sys 
module_path = os.path.abspath(os.path.join('..'))
if module_path not in sys.path:
    sys.path.append(module_path)
from event_model import event_model as em


import event_model as em

#### Interesting Implementations
# https://github.com/andreasfelix/hopfieldnetwork/tree/9e06c43c71c858afe4d8fc74f4c11c63ef83c22c

class smallHoppfiedNetwork():
    # takes 3 modules
    def __init__(self, m1, m2, m3) -> None:
        N1 = None # all connections between m1 & m2 -> dtype np.array size |m1| * |m2|  1000
        N2 = None # all connections btw m2, m3 -> dtype np.array size |m2| * |m3| 1000
        
        N1_info = None
        N2_info = None
        
        weights = None # dim |N1| * |N2| 

    # fill the weight matrix based on the geometric properties of the hits / segments
    def init_weights(self):
        pass
    
    def energy(self):
        pass

    def update(self):
        pass

    def wibble_wabble(self):
        pass

    def tracks(self):
        pass


if __name__ == "__main__":
    f = open("./events/bsphiphi/velo_event_12.json")
    json_data = json.loads(f.read())
    event = event.em.event(json_data)
    print(event)
    pass



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

