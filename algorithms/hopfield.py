import numpy as np

#### Interesting Implementations
# https://github.com/andreasfelix/hopfieldnetwork/tree/9e06c43c71c858afe4d8fc74f4c11c63ef83c22c



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

