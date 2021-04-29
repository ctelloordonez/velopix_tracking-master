# Hopfield v2 -> all even mudules
############################### DEPENDENCIES ##################################
import numpy as np
import json
import os
import sys
import numpy as np
import inspect
import os.path
from math import pi, atan, sin, sqrt, tanh

from numpy.core.fromnumeric import shape

########################## LOCAL DEPENDENCIES #################################
filename = inspect.getframeinfo(inspect.currentframe()).filename
file_path = os.path.dirname(os.path.abspath(filename))
project_root = os.path.dirname(file_path)

if project_root not in sys.path:
    sys.path.append(project_root)
from event_model import event_model as em
import data_analysis.event_generator as eg

########################### GLOBAL VARIABLES ##################################

############################### BODY ##########################################


class Hopfield:
    def __init__(self, modules: list, parameters: dict):
        # set self variables, such as the maximum size
        self.p = parameters
        self.m = modules
        self.modules_count = len(modules)
        self.hit_counts = [len(module.hits()) for module in self.m]
        self.neuron_count = [
            self.hit_counts[i] * self.hit_counts[i + 1]
            for i in range(self.modules_count - 1)
        ]
        self.max_neurons = max(self.neuron_count)

        self.init_neurons()
        self.init_weights()

    def init_neurons(self, unit=True):
        self.N = np.ones(shape=(self.modules_count - 1, self.max_neurons))
        for idx, nc in enumerate(self.neuron_count):
            self.N[idx, nc:] = 0

        self.N_info = np.zeros(shape=(self.modules_count - 1, self.max_neurons, 3))
        for idx in range(self.modules_count - 1):
            m1 = self.m[idx]
            m2 = self.m[idx + 1]

            for i, hit1 in enumerate(m1.hits()):
                for j, hit2 in enumerate(m2.hits()):
                    n_idx = i * self.hit_counts[idx + 1] + j
                    # maybe we can check these angles again
                    angle_xz = atan((hit2.x - hit1.x) / (hit2.z - hit1.z))
                    angle_yz = atan((hit2.y - hit1.y) / (hit2.z - hit1.z))
                    norm_dist = sqrt(
                        (hit2.y - hit1.y) ** 2 + (hit2.x - hit1.x) ** 2
                    ) / sqrt((hit2.z - hit1.z) ** 2)
                    self.N_info[idx, n_idx, 0] = angle_xz
                    self.N_info[idx, n_idx, 1] = angle_yz
                    self.N_info[idx, n_idx, 2] = norm_dist

    def init_weights(self):
        #### get params from the dict #######
        alpha = self.p["ALPHA"]
        beta = self.p["BETA"]
        gamma = self.p["GAMMA"]
        #####################################

        self.W = np.zeros(
            shape=(self.modules_count - 2, self.max_neurons, self.max_neurons,)
        )

        # loops neuron neuron weight matrices
        for w_idx in range(self.modules_count - 2):
            # loops hits of the module connecting the neuron layers
            for con_idx in range(self.hit_counts[w_idx + 1]):
                for i in range(self.hit_counts[w_idx]):
                    ln_idx = i * self.hit_counts[w_idx + 1] + con_idx  # left_neuron_idx
                    for j in range(self.hit_counts[w_idx + 2]):
                        # went out of bounds...
                        rn_idx = con_idx * self.hit_counts[w_idx + 2] + j

                        theta = abs(
                            self.N_info[w_idx, ln_idx, 0]
                            - self.N_info[w_idx + 1, rn_idx, 0]
                        )
                        phi = abs(
                            self.N_info[w_idx, ln_idx, 1]
                            - self.N_info[w_idx + 1, rn_idx, 1]
                        )
                        self.W[w_idx, ln_idx, rn_idx] = (
                            alpha
                            * ((1 - sin(theta)) ** beta)
                            * ((1 - sin(phi)) ** gamma)
                        )

    def update(self):
        b = self.p["B"]
        t = self.p["T"]
        # I think here we need to look at the first neurons,
        # then all the neurons in beween as they are dependent on two other layers of neurons
        # then the last layer of the neurons as they only dep on one side

        for idx in range(self.modules_count - 2):
            if idx == 0:
                pass
            if idx == self.modules_count - 1:
                pass
            c1 = self.hit_counts[idx]
            c2 = self.hit_counts[idx + 1]

            # loop over all neurons in the current layer we are looking at
            for i in range(c1 * c2):
                update = 0
                if idx > 0:
                    update += self.W[idx, i, :] @ self.N[idx + 1, :]
                if idx < c1 * c2 - 1:
                    update += self.N[idx, :].T @ self.W[idx, :, i]

                # left module and right module hit id -> current neuron connects hit lm_id with hit rn_id
                lm_id = i // c2
                rm_id = i % c2

                # there can be a lot improved runtime wise with storing the sums and adj
                # but too complicated for now
                # all segments mapping to the hit in m1 -> the left module
                m1h = np.sum(self.N[idx, lm_id * c2 : (lm_id + 1) * c2])

                # all segments mapping to the hit in m2 - the right module
                m2h = np.sum(self.N[idx, : c1 * c2].reshape(c2, c1)[rm_id, :])

                # we need to substact the neuron of the segment 2 times because we add it 2 times
                pen = m1h + m2h - 2 * self.N[idx, i]

                self.N[idx, i] = 0.5 * (1 + tanh(update / (t) - b * pen / (t)))

    def energy(self):
        b = self.p["B"]

        E = 0
        bifurc_pen = 0
        for idx in range(self.modules_count - 2):
            c1 = self.hit_counts[idx]
            c2 = self.hit_counts[idx + 1]
            c3 = self.hit_counts[idx + 2]

            N1_pen = self.N[idx, : c1 * c2].reshape(c2, c1)
            N2_pen = self.N[idx + 1, : c2 * c3].reshape(c2, c3)
            bifurc_pen = (
                np.sum(np.trace(N1_pen @ N1_pen.T))
                + np.sum(np.trace(N2_pen @ N2_pen.T))
                - np.sum(self.N[idx, :] * self.N[idx, :])
                - np.sum(self.N[idx + 1, :] * self.N[idx + 1, :])
            )

            E += (
                -0.5 * (self.N[idx, :].T @ self.W[idx, :, :] @ self.N[idx + 1, :])
                + b * bifurc_pen
            )
        return E

    def converge(self):
        pass

    def network_stats(self):
        #  well this could actually be the __repr__ function of our class
        pass

    def plot_network_results(self):
        pass


def prepare_instance(even=True, num_modules=10, plot_events=False, num_tracks=3):
    if num_modules > 26:
        num_modules = 26
    elif num_modules < 3 or type(num_modules) != int:
        num_modules = 3

    # Generating a test event to work with
    tracks = eg.generate_test_tracks(
        allowed_modules=[i * 2 for i in range(num_modules)],
        num_test_events=1,
        num_tracks=num_tracks,
        reconstructable_tracks=True,
    )[0]

    modules = eg.tracks_to_modules(tracks)
    if plot_events:
        eg.plot_tracks_and_modules(tracks, modules)
    return modules


def load_event():
    # loading some event
    # event_path = os.path.join(project_root, "events/bsphiphi/velo_event_12.json")
    # f = open(event_path)
    # json_data = json.loads(f.read())
    # event = em.event(json_data)
    pass


if __name__ == "__main__":
    #################### PARAMETERS #######################
    ### WHEIGHTS ###
    parameters = {
        "ALPHA": 1,
        "BETA": 10,
        "GAMMA": 10,
        #### UPDATE ###
        "T": 1,
        "B": 1,
        #### Threshold ###
        "maxActivation": False,
        "THRESHOLD": 0.2,
    }
    ###########
    #######################################################

    modules = prepare_instance(num_modules=10, plot_events=True, num_tracks=20)
    my_instance = Hopfield(modules=modules, parameters=parameters)
    for i in range(10):
        print(my_instance.energy())
        my_instance.update()
