# Hopfield v2 -> all even mudules
############################### DEPENDENCIES ##################################
import numpy as np
import json
import os
import sys
import numpy as np
import inspect
import os.path
import matplotlib.pyplot as plt
from math import pi, atan, sin, sqrt, tanh, cosh, exp
from visual.color_map import Colormap
import seaborn as sns

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
        self.N = None
        self.N_info = None
        self.modules_count = len(modules)
        self.hit_counts = [len(module.hits()) for module in self.m]
        self.neuron_count = [
            self.hit_counts[i] * self.hit_counts[i + 1]
            for i in range(self.modules_count - 1)
        ]
        self.max_neurons = max(self.neuron_count)

        self.init_neurons()
        self.init_weights()
        self.extracted_tracks = []
        self.extracted_track_states = []

    def init_neurons(self, unit=True):
        self.N = np.ones(shape=(self.modules_count - 1, self.max_neurons))
        for idx, nc in enumerate(self.neuron_count):
            self.N[idx, nc:] = 1

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
                    self.N_info[idx, n_idx, 0] = abs(angle_xz)
                    self.N_info[idx, n_idx, 1] = abs(angle_yz)
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

                        # Constant term from the other group
                        constant = self.N_info[w_idx, ln_idx, 2] - self.N_info[w_idx + 1, rn_idx, 2]
                        constant = tanh(constant) * (self.p["narrowness"] + 1)  # tanh to force between -1 and 1
                        constant = (-2 * constant**2) + 1  # this should be high if both terms are similar and low/penalizing if both are not similar
                        constant = constant * self.p["constant_factor"]
                        constant = min(max(constant, -self.p["constant_factor"]), self.p["constant_factor"])

                        theta = abs(
                            self.N_info[w_idx, ln_idx, 0]
                            - self.N_info[w_idx + 1, rn_idx, 0]
                        )
                        phi = abs(
                            self.N_info[w_idx, ln_idx, 1]
                            - self.N_info[w_idx + 1, rn_idx, 1]
                        )
                        self.W[w_idx, ln_idx, rn_idx] = max(
                            alpha
                            * ((1 - sin(theta)) ** beta)
                            * ((1 - sin(phi)) ** gamma)
                            + constant, 0)

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
                if idx < self.modules_count - 1:
                    update += self.N[idx, :].T @ self.W[idx, :, i]
                if 0 < idx < self.modules_count - 1:
                    update /= 2

                # left module and right module hit id -> current neuron connects hit lm_id with hit rn_id
                lm_id = i // c2
                rm_id = i % c2

                # there can be a lot improved runtime wise with storing the sums and adj
                # but too complicated for now
                # all segments mapping to the hit in m1 -> the left module
                m1h = np.sum(self.N[idx, lm_id * c2: (lm_id + 1) * c2])

                # all segments mapping to the hit in m2 - the right module
                m2h = np.sum(self.N[idx, : c1 * c2].reshape(c2, c1)[rm_id, :])

                # we need to subtract the neuron of the segment 2 times because we add it 2 times
                pen = m1h + m2h - 2 * self.N[idx, i]

                self.N[idx, i] = 0.5 * (1 + tanh(update / t - b * pen / t))

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
        # Basically keep updating until the difference in Energy between timesteps is lower than 0.0005 (Based on Stimfple-Abele)
        # Passaleva uses a different kind of convergence i think (4)
        energies = [self.energy()] # store all energies (not fastest but maybe nice for visualisations)
        self.update()
        energies.append(self.energy())
        t = 0 # timesteps
        while abs(abs(energies[-2]) - abs(energies[-1])) >= self.p['convergence_threshold']:
            self.update()
            energies.append(self.energy())
            t += 1
            self.p["T"] = self.p["T_decay"](self.p["T"])
        

        print("Network Converged after " + str(t) + " steps")
        print("Energy = " + str(energies[-1]))

    def tracks(self, activation_threshold : list=None):
        if self.extracted_tracks:
            return self.extracted_tracks
        # What the papers say:  The answer is given by the final set of active Neurons
        #                       All sets of Neurons connected together are considered as track candidates
        #                       
        # IDEA: All neurons that share a hit and are both connected are track candidates
        if activation_threshold is None:
            activation_threshold = [0.1, 0.1]
        global_candidates = []
        global_candidate_states = []

        for idx in range(self.modules_count - 2):
            candidates = []
            candidate_states = []
            l1 = self.hit_counts[idx]  # number of hits in module 1
            l2 = self.hit_counts[idx + 1]
            l3 = self.hit_counts[idx + 2]
            # for i, state1 in enumerate(self.N[idx,:]):
            #     for j, state2 in enumerate(self.N[idx+1,:]):
            #         if state1 < activation_threshold[0] or state2 < activation_threshold[1]:
            #             continue
            #         # segments are connected via hit in m2
            #         if int(i % l2) == int(j // l3):
            #             candidates.append(em.track([self.m[idx].hits()[int(i // l2)],
            #                                         self.m[idx + 1].hits()[int(i % l2)],
            #                                         self.m[idx + 2].hits()[int(j % l3)]]))

            if self.p["maxActivation"]:
                candidates = []
                thresh = self.p["THRESHOLD"]
                n1_transform = self.N[idx, :l2*l1].reshape(l2, l1).copy()
                n2_transform = self.N[idx+1, :l3*l2].reshape(l3, l2).copy()
                for con in range(l2):  # loop over the connection hits in module 2
                    h1_idx = np.argmax(n1_transform[con, :])
                    h3_idx = np.argmax(n2_transform[:, con])
                    if n1_transform[con, h1_idx] < thresh or n2_transform[h3_idx, con] < thresh:
                        continue
                    candidates.append(em.track([self.m[idx].hits()[h1_idx],
                                                self.m[idx + 1].hits()[con],
                                                self.m[idx + 2].hits()[h3_idx]]))
                    candidate_states.append(n1_transform[con, h1_idx])
                    candidate_states.append(n2_transform[h3_idx, con])
                    # n1_transform[:, h1_idx] = 0  # set this hit to 0 so it's not chosen again
                    # n2_transform[h3_idx, :] = 0
                    
            global_candidates += candidates
            global_candidate_states += candidate_states

        self.extracted_tracks = global_candidates
        self.extracted_track_states = global_candidate_states

        return global_candidates

    def network_stats(self):
        #  well this could actually be the __repr__ function of our class
        pass

    def plot_network_results(self, show_states=False):
        if show_states:
            # Creates a colormap from blue to red for small to large values respectively
            c_map = Colormap(0, 1, 2 / 3.0, 0)
            colors = []
            [colors.append(c_map.get_color_rgb(v)) for v in self.extracted_track_states]
            print(c_map.get_color_rgb(0))
            eg.plot_tracks_and_modules(self.extracted_tracks, self.m,
                                       colors=colors,
                                       title="Hopfield Output with states")
        else:
            eg.plot_tracks_and_modules(self.extracted_tracks, self.m, title="Hopfield Output")


def prepare_instance(even=True, num_modules=10, plot_events=False, num_tracks=3,
                     save_to_file: str = None):
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

    if save_to_file:
        eg.write_tracks(tracks, save_to_file)

    modules = eg.tracks_to_modules(tracks)
    if plot_events:
        eg.plot_tracks_and_modules(tracks, modules, title="Generated Instance")
    return modules


def load_instance(file_name, plot_events=False):
    tracks = eg.read_tracks(file_name)
    modules = eg.tracks_to_modules(tracks)
    if plot_events:
        eg.plot_tracks_and_modules(tracks, modules, title="Generated Instance")
    return modules


if __name__ == "__main__":
    #################### PARAMETERS #######################
    parameters = {
        ### WEIGHTS ###
        "ALPHA": 1,
        "BETA": 10,
        "GAMMA": 10,
        "narrowness": 100,
        "constant_factor": 1,
        #### UPDATE ###
        "T": 1,
        "B": 0,
        "T_decay": lambda t: max(0.00001, t * 1),
        #### THRESHOLD ###
        "maxActivation": True,
        "THRESHOLD": 0,
        ##### CONVERGENCE ###
        "convergence_threshold": 0.0005
    }
    ###########
    #######################################################

    modules = load_instance("hopfield_test.txt", plot_events=True)
    # modules = prepare_instance(num_modules=10, plot_events=True, num_tracks=2,
    #                            save_to_file="hopfield_test.txt")
    my_instance = Hopfield(modules=modules, parameters=parameters)
    for i in range(len(my_instance.W)):
        sns.heatmap(my_instance.W[i])
        plt.show()
    for i in range(3):
        print("Hello"+"?"*i, my_instance.N_info[:,:,i])
    my_instance.converge()
    print(my_instance.N)

    print("Number of Track elements: " + str(len(my_instance.tracks())))
    print(my_instance.extracted_track_states)


    my_instance.plot_network_results(show_states=True)
