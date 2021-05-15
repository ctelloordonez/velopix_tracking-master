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
import seaborn as sns
from numpy.core.fromnumeric import shape
import random

########################## LOCAL DEPENDENCIES #################################
filename = inspect.getframeinfo(inspect.currentframe()).filename
file_path = os.path.dirname(os.path.abspath(filename))
project_root = os.path.dirname(file_path)

if project_root not in sys.path:
    sys.path.append(project_root)
from event_model import event_model as em
import data_analysis.event_generator as eg
from visual.color_map import Colormap

########################### GLOBAL VARIABLES ##################################

############################### BODY ##########################################


class Hopfield:
    def __init__(self, modules: list, parameters: dict):
        # set self variables, such as the maximum size
        self.p = parameters
        self.m = modules
        self.N = None
        self.N_info = None
        self.sync_rounds = parameters["sync_rounds"]
        self.modules_count = len(modules)
        self.hit_counts = [len(module.hits()) for module in self.m]
        self.neuron_count = [
            self.hit_counts[i] * self.hit_counts[i + 1]
            for i in range(self.modules_count - 1)
        ]
        self.flips = 0
        self.max_neurons = max(self.neuron_count)

        self.init_neurons()
        self.init_weights()
        self.extracted_tracks = []
        self.extracted_track_states = []
        self.energies = []

    def init_neurons(self, unit=True):
        # cosider hits in 2 modules
        # the neurons in N are ordered h1,1-h2,1; h1,1-h2,2; h1,1-h2,3 etc
        if self.p["random_neuron_init"]:
            self.N = np.random.uniform(size=(self.modules_count - 1, self.max_neurons))
        else:
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
            for con_idx in range(self.hit_counts[w_idx + 1]):  # m2
                for i in range(self.hit_counts[w_idx]):  # m1
                    ln_idx = i * self.hit_counts[w_idx + 1] + con_idx  # left_neuron_idx
                    for j in range(self.hit_counts[w_idx + 2]):  # m3
                        rn_idx = con_idx * self.hit_counts[w_idx + 2] + j

                        # Constant term from the other group
                        constant = (
                            self.N_info[w_idx, ln_idx, 2]
                            - self.N_info[w_idx + 1, rn_idx, 2]
                        )
                        constant = tanh(constant) * (
                            self.p["narrowness"] + 1
                        )  # tanh to force between -1 and 1
                        constant = (
                            -2 * constant ** 2
                        ) + 1  # this should be high if both terms are similar and low/penalizing if both are not similar
                        constant = constant * self.p["constant_factor"]
                        constant = min(
                            max(constant, -self.p["constant_factor"]),
                            self.p["constant_factor"],
                        )

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
                            + constant,
                            0,
                        )

    def update(self, synchronous=False):
        b = self.p["B"]
        t = self.p["T"]
        # I think here we need to look at the first neurons,
        # then all the neurons in between as they are dependent on two other layers of neurons
        # then the last layer of the neurons as they only dep on one side

        N_sync = None
        if synchronous:
            N_sync = self.N.copy()

        update_list = []
        for idx in range(self.modules_count - 1):
            c1 = self.hit_counts[idx]
            c2 = self.hit_counts[idx + 1]
            for i in range(c1 * c2):
                update_list.append((idx, i))

        if self.p["randomized_updates"]:
            random.shuffle(update_list)

        for idx, i in update_list:
            c1 = self.hit_counts[idx]
            c2 = self.hit_counts[idx + 1]

            update = 0
            if idx > 0:
                update += self.N[idx - 1, :].T @ self.W[idx - 1, :, i]

            if idx < self.modules_count - 2:
                update += self.W[idx, i, :] @ self.N[idx + 1, :]

            if 0 < idx < self.modules_count - 1:
                update /= 2

            # left module and right module hit id -> current neuron connects hit lm_id with hit rn_id
            lm_id = i // c2
            rm_id = i % c2

            # there can be a lot improved runtime wise with storing the sums and adj
            # but too complicated for now
            # all segments mapping to the hit in m1 -> the left module
            m1h = np.sum(self.N[idx, lm_id * c2 : (lm_id + 1) * c2])

            # all segments mapping to the hit in m2 - the right module
            m2h = np.sum(
                self.N[idx, : c1 * c2].reshape(c2, c1)[rm_id, :]
            )  # correct as well...

            # we need to subtract the neuron of the segment 2 times because we add it 2 times
            pen = m1h + m2h - 2 * self.N[idx, i]

            if synchronous:
                N_sync[idx, i] = 0.5 * (1 + tanh(update / t - b * pen / t))
            else:
                _update = 0.5 * (1 + tanh(update / t - b * pen / t))
                if self.p["binary_states"]:
                    if random.random() < _update:
                        self.N[idx, i] = 1
                    else:
                        self.N[idx, i] = 0
                else:
                    self.N[idx, i] = _update
                # if np.random.random() < t:
                #    self.flips += 1
                #    self.N[idx, i] = 1 - self.N[idx, i]
            if idx == 2 and i > 3:
                pass
        if synchronous:
            self.N = N_sync

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
        self.energies = [
            self.energy()
        ]  # store all energies (not fastest but maybe nice for visualisations)
        t = 0  # timesteps
        # print(f"N at iteration{t}:", np.round(my_instance.N, 1))
        self.update(synchronous=t < self.p["sync_rounds"])
        t += 1
        self.energies.append(self.energy())
        while (
            abs(abs(self.energies[-2]) - abs(self.energies[-1]))
            >= self.p["convergence_threshold"]
        ):
            self.update(synchronous=t < self.p["sync_rounds"])
            self.energies.append(self.energy())
            # print(f"N at iteration{t}:", np.round(my_instance.N, 1))
            t += 1
            self.p["T"] = self.p["T_decay"](self.p["T"])
            # XXX: added b decay
            self.p["B"] = self.p["B_decay"](t)

            # print("T: " + str(t) + " Flips : " + str(self.flips))

        print("Network Converged after " + str(t) + " steps")
        print("Energy = " + str(self.energies[-1]))
        return self.N

    def bootstrap_converge(self, bootstraps=50):
        states_list = []
        for i in range(bootstraps):

            if self.p["random_neuron_init"]:
                # We only need to reinitialize if we randomly initialize
                self.init_neurons()

            states = self.converge()

            states_list.append(states)

        stacked_states = np.stack(states_list, axis=2)

        self.N = np.mean(stacked_states, axis=2)

    def tracks(self, activation_threshold: list = None):
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

            if self.p["maxActivation"]:
                candidates = []
                thresh = self.p["THRESHOLD"]
                if l2 != l1:
                    pass
                print(self.N[idx, : l2 * l1])
                n1_transform = self.N[idx, : l2 * l1].reshape(l1, l2).T.copy()
                n2_transform = self.N[idx + 1, : l3 * l2].reshape(l2, l3).T.copy()
                for con in range(l2):  # loop over the connection hits in module 2
                    h1_idx = np.argmax(n1_transform[con, :])
                    h3_idx = np.argmax(n2_transform[:, con])
                    if (
                        n1_transform[con, h1_idx] < thresh
                        or n2_transform[h3_idx, con] < thresh
                    ):
                        continue
                    candidates.append(
                        em.track(
                            [
                                self.m[idx].hits()[h1_idx],
                                self.m[idx + 1].hits()[con],
                                self.m[idx + 2].hits()[h3_idx],
                            ]
                        )
                    )
                    candidate_states.append(n1_transform[con, h1_idx])
                    candidate_states.append(n2_transform[h3_idx, con])
                    n1_transform[
                        :, h1_idx
                    ] = 0  # set this hit to 0 so it's not chosen again
                    n2_transform[h3_idx, :] = 0

            global_candidates += candidates
            global_candidate_states += candidate_states

        self.extracted_tracks = global_candidates
        self.extracted_track_states = global_candidate_states

        return global_candidates

    def show_all_tracks(self, colors=False):
        # Creates a colormap from blue to red for small to large values respectively
        c_map = Colormap(0, 1, 2 / 3.0, 0)
        c = []
        tracks = []
        for idx in range(self.modules_count - 1):
            m1 = self.m[idx]
            m2 = self.m[idx + 1]

            for i, hit1 in enumerate(m1.hits()):
                for j, hit2 in enumerate(m2.hits()):
                    if colors:
                        n_idx = i * self.hit_counts[idx + 1] + j
                        c.append(c_map.get_color_rgb(self.N[idx, n_idx]))
                    tracks.append(em.track([hit1, hit2]))
        eg.plot_tracks_and_modules(tracks, self.m, colors=c)

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
            eg.plot_tracks_and_modules(
                self.extracted_tracks,
                self.m,
                colors=colors,
                title="Hopfield Output with states",
            )
        else:
            eg.plot_tracks_and_modules(
                self.extracted_tracks, self.m, title="Hopfield Output"
            )


def prepare_instance(
    even=True, num_modules=10, plot_events=False, num_tracks=3, save_to_file: str = None
):
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
    for t in tracks:
        print([hit.y for hit in t.hits])
        print()
    modules = eg.tracks_to_modules(tracks)
    if plot_events:
        eg.plot_tracks_and_modules(tracks, modules, title="Generated Instance")
    return modules


if __name__ == "__main__":
    #################### PARAMETERS #######################
    parameters = {
        ### NEURONS ###
        "random_neuron_init": True,
        "binary_states": False,
        ### WEIGHTS ###
        "ALPHA": 1,
        "BETA": 10,
        "GAMMA": 10,
        "narrowness": 200,
        "constant_factor": 0.9,
        #### UPDATE ###
        "T": 5,
        "B": 0.1,
        "B_decay": lambda t: max(0.1, t * 0.01),
        "T_decay": lambda t: max(0.00001, t * 0.8),
        "sync_rounds": 0,
        "randomized_updates": True,
        #### THRESHOLD ###
        "maxActivation": True,
        "THRESHOLD": 0.2,
        ##### CONVERGENCE ###
        "convergence_threshold": 0.00005,
    }
    ###########
    #######################################################

    # modules = load_instance("test.txt", plot_events=False)
    modules = prepare_instance(
        num_modules=3, plot_events=True, num_tracks=3, save_to_file="test.txt"
    )
    for m in modules:
        m.hits()
        print([hit.y for hit in m.hits()])

    my_instance = Hopfield(modules=modules, parameters=parameters)
    # for i in range(len(my_instance.W)):
    #     sns.heatmap(my_instance.W[i])
    #     plt.show()
    # for i in range(3):
    #     print(my_instance.N_info[:,:,i])

    print(np.round(my_instance.N, 2))
    my_instance.converge()

    print("Number of Track elements: " + str(len(my_instance.tracks())))

    # plt.plot(my_instance.energies)
    # plt.show()

    print(np.round(my_instance.N, 1))

    print(np.shape(my_instance.N))

    my_instance.bootstrap_converge(bootstraps=50)
    print(my_instance.flips)
    my_instance.plot_network_results(show_states=True)

    my_instance.show_all_tracks(colors=True)

