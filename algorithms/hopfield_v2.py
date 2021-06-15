# Hopfield v2 -> all even mudules
############################### DEPENDENCIES ##################################
import json
import os
import sys
import contextlib
import io
import numpy as np
import inspect
import os.path
import matplotlib.pyplot as plt
from math import pi, atan, sin, sqrt, tanh, cosh, exp, ceil
import seaborn as sns
from numpy.core.fromnumeric import shape
import random
import time
import math

########################## LOCAL DEPENDENCIES #################################
filename = inspect.getframeinfo(inspect.currentframe()).filename
file_path = os.path.dirname(os.path.abspath(filename))
project_root = os.path.dirname(file_path)

if project_root not in sys.path:
    sys.path.append(project_root)
from event_model import event_model as em
from validator import validator_lite as vl
import data_analysis.event_generator as eg
from visual.color_map import Colormap

########################### GLOBAL VARIABLES ##################################


########################### CONTEXTS ##################################
@contextlib.contextmanager
def nostdout():
    save_stdout = sys.stdout
    sys.stdout = io.BytesIO()
    yield
    sys.stdout = save_stdout


########################### HELPER FUNCTIONS ##################################
def get_polar_coordinates(x, y):
    r = math.sqrt(x ** 2 + y ** 2)
    phi = math.atan2(x, y)
    if phi < 0:
        phi = math.pi - phi
    return r, phi


############################### BODY ##########################################


### DESIRED INPUT for VALIDATOR###
# what the validation function needs
# json_data_all_events => [json_data_event_1, json_data_event_2, ..., json_data_event_n]
# all_tracks => [all_tracks_event_1, all_tracks_event_2, ..., all_tracks_event_n]

# all_tracks_event_j => a list of all tracks in event j. Tracks are track objects.
###


class Hopfield:
    def __init__(self, modules: list, parameters: dict, tracks: list = None):
        # set self variables, such as the maximum size
        self.p = parameters
        self.m = modules
        self.start_T = self.p["T"]
        self.start_B = self.p["B"]
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

        self.init_neurons(tracks=tracks)
        self.init_weights()
        self.extracted_hits = set()
        self.extracted_tracks = []
        self.extracted_track_states = []
        self.energies = []

    def init_neurons(self, tracks: list = None):
        # cosider hits in 2 modules
        # the neurons in N are ordered h1,1-h2,1; h1,1-h2,2; h1,1-h2,3 etc
        if self.p["random_neuron_init"]:
            self.N = np.random.uniform(size=(self.modules_count - 1, self.max_neurons))
        else:
            self.N = np.ones(shape=(self.modules_count - 1, self.max_neurons))
        if tracks:
            self.N = np.zeros(shape=(self.modules_count - 1, self.max_neurons))
        for idx, nc in enumerate(self.neuron_count):
            self.N[idx, nc:] = 0
        self.N_info = np.zeros(shape=(self.modules_count - 1, self.max_neurons, 4))
        for idx in range(self.modules_count - 1):
            m1 = self.m[idx]
            m2 = self.m[idx + 1]

            for i, hit1 in enumerate(m1.hits()):
                for j, hit2 in enumerate(m2.hits()):
                    n_idx = i * self.hit_counts[idx + 1] + j
                    if tracks:
                        for t in tracks:
                            if hit1 in t and hit2 in t:
                                self.N[idx, n_idx] = 1
                    # maybe we can check these angles again
                    angle_xz = atan((hit2.x - hit1.x) / (hit2.z - hit1.z))
                    angle_yz = atan((hit2.y - hit1.y) / (hit2.z - hit1.z))
                    norm_dist = sqrt(
                        (hit2.y - hit1.y) ** 2 + (hit2.x - hit1.x) ** 2
                    ) / sqrt((hit2.z - hit1.z) ** 2)

                    _, r_hit1 = get_polar_coordinates(hit1.x, hit1.y)
                    _, r_hit2 = get_polar_coordinates(hit2.x, hit2.y)
                    monotone_dist = (r_hit2 - r_hit1) / (hit2.z - hit1.z)

                    self.N_info[idx, n_idx, 0] = abs(angle_xz)
                    self.N_info[idx, n_idx, 1] = abs(angle_yz)
                    self.N_info[idx, n_idx, 2] = norm_dist
                    self.N_info[idx, n_idx, 3] = monotone_dist

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

                        # monotone constant
                        monotone_constant = (
                            self.N_info[w_idx, ln_idx, 3]
                            - self.N_info[w_idx + 1, rn_idx, 3]
                        )
                        monotone_constant = tanh(monotone_constant) * (
                            self.p["narrowness"] + 1
                        )  # tanh to force between -1 and 1
                        monotone_constant = (
                            -2 * monotone_constant ** 2
                        ) + 1  # this should be high if both terms are similar and low/penalizing if both are not similar
                        monotone_constant = (
                            monotone_constant * self.p["monotone_constant_factor"]
                        )
                        monotone_constant = min(
                            max(monotone_constant, -self.p["monotone_constant_factor"]),
                            self.p["monotone_constant_factor"],
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
                            + constant
                            + monotone_constant,
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

            # XXX: just realised we are counting most penalties too often!! -> fixed
            f1 = 0.5
            f2 = 0.5
            if idx == 0:
                f1 = 1
            if idx == self.modules_count - 3:
                f2 = 1
            N1_pen = self.N[idx, : c1 * c2].reshape(c2, c1)
            N2_pen = self.N[idx + 1, : c2 * c3].reshape(c2, c3)
            bifurc_pen = (
                np.sum(np.trace(N1_pen @ N1_pen.T)) * f1
                + np.sum(np.trace(N2_pen @ N2_pen.T)) * f2
                - np.sum(self.N[idx, :] * self.N[idx, :]) * f1
                - np.sum(self.N[idx + 1, :] * self.N[idx + 1, :]) * f2
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
        self.p["T"] = self.start_T
        # self.p["B"] = self.start_B
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
            if not self.p["decay_off"]:
                self.p["T"] = self.p["T_decay"](self.p["T"])
                # XXX: added b decay
                self.p["B"] = self.p["B_decay"](t)
            else:
                pass

            # print("T: " + str(t) + " Flips : " + str(self.flips))

        # print("Network Converged after " + str(t) + " steps")
        # print("Energy = " + str(self.energies[-1]))
        return self.N, self.energies[-1]

    def bootstrap_converge(self, bootstraps=50, method="mean"):
        start_time = time.time()
        states_list = []
        energy_list = []
        for i in range(bootstraps):

            if self.p["random_neuron_init"]:
                # We only need to reinitialize if we randomly initialize
                self.init_neurons()

            states, energy = self.converge()
            print("energy: " + str(energy))

            states_list.append(states)
            energy_list.append(energy)
            # print(f"Finished {i+1}/{bootstraps} iterations")

        if method == "minimum":
            self.N = states_list[np.argmax(energy_list)]
            energy_list = [np.amax(energy_list)]
        else:
            stacked_states = np.stack(states_list, axis=2)
            self.N = np.mean(stacked_states, axis=2)

        end_time = time.time() - start_time
        print(
            "[HOPFIELD] converged network by %s after %i mins %.2f seconds; (energy: %.2f)"
            % (method, end_time // 60, end_time % 60, np.mean(energy_list))
        )

    def tracks(self):
        # What the papers say:  The answer is given by the final set of active Neurons
        #                       All sets of Neurons connected together are considered as track candidates
        #
        # IDEA: All neurons that share a hit and are both connected are track candidates
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

                n1_transform = self.N[idx, : l2 * l1].reshape(l1, l2).T.copy()
                n2_transform = self.N[idx + 1, : l3 * l2].reshape(l2, l3).T.copy()

                for con in range(l2):  # loop over the connection hits in module 2
                    # XXX i try swapping these....
                    h1_idx = np.argmax(n1_transform[con, :])
                    h3_idx = np.argmax(n2_transform[:, con])

                    if (
                        n1_transform[con, h1_idx] < thresh
                        or n2_transform[h3_idx, con] < thresh
                    ):
                        continue

                    hit1 = self.m[idx].hits()[h1_idx]
                    hit2 = self.m[idx + 1].hits()[con]
                    hit3 = self.m[idx + 2].hits()[h3_idx]
                    candidates.append(em.track([hit1, hit2, hit3]))
                    self.extracted_hits.add(hit1)
                    self.extracted_hits.add(hit2)
                    self.extracted_hits.add(hit3)
                    # if we get the same state
                    candidate_states.append(n1_transform[con, h1_idx])
                    candidate_states.append(n2_transform[h3_idx, con])

                    # XXX
                    # this prevents the display of bifurcation?!
                    # n1_transform[
                    #    :, h1_idx
                    # ] = 0  # set this hit to 0 so it's not chosen again
                    # n2_transform[h3_idx, :] = 0

            global_candidates += candidates
            global_candidate_states += candidate_states

        self.extracted_tracks = global_candidates
        # self.process_triplet_tracks()
        self.extracted_track_states = global_candidate_states

        return global_candidates

    def full_tracks(self):
        # this will deal with stange angles!!!
        global_candidates = []
        global_candidate_states = []
        global_candidate_info = []
        # under the assumption that we removed bifuration completely
        # init this active tracks with all active neurons in layer 1! -> key is the right hit
        tracks = {}

        for idx in range(self.modules_count - 1):
            tracks_2 = {}
            l1 = self.hit_counts[idx]  # number of hits in module 1
            l2 = self.hit_counts[idx + 1]

            tr = self.p["THRESHOLD"]
            for segment in range(l1 * l2):
                if self.N[idx, segment] < tr:
                    continue

                r_hit = self.m[idx + 1].hits()[segment % l2]
                l_hit = self.m[idx].hits()[segment // l2]

                if l_hit in tracks.keys():
                    (track, states, angle) = tracks[l_hit]

                    track = track + [r_hit]
                    states = states + [self.N[idx, segment]]
                    info = angle + [self.N_info[idx, segment, :]]
                    del tracks[l_hit]

                    self.extracted_hits.add(r_hit)
                    tracks_2[r_hit] = (track, states, info)

                else:
                    track = [l_hit, r_hit]
                    states = [self.N[idx, segment]]
                    info = [self.N_info[idx, segment, :]]

                    tracks_2[r_hit] = (track, states, info)
                    self.extracted_hits.add(r_hit)
                    self.extracted_hits.add(l_hit)

            for _, value in tracks.items():
                (track, states, info) = value
                global_candidates = global_candidates + [track]
                global_candidate_states = global_candidate_states + [states]
                global_candidate_info = global_candidate_info + [info]
            tracks = tracks_2

        for _, value in tracks.items():
            (track, states, info) = value
            global_candidates = global_candidates + [track]
            global_candidate_states = global_candidate_states + [states]
            global_candidate_info = global_candidate_info + [info]

        # here comes the function of 'pruning...' maybe i need to store more info for doing that!!!
        global_candidates = self.prune_tracks(global_candidates, global_candidate_info)

        global_candidates = [em.track(hits) for hits in global_candidates]
        self.extracted_tracks = global_candidates
        self.extracted_track_states = global_candidate_states

        return global_candidates

    # tr 0.1 seems decent for sum of info differences...
    # we could look more carully into a criterion for this on big instances but here is fine...
    def prune_tracks(self, tracks, track_infos, tr=0.01):
        out_tracks = []
        for track, info in zip(tracks, track_infos):
            num_hits = len(track)
            if num_hits < 3:  # sorting out the tracks that are not relevant
                continue
            # only if len> 6 need to think about splitting!!!
            cand = [track[0], track[1]]
            cand_info = info[0]
            for idx in range(1, num_hits - 1):
                # print(abs(sum(cand_info - info[idx])))
                if abs(sum(cand_info - info[idx])) < tr:
                    cand = cand + [track[idx + 1]]
                else:
                    if len(cand) > 2:
                        out_tracks = out_tracks + [cand]
                    cand = [track[idx], track[idx + 1]]
                cand_info = info[idx]

            if len(cand) > 2:
                out_tracks = out_tracks + [cand]
        return out_tracks

    def mark_bifurcation(self, zero, max_activation):
        if max_activation:
            zero = False
        tr = self.p["THRESHOLD"]
        self.N[self.N <= tr] = 0

        # search for bifurcation neurons

        for idx in range(self.modules_count - 1):
            # so basically we visit all neurons in one layer and check for neurons where the activation is bigger than tr
            # then we check all adjacent neurons for activation and check how many exceid the treshold
            # for each segment we look wether there is bifurcation on the left or right hit

            c1 = self.hit_counts[idx]
            c2 = self.hit_counts[idx + 1]

            for segment in range(c1 * c2):
                if self.N[idx, segment] < tr:
                    continue

                r_hit = segment % c2
                l_hit = segment // c2

                # left-right bifurction
                activation_mask = self.N[idx, : c1 * c2].reshape(c1, c2)[:, r_hit] > tr
                if sum(activation_mask) > 1:  # we have bifuct into the right hit
                    affected_neurons = []
                    for i in range(c1):  # loop over all nerons affected by bifurc
                        if activation_mask[i]:
                            # well here are the bifurcation things detected. here we would need to come up with a smart way to resolve it
                            if zero:
                                self.N[idx, (i * c2) + r_hit] = 0
                            else:
                                affected_neurons = affected_neurons + [(i * c2) + r_hit]

                    if max_activation:
                        # check to the left or max score
                        max_activation = self.N[idx, affected_neurons[0]]
                        max_id = 0
                        for e in affected_neurons:
                            if self.N[idx, e] > max_activation:
                                max_id = e
                                max_activation = self.N[idx, e]

                            self.N[idx, e] = 0

                        self.N[idx, max_id] = 1
                        # simple rule -> when bifurc is detected on right side -> we look next active neurons going out
                        # and promote the ones where the weight is high... (angle diff is low)
                        # if next neuron layer exist!!!
                        # also clean bifurc there if active

                    # print(
                    #    f"""neuron_layer {idx}: detected bifurc in right_hit {r_hit}, y: {self.m[idx+1].hits()[r_hit].y}
                    #    from y: {self.m[idx].hits()[l_hit].y}, activation: {self.N[idx, segment]}"""
                    # )
                    # i set all affected neurons to 1 and all others to 0 to see if it works
                    # deceide here which one survives (higher activation)
                # right-left bifurcation
                activation_mask = self.N[idx, : c1 * c2].reshape(c1, c2)[l_hit, :] > tr
                if sum(activation_mask) > 1:
                    affected_neurons = []
                    for i in range(c2):
                        if activation_mask[i]:
                            if zero:
                                self.N[idx, (l_hit * c2) + i] = 0
                            else:
                                affected_neurons = affected_neurons + [(l_hit * c2) + i]

                    if max_activation:
                        # check to the left or max score
                        max_activation = self.N[idx, affected_neurons[0]]
                        max_id = 0
                        for e in affected_neurons:
                            if self.N[idx, e] > max_activation:
                                max_id = e
                                max_activation = self.N[idx, e]
                            self.N[idx, e] = 0
                        self.N[idx, max_id] = 1

        # converged, averaged neuron state
        # what do we want to do -> search all neurons wether there is bifurcation
        # how to search this, by the indices... and then we store it as a combination of hit id_s, and which side of this element the bifurcation occurs
        # bifiurcation is stored as a list of hit ids, with

    def process_triplet_tracks(self):
        total_hits = sorted(self.extracted_hits, key=lambda hit: -hit.y)
        total_hits = sorted(total_hits, key=lambda hit: hit.module_number)
        print(len(total_hits))
        total_tracks = []
        while len(total_hits) != 0:
            c_hit = total_hits[0]
            c_track = [c_hit]
            total_hits.remove(c_hit)
            while c_hit is not None:
                candidate_tracks = self.tracks_with_hit(c_hit)
                candidate_hits = set()
                for track in candidate_tracks:
                    index = track.hits.index(c_hit)
                    if index < 2:
                        candidate_hits.add(track.hits[index + 1])
                if len(candidate_hits) == 0:
                    break
                c_hit = list(candidate_hits)[
                    0
                ]  # Should really use the weight matrix...
                if c_hit not in total_hits:
                    break
                c_track.append(c_hit)
                total_hits.remove(c_hit)
            total_tracks.append(em.track(c_track))
        eg.plot_tracks_and_modules(total_tracks, self.m)

    def show_all_tracks(self, threshold=None, show_states=False):
        # Creates a colormap from blue to red for small to large values respectively
        c_map = Colormap(0, 1, 2 / 3.0, 0)
        c = []
        tracks = []
        for idx in range(self.modules_count - 1):
            m1 = self.m[idx]
            m2 = self.m[idx + 1]

            for i, hit1 in enumerate(m1.hits()):
                for j, hit2 in enumerate(m2.hits()):
                    n_idx = i * self.hit_counts[idx + 1] + j
                    if threshold:
                        if self.N[idx, n_idx] >= threshold:
                            tracks.append(em.track([hit1, hit2]))
                            if show_states:
                                c.append(c_map.get_color_rgb(self.N[idx, n_idx]))
                        continue
                    if show_states:
                        c.append(c_map.get_color_rgb(self.N[idx, n_idx]))
                    tracks.append(em.track([hit1, hit2]))
        eg.plot_tracks_and_modules(tracks, self.m, colors=c)

    def tracks_with_hit(self, hit):
        return [track for track in self.extracted_tracks if hit in track.hits]

    def print_neurons(self):
        n = len(self.N)
        for i in range(n):
            m = int(sqrt(len(self.N[i])))
            for j in range(m):
                print(f"m{i+1}h{j+1}: {self.N[i, (j*m):((j+1)*m)]}")

    def plot_network_results(self, show_states=False):
        if show_states:
            # Creates a colormap from blue to red for small to large values respectively
            c_map = Colormap(0, 1, 2 / 3.0, 0)
            colors = []
            [colors.append(c_map.get_color_rgb(v)) for v in self.extracted_track_states]
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
    return modules, tracks


def load_instance(file_name, plot_events=False):
    tracks = eg.read_tracks(file_name)
    modules = eg.tracks_to_modules(tracks)
    if plot_events:
        eg.plot_tracks_and_modules(tracks, modules, title="Generated Instance")
    return modules, tracks


def load_event(file_name, plot_event=False):
    f = open(file_name)
    json_data_event = json.loads(f.read())

    ev = em.event(json_data_event, read_tracks=True)

    modules = ev.modules
    tracks = ev.real_tracks

    if plot_event:
        eg.plot_tracks_and_modules(tracks, modules, title="Loaded Event")

    modules_even = []
    modules_odd = []

    for i in range(len(modules)):
        if i % 2 == 0:
            modules_even.append(modules[i])
        else:
            modules_odd.append(modules[i])

    return json_data_event, (modules_even, modules_odd)


def evaluate_events(
    file_name, parameters, nr_events=1, plot_event=False, bootstrap_method="mean"
):

    json_data_all_events = []
    all_tracks = []

    for i in range(nr_events):

        print("[INFO] Evaluate Event: %s" % file_name + str(i))
        json_data_event, modules = load_event(
            file_name + str(i) + ".json", plot_event=True
        )

        start_time = time.time()
        even_hopfield = Hopfield(modules=modules[0], parameters=parameters)
        odd_hopfield = Hopfield(modules=modules[1], parameters=parameters)
        end_time = time.time() - start_time
        print(
            "[INFO] Hopfield Networks initialized in %i mins %.2f seconds"
            % (end_time // 60, end_time % 60)
        )

        even_hopfield.bootstrap_converge(bootstraps=4, method=bootstrap_method)
        odd_hopfield.bootstrap_converge(bootstraps=4, method=bootstrap_method)

        start_time = time.time()
        even_hopfield.mark_bifurcation(zero=False, max_activation=True)
        odd_hopfield.mark_bifurcation(zero=False, max_activation=True)
        even_tracks = even_hopfield.full_tracks()
        odd_tracks = odd_hopfield.full_tracks()
        event_tracks = even_tracks + odd_tracks
        end_time = time.time() - start_time
        print(
            "[INFO] tracks extracted in %i mins %.2f seconds"
            % (end_time // 60, end_time % 60)
        )

        json_data_all_events.append(json_data_event)
        all_tracks.append(event_tracks)

        if plot_event:
            even_hopfield.plot_network_results()
            odd_hopfield.plot_network_results()

    start_time = time.time()
    velo, long, long5, long_strange, \
    long_strange5, long_fromb, long_fromb5 = vl.validate_print(json_data_all_events, all_tracks, return_data=True)
    end_time = time.time() - start_time

    # we could check how many tracks acutally cross the detector sides i guess to identify where some clones come from...
    print(
        "[INFO] validation excecuted in %i mins %.2f seconds"
        % (end_time // 60, end_time % 60)
    )


def mse(network, tracks):
    true_network = Hopfield(modules=network.m, parameters=network.p, tracks=tracks)
    return ((network.N - true_network.N) ** 2).mean(axis=None)


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
        "monotone_constant_factor": 0.0,
        #### UPDATE ###
        "T": 5,
        "B": 0.2,
        "B_decay": lambda t: max(0.1, t * 0.04),
        "T_decay": lambda t: max(1e-8, t * 0.01),
        "decay_off": False,
        "sync_rounds": 0,
        "randomized_updates": True,
        #### THRESHOLD ###
        "maxActivation": True,
        "THRESHOLD": 0.3,
        ##### CONVERGENCE ###
        "convergence_threshold": 0.00000005,
    }
    ###########
    #######################################################

    # modules, tracks = load_instance("test.txt", plot_events=True)
    # modules, tracks = prepare_instance(
    #    num_modules=26, plot_events=True, num_tracks=30, save_to_file="test.txt"
    # )

    evaluate_events(
        project_root + "/events/small_dataset/velo_event_",
        parameters,
        nr_events=1,
        plot_event=True,
    )
    exit()
    my_instance = Hopfield(modules=modules, parameters=parameters)
    # for i in range(len(my_instance.W)):
    #     sns.heatmap(my_instance.W[i])
    #     plt.show()
    # for i in range(3):
    #     print(my_instance.N_info[:,:,i])

    # print(np.round(my_instance.N, 2))
    np.random.seed(20)
    my_instance.converge()

    # plt.plot(my_instance.energies)
    # plt.show()

    my_instance.bootstrap_converge(bootstraps=1, method="mean")
    # print("Converged:", my_instance.energies[-1])
    # my_instance.tracks()
    # my_instance.plot_network_results(show_states=True)
    my_instance.mark_bifurcation(zero=False, max_activation=True)
    my_instance.full_tracks()
    my_instance.plot_network_results(show_states=False)
    # my_instance.tracks()
    # my_instance.plot_network_results(show_states=True)

    # true_instance = Hopfield(modules=modules, parameters=parameters, tracks=tracks)
    # print("True:", true_instance.energy())
    # true_instance.tracks()
    # true_instance.plot_network_results()

    # THIS CODE IS FOR COMPUTING THE MSES FOR DIFFERENT BOOTSTRAPS

    # total_mse = []
    # for i in range(50):
    #     mses = []
    #     for j in range(ceil(50/(i+1))):
    #         instance = Hopfield(modules=modules, parameters=parameters)
    #         instance.bootstrap_converge(bootstraps=i+1)
    #         mses.append(mse(instance, tracks))
    #     total_mse.append(sum(mses)/len(mses))
    #
    # plt.plot(total_mse)
    # plt.show()

