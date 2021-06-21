import math
import numpy as np
from event_model import event_model as em


class SearchByPhi:
    def __init__(self, event):
        self.hits = sort_by_phi(event.hits) # sort hits by phi
        self.acceptance = 0.015
        self.min_track_length = 3

    def solve(self):
        tracks = [] # list of tracks found
        current_track = [] # list of hits representing the track currently under consideration
        current_track.append(self.hits[0]) # initialize by considering the first hit
        self.hits = self.hits[1:] # Take the first hit out of consideration for forming a track with itself
        skipped = 0 # variable to keep track of the number of modules skipped
        for h in self.hits: # loop over the hits in the event
            if abs(math.atan2(h.y, h.x) - math.atan2(current_track[0].y, current_track[0].x)) < self.acceptance: # if the hits polar angle is withing range of the track
                # We skipped 1 module in the track already, hit is in module after the first hit of the track
                if skipped == 1 and h.module_number == current_track[0].module_number + 2:
                    current_track.insert(0, h) # add the hit after the first track hit
                # We skipped 1 module in the track already, hit is in module before the first hit of the track
                elif skipped == 1 and h.module_number == current_track[-1].module_number - 2:
                    current_track.append(h) # add the hit at the end of the track
                # We did not yet skip a module in the current track,  hit is in module before the first hit of the track
                elif h.module_number == current_track[0].module_number - 2:
                    current_track.insert(0, h) # add the hit after the first track hit
                # We did not yet skip a module in the current track,  hit is in module after the last hit of the track
                elif h.module_number == current_track[-1].module_number + 2:
                    current_track.append(h) # add the hit at the end of the track
                # We did not yet skip a module in the current track,  hit is 2 modules before the first hit of the track
                elif skipped <= 0 and h.module_number == current_track[0].module_number - 4:
                    current_track.insert(0, h) # add the hit
                    skipped += 1 # increase skipped module counter
                # We did not yet skip a module in the current track,  hit is 2 modules after the last hit of the track
                elif skipped <= 0 and h.module_number == current_track[-1].module_number + 4:
                    current_track.append(h) # add the hit, skip next module counter
                    skipped += 1 # increase skipped module counter
                else:
                    skipped = 0 # reset skipped modules counnter
                    if len(current_track) > 1: # if the current track has a minimal length
                        tracks.append(em.track(current_track)) # add it to the list of tracks found
                    current_track = [] # Reset current track
                    current_track.append(h) # Add the hit as the initial hit of the track
            else: # if hit h is not sufficiently close to the current track
                skipped = 0 # reset skipped counter
                if len(current_track) >= self.min_track_length: # If the current track has a minimal number of hits
                    tracks.append(em.track(current_track)) # add it to the list of tracks found
                current_track = [] # Reset current track
                current_track.append(h) # Add the hit as the initial hit of the track
        if len(current_track) >= self.min_track_length:
            tracks.append(em.track(current_track)) # append final track if longer than 1

        return tracks

    def solve2(self):
        grouped_by_phi = []
        g_track = []
        for h in self.hits:
            if len(g_track) <= 0:
                g_track.append(h)
            elif abs(math.atan2(h.y, h.x) - math.atan2(g_track[0].y, g_track[0].x)) < self.acceptance: # if the hits polar angle is withing range of the track
                g_track.append(h)
            else:
                grouped_by_phi.append(g_track)
                g_track = []
        if len(g_track) > 0:
            grouped_by_phi.append(g_track)

        tracks = []
        for g_track in grouped_by_phi:
            sorted_by_module = sort_by_module(g_track)
            current_track = []
            skipped = 0
            for h in sorted_by_module:
                if len(current_track) <= 0:
                    current_track.append(h)
                    continue
                if h.module_number == current_track[-1].module_number:
                    continue
                if h.module_number == current_track[-1].module_number + 2 or \
                        h.module_number == current_track[-1].module_number + 1:
                    current_track.append(h)
                elif skipped < 1 and (h.module_number == current_track[-1].module_number + 4 or
                                      h.module_number == current_track[-1].module_number + 3):
                    skipped += 1
                    current_track.append(h)
                else:
                    if len(current_track) >= self.min_track_length:
                        tracks.append(em.track(current_track))
                    skipped = 0
                    current_track = [h]

            if len(current_track) >= self.min_track_length:
                tracks.append(em.track(current_track))

        return tracks

    def solve3(self):
        grouped_by_phi = []
        g_track = []
        g_track_phi_sum = 0
        for h in self.hits:
            if len(g_track) <= 0:
                g_track.append(h)
                g_track_phi_sum += math.atan2(h.y, h.x)
            elif abs(math.atan2(h.y, h.x) - (g_track_phi_sum / len(g_track))) < self.acceptance: # if the hits polar angle is withing range of the track
                g_track.append(h)
                g_track_phi_sum += math.atan2(h.y, h.x)
            else:
                grouped_by_phi.append(g_track)
                g_track = []
                g_track_phi_sum = 0
        if len(g_track) > 0:
            grouped_by_phi.append(g_track)

        tracks = []
        for g_track in grouped_by_phi:
            sorted_by_module = sort_by_module(g_track)
            current_track = []
            skipped = 0
            for h in sorted_by_module:
                if len(current_track) <= 0:
                    current_track.append(h)
                    continue
                if h.module_number == current_track[-1].module_number:
                    continue
                if h.module_number == current_track[-1].module_number + 2 or \
                        h.module_number == current_track[-1].module_number + 1:
                    current_track.append(h)
                elif skipped < 1 and (h.module_number == current_track[-1].module_number + 4 or
                                      h.module_number == current_track[-1].module_number + 3):
                    skipped += 1
                    current_track.append(h)
                else:
                    if len(current_track) >= self.min_track_length:
                        tracks.append(em.track(current_track))
                    skipped = 0
                    current_track = [h]

            if len(current_track) >= self.min_track_length:
                tracks.append(em.track(current_track))

        return tracks

    def solve4(self):
        grouped_by_phi = []
        g_track = []
        g_track_phi_sum = 0
        for h in self.hits:
            if len(g_track) <= 0:
                g_track.append(h)
                g_track_phi_sum += math.atan2(h.y, h.x)
            elif check_acceptance(h, g_track, g_track_phi_sum, self.acceptance): # if the hits polar angle is withing range of the track
                g_track.append(h)
                g_track_phi_sum += math.atan2(h.y, h.x)
            else:
                grouped_by_phi.append(g_track)
                g_track = []
                g_track_phi_sum = 0
        if len(g_track) > 0:
            grouped_by_phi.append(g_track)

        tracks = []
        for g_track in grouped_by_phi:
            sorted_by_module = sort_by_module(g_track)
            current_track = []
            skipped = 0
            for h in sorted_by_module:
                if len(current_track) <= 0:
                    current_track.append(h)
                    continue
                if h.module_number == current_track[-1].module_number:
                    continue
                if h.module_number == current_track[-1].module_number + 2 or \
                        h.module_number == current_track[-1].module_number + 1:
                    current_track.append(h)
                elif skipped < 1 and (h.module_number == current_track[-1].module_number + 4 or
                                      h.module_number == current_track[-1].module_number + 3):
                    skipped += 1
                    current_track.append(h)
                else:
                    if len(current_track) >= self.min_track_length:
                        tracks.append(em.track(current_track))
                    skipped = 0
                    current_track = [h]

            if len(current_track) >= self.min_track_length:
                tracks.append(em.track(current_track))

        return tracks


def check_acceptance(hit, track, t_phi_sum, acceptance):
    deviation_sum = 0
    mean = t_phi_sum / len(track)
    for h in track:
        deviation_sum += (math.atan2(h.y, h.x) - mean)**2

    variance = deviation_sum / len(track)
    s_deviation = math.sqrt(variance)

    candiate_deviation = (math.atan2(hit.y, hit.x) - mean)**2
    return math.sqrt((candiate_deviation - s_deviation)**2) < 0.002
    # rss = 0
    # for h in track:
    #     rss += (math.atan2(hit.y, hit.x) - math.atan2(h.y, h.x))**2
    # mrss = rss / len(track)
    # return math.sqrt(mrss) < acceptance


def sort_by_phi(hits):
    phis = []
    for h in hits:
        phis.append(math.atan2(h.y, h.x))
    sorted_index = np.argsort(phis)
    x = np.array(hits)[sorted_index]
    return x


def sort_by_module(hits):
    modules = []
    for h in hits:
        modules.append(h.module_number)
    sorted_index = np.argsort(modules)
    x = np.array(hits)[sorted_index]
    return x
