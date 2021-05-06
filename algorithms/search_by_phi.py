import math
import numpy as np
from event_model import event_model as em


class SearchByPhi:
    def __init__(self, event):
        self.hits = sort_by_phi(event.hits) # sort hits by phi

    def solve(self):
        tracks = [] # list of tracks found
        current_track = [] # list of hits representing the track currently under consideration
        current_track.append(self.hits[0]) # initialize by considering the first hit
        self.hits = self.hits[1:] # Take the first hit out of consideration for forming a track with itself
        skipped = 0 # variable to keep track of the number of modules skipped
        for h in self.hits: # loop over the hits in the event
            if abs(math.atan2(h.y, h.x) - math.atan2(current_track[-1].y, current_track[-1].x)) < 0.05: # if the hits polar angle is withing range of the track
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
                if len(current_track) > 1: # If the current track has a minimal number of hits
                    tracks.append(em.track(current_track)) # add it to the list of tracks found
                current_track = [] # Reset current track
                current_track.append(h) # Add the hit as the initial hit of the track
        if len(current_track) > 1:
            tracks.append(em.track(current_track)) # append final track if longer than 1

        return tracks


def sort_by_phi(hits):
    phis = []
    for h in hits:
        phis.append(math.atan2(h.y, h.x))
    sorted_index = np.argsort(phis)
    x = np.array(hits)[sorted_index]
    return x
