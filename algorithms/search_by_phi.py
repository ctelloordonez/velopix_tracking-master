import math
import numpy as np
from event_model import event_model as em


class SearchByPhi:
    def __init__(self, hits):
        self.hits = sort_by_phi(hits)

    def solve(self):
        tracks = []
        current_track = []
        current_track.append(self.hits[0])
        self.hits = self.hits[1:]
        skipped = 0
        for h in self.hits:
            if abs(math.atan2(h.y, h.x) - math.atan2(current_track[0].y, current_track[0].x)) < 0.01:
                if skipped == 1 and h.module_number == current_track[0].module_number + 2:
                    current_track.insert(0, h)
                elif skipped == 1 and h.module_number == current_track[-1].module_number - 2:
                    current_track.append(h)
                elif h.module_number == current_track[0].module_number - 2:
                    current_track.insert(0, h)
                elif h.module_number == current_track[-1].module_number + 2:
                    current_track.append(h)
                elif skipped <= 0 and h.module_number == current_track[0].module_number - 4:
                    current_track.insert(0, h)
                    skipped += 1
                elif skipped <= 0 and h.module_number == current_track[-1].module_number + 4:
                    current_track.append(h)
                    skipped += 1
                else:
                    skipped = 0
                    if len(current_track) > 1:
                        tracks.append(em.track(current_track))
                    current_track = []
                    current_track.append(h)
            else:
                skipped = 0
                if len(current_track) > 1:
                    tracks.append(em.track(current_track))
                current_track = []
                current_track.append(h)
        tracks.append(em.track(current_track))


                # if h.module_number == current_track[0].module_number - 2:
                #     if abs(math.atan2(h.y, h.x) - math.atan2(current_track[0].y, current_track[0].x)) < 0.1:
                #         current_track.insert(0, h)
                #
                # elif h.module_number == current_track[-1].module_number + 2:
                #     if abs(math.atan2(h.y, h.x) - math.atan2(current_track[-1].y, current_track[-1].x) < 0.1):
                #         current_track.append(h)
                #
                # elif skipped <= 0:
                #     print('hi')
                #     if h.module_number == current_track[0].module_number - 4:
                #         if abs(math.atan2(h.y, h.x) - math.atan2(current_track[0].y, current_track[0].x)) < 0.1:
                #             current_track.insert(0, h)
                #             skipped += 1
                #     elif h.module_number == current_track[-1].module_number + 4:
                #         if abs(math.atan2(h.y, h.x) - math.atan2(current_track[-1].y, current_track[-1].x) < 0.1):
                #             current_track.append(h)
                #             skipped += 1
                #     else:
                #         skipped = 0
                #         tracks.append(em.track(current_track))
                #         current_track = []
                #         current_track.append(h)
                # else:
                #     skipped = 0
                #     tracks.append(em.track(current_track))
                #     current_track = []
                #     current_track.append(h)

        return tracks


def sort_by_phi(hits):
    phis = []
    for h in hits:
        phis.append(math.atan2(h.y, h.x))
    sorted_index = np.argsort(phis)
    x = np.array(hits)[sorted_index]
    # for i in x:
    #     print(f'{math.atan2(i.y, i.x)}, {i.module_number}')
    return x
