import copy
import math

from event_model.event_model import hit, track


class clustering:
    def __init__(self, K=4):
        self.bins_count = K
        self.bins = [[] for i in range(self.bins_count)]

    def get_bins(self, event):
        for hit in event.hits:
            hit_r, hit_phi = get_polar_coordinates(hit.x, hit.y)
            for b, bin in enumerate(self.bins):
                if b * (2 * math.pi / self.bins_count) < abs(hit_phi) <= (b+1) * (2 * math.pi / self.bins_count):
                    bin.append(hit)
        return self.bins


def get_polar_coordinates(x, y):
    r = math.sqrt(x ** 2 + y ** 2)
    phi = math.atan2(x, y)
    if phi < 0:
        phi = math.pi - phi
    return r, phi


def tracks_from_data(json_data):
    tracks = []
    hits = []
    for hid, (x, y, z) in enumerate(zip(json_data["x"], json_data["y"], json_data["z"])):
        hits.append(hit(x, y, z, hid))
    description = json_data["montecarlo"]["description"]
    particles = json_data["montecarlo"]["particles"]
    for p in particles:
        d = {description[i]: p[i] for i in range(len(description))}
        trackhits = [hits[hit_number] for hit_number in d["hits"]]
        tracks.append(track(trackhits))
    return tracks
