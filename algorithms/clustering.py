import math


class Clustering:
    def __init__(self, K=4):
        self.bins_count = K
        self.bins = [[] for i in range(self.bins_count)]

    def get_angle_bins(self, event, x=0, y=1):
        for hit in event.hits:
            hit_r, hit_phi = get_polar_coordinates(hit[x], hit[y])
            if hit_phi < 0:
                hit_phi = math.pi - hit_phi
            for b, cluster_bin in enumerate(self.bins):
                if b * (2 * math.pi / self.bins_count) < hit_phi <= (b+1) * (2 * math.pi / self.bins_count):
                    cluster_bin.append(hit)
        return self.bins


def get_polar_coordinates(x, y):
    r = math.sqrt(x ** 2 + y ** 2)
    phi = math.atan2(x, y)
    return r, phi

