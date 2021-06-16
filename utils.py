import math


def get_polar_coordinates(x, y):
    r = math.sqrt(x ** 2 + y ** 2)
    phi = math.atan2(x, y)
    return r, phi


def get_angle_bins(event, k=4, x=0, y=1):
    bins = [[] for i in range(k)]
    for hit in event.hits:
        hit_r, hit_phi = get_polar_coordinates(hit[x], hit[y])
        if hit_phi < 0:
            hit_phi = math.pi - hit_phi
        for b, cluster_bin in enumerate(bins):
            if b * (2 * math.pi / k) < hit_phi <= (b+1) * (2 * math.pi / k):
                cluster_bin.append(hit)
    return bins


def get_space_bins(event, k=4, coordinate=2):
    bins = [[] for i in range(k)]
    for hit in event.hits:
        if coordinate == 2:
            offset = 400
        else:
            offset = 40
        value = hit[coordinate] + offset
        if coordinate == 2:
            space_range = 1200
        else:
            space_range = 80
        for b, cluster_bin in enumerate(bins):
            if b * (space_range / k) < value <= (b+1) * (space_range / k):
                cluster_bin.append(hit)
    return bins
