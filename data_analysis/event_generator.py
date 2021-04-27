import os
import json
import event_model.event_model as em
import random

import sys
from collections import defaultdict
import matplotlib.pyplot as plt


project_root = os.path.abspath(os.path.join(__file__,'..','..'))
if project_root not in sys.path:
    sys.path.append(project_root)


def get_events_from_folder(data_set_folder):
    events = []
    # for (dirpath, dirnames, filenames) in os.walk(f"../events/{data_set_folder}"):
    for (dirpath, dirnames, filenames) in os.walk(os.path.join(project_root, f"events/{data_set_folder}")):
        for i, filename in enumerate(filenames):
            # Get an event
            # print(f'opening: {filename}')
            f = open(os.path.realpath(os.path.join(dirpath, filename)))
            json_data = json.loads(f.read())
            event = em.event(json_data, read_tracks=True)
            f.close()
            # print(f'closing : {filename}')

            events.append(event)
    return events


def track_to_string(track):
    return ", ".join([hit_to_string(hit) for hit in track]) + "\n"


def hit_to_string(hit):
    return str(hit.id) + " " + str(hit.module_number) + " " + \
           str(hit.x) + " " + str(hit.y) + " " + str(hit.z)


def write_tracks(tracks, file_name):
    f = open(file_name, 'w')
    lines = [track_to_string(track) for track in tracks]
    f.writelines(lines)
    f.close()


def read_tracks(file_name):
    f = open(file_name, 'r')
    tracks = []
    while True:
        line = f.readline().split("\n")[0]
        if not line:
            break
        hit_strings = [hit.split(" ") for hit in line.split(", ")]
        hits = [em.hit(float(hit[2]), float(hit[3]), float(hit[4]),
                       int(hit[0]), int(hit[1])) for hit in hit_strings]
        tracks.append(em.track(hits))
    f.close()
    return tracks


def plot_tracks(tracks):
    [plt.plot([hit.z for hit in track.hits],
              [hit.y for hit in track.hits], color='b') for track in tracks]
    plt.show()


def plot_modules(modules):
    [plt.scatter([hit.z for hit in module.hits()],
                 [hit.y for hit in module.hits()]) for module in modules]
    min_y = min([min(module.hits(), key=lambda hit: hit.y).y for module in modules])
    max_y = max([max(module.hits(), key=lambda hit: hit.y).y for module in modules])
    [plt.plot([sum([hit.z for hit in module.hits()])/len(module.hits())]*2,
              [min_y, max_y]) for module in modules]
    plt.show()


def plot_tracks_and_modules(tracks, modules):
    [plt.plot([hit.z for hit in track.hits],
              [hit.y for hit in track.hits], color='black') for track in tracks]
    plot_modules(modules)


def tracks_to_modules(tracks):
    modules = defaultdict(list)
    [[modules[hit.module_number].append(hit) for hit in track.hits] for track in tracks]
    return sorted([em.module(key, set([hit.z for hit in modules[key]]),
                   0, len(modules[key]), modules[key]) for key in modules.keys()],
                  key=lambda module: module.module_number)


def generate_test_tracks(allowed_modules: list = range(52), num_tracks=10, num_test_events=1, dataset="small_dataset"):
    events = random.sample(get_events_from_folder(dataset), num_test_events)
    total_tracks = []
    for event in events:
        tracks = [em.track([hit for hit in track.hits if hit.module_number in allowed_modules])
                  for track in event.real_tracks]
        tracks = [track for track in tracks if len(track.hits) > 0]
        total_tracks.append(random.sample(tracks, num_tracks))
    return total_tracks


if __name__ == '__main__':
    data_set = "small_dataset"
    events_tracks = []
    test_tracks = generate_test_tracks(allowed_modules=[2, 4, 6, 8])[0]
    modules = tracks_to_modules(test_tracks)
    plot_tracks_and_modules(test_tracks, modules)

    # events = get_events_from_folder(data_set)
    # write_tracks(random.sample(random.choice(events).real_tracks, 20), "test.txt")
    # [plot_tracks(read_tracks(f"event_test_cases/case{i+1}.txt")) for i in range(10)]
    # print(tracks_to_modules(read_tracks("event_test_cases/case10.txt")))


