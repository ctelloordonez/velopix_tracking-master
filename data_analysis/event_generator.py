import os
import json
import event_model.event_model as em

import sys
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
    [plt.plot([hit.z for hit in track.hits], [hit.y for hit in track.hits], color='b') for track in tracks]
    plt.show()


if __name__ == '__main__':
    data_set = "small_dataset"
    events_tracks = []

    # events = get_events_from_folder(data_set)
    # write_tracks(random.sample(random.choice(events).real_tracks, 20), "test.txt")
    [plot_tracks(read_tracks(f"event_test_cases/case{i+1}.txt")) for i in range(10)]


