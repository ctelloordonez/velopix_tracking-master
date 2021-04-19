
import json
import os
from operator import attrgetter
import random
import sys
import event_model.event_model as em

# this gets the full path of the project root on your pc and adds it to the path
project_root = os.path.abspath(os.path.join(__file__,'..','..'))
if project_root not in sys.path:
    sys.path.append(project_root)


# inherit classes and alter
def get_events_data_from_folder(data_set_folder, num_events = -1, shuffle = False):
    events = []
    # for (dirpath, dirnames, filenames) in os.walk(f"../events/{data_set_folder}"):
    for (dirpath, dirnames, filenames) in os.walk(os.path.join(project_root, data_set_folder)):
        if shuffle:
            random.shuffle(filenames)
        if num_events > 0 or num_events < len(filenames):
            filenames = filenames[:num_events]
        for i, filename in enumerate(filenames):
            # Get an event
            print(f'opening: {filename}')
            f = open(os.path.realpath(os.path.join(dirpath, filename)))
            json_data = json.loads(f.read())
            event = em.event(json_data,read_tracks=True)
            events.append(event)
            f.close()
            print(f'closing : {filename}')
    return events


def sort_by_x():
    dataset = 'bsphiphi'
    events = get_events_data_from_folder(dataset_folder[dataset], num_events=10, shuffle=False)
    events_sorted = []
    for event in events:
        events_sorted.append(sorted(event.hits, key=attrgetter('x', 'y')))
    return events_sorted


def sort_by_y():
    dataset = 'bsphiphi'
    events = get_events_data_from_folder(dataset_folder[dataset], num_events=10, shuffle=False)
    events_sorted = []
    for event in events:
        events_sorted.append(sorted(event.hits, key=attrgetter('y', 'x')))
    return events_sorted


def sort_by_z():
    dataset = 'bsphiphi'
    events = get_events_data_from_folder(dataset_folder[dataset], num_events=10, shuffle=False)
    events_sorted = []
    for event in events:
        events_sorted.append(sorted(event.hits, key=attrgetter('z', 'x')))
    return events_sorted


def sort_by_pol_r():
    dataset = 'bsphiphi'
    events = get_events_data_from_folder(dataset_folder[dataset], num_events=10, shuffle=False)
    events_sorted=[]
    for event in events:
        for _ in event.hits:
            _.update_polar()
        print(sorted(event.hits, key=attrgetter('pol_r','z')))
        events_sorted.append(sorted(event.hits, key=attrgetter('pol_r','z')))
    return events_sorted


def sort_by_pol_phi():
    dataset = 'bsphiphi'
    events = get_events_data_from_folder(dataset_folder[dataset], num_events=10, shuffle=False)
    events_sorted=[]
    for event in events:
        for _ in event.hits:
            _.update_polar()
        print(sorted(event.hits, key=attrgetter('pol_phi','z')))
        events_sorted.append(sorted(event.hits, key=attrgetter('pol_phi','z')))
    return events_sorted

