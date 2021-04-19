
import json
import os
from operator import attrgetter
import random
import math
import sys

# this gets the full path of the project root on your pc and adds it to the path
project_root = os.path.abspath(os.path.join(__file__,'..','..'))
if project_root not in sys.path:
    sys.path.append(project_root)

import event_model.event_model as em

# resource switch
dataset_folder={
    'small':'events/small_dataset',
    'bsphiphi': 'events/bsphiphi',
    'minibias': 'events/minibias'
}

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


def sort_by_x(e):
  return e[0]


def sort_by_y(e):
  return e[1]


def sort_by_z(e):
  return e[2]


def sort_by_angle(e):
  origin=(0,0,0)
  return math.atan2(origin[0],e[0])


# Iterate all events

""""
for (dirpath, dirnames, filenames) in os.walk("events"):
  print((dirpath, dirnames, filenames))
  for i, filename in enumerate(filenames):
    # Get an event
    f = open(os.path.realpath(os.path.join(dirpath, filename)))
    json_data = json.loads(f.read())
    event = em.event(json_data)
    print('hallo')
    f.close()
    print('halllo')
    print(event.hits.sort(key=sort_by_angle))

"""
dataset = 'bsphiphi'
events = get_events_data_from_folder(dataset_folder[dataset], num_events=10, shuffle=False)
for event in events:

  print(event.hits[0:10])
  print(sorted(event.hits, key=attrgetter('x','y')))

print("okd")
