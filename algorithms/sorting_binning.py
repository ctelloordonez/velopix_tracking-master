
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

#method that takes as input the parameter on which it has to sort hits
#the possible parameters are 'x','y','z','pol_r','pol_phi','theta','psi'
#this happens in ascending order

def sort_by(parameter):
    dataset = 'bsphiphi'
    events = get_events_data_from_folder(dataset_folder[dataset], num_events=10, shuffle=False)
    events_sorted = []
    for event in events:
        if parameter == 'pol_r' or parameter == 'pol_phi':
            for _ in event.hits:
                _.update_polar()
        if parameter == 'psi':
            for _ in event.hits:
                _.update_psi()
        if parameter == 'theta':
            for _ in event.hits:
                _.update_theta()
        events_sorted.append(sorted(event.hits, key=attrgetter(parameter)))

    return events_sorted

#method that takes as input the required amount of bins,
#and the parameter on which it sorts and bins the hits
#this method starts by calling the sorting method above

def get_bins(bins_count,parameter):
    events= sort_by(parameter)
    all_events_bins=[]
    for event in events:
        number_of_hits_of_event=len(event)
        bins= [[] for i in range(bins_count)]
        for index in range(len(event)):
            for b,bin in enumerate(bins):
                if b*number_of_hits_of_event/bins_count<index<(b + 1)*number_of_hits_of_event/bins_count:
                    bin.append(event[index])
        all_events_bins.append(bins)
    return all_events_bins

#this method bins the hits by module.
#this method starts by calling the sorting method above

def get_bins_by_module():
    events= sort_by('z')
    all_events_bins=[]
    for event in events:
        print(event)
        number_of_hits_of_event=len(event)
        bins= [[] for i in range(52)]
        for hit in event:
            module = hit.module_number
            print(module)
            bins[module].append(hit)
        all_events_bins.append(bins)
        print(all_events_bins)
    return all_events_bins


get_bins_by_module()


'''
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
        events_sorted.append(sorted(event.hits, key=attrgetter('pol_r','z')))
    return events_sorted


def sort_by_pol_phi():
    dataset = 'bsphiphi'
    events = get_events_data_from_folder(dataset_folder[dataset], num_events=10, shuffle=False)
    events_sorted=[]
    for event in events:
        for _ in event.hits:
            _.update_polar()
        events_sorted.append(sorted(event.hits, key=attrgetter('pol_phi','z')))
    return events_sorted
'''
