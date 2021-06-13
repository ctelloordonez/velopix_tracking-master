#!/usr/bin/python3

import event_model.event_model as em
import numpy as np
import validator.validator_lite as vl
import json

# Solvers
from algorithms.graph_dfs import graph_dfs
from algorithms.track_following import track_following
from algorithms.search_by_phi import SearchByPhi
from visual.print_phi import print_event_2d_phi, print_event_3d_phi, \
    print_event_2d_2phi, print_event_2d_phi_r, print_event_2d_2phi_r, print_event_3d_3phi, plot_projection, \
    plot_phi_by_module
from data_analysis.events_analysis import tracks_from_data

solutions = {}

# Get an event
f = open("../events/minibias/velo_event_0.json")
json_data = json.loads(f.read())
event = em.event(json_data)
f.close()

for s_number in range(0, event.number_of_modules):
    event.modules[s_number].z = event.modules[s_number].hits()[0].z

tracks = tracks_from_data(json_data, only_reconstructible=True)
ghost = vl.ghost_from_reconstruction(json_data, SearchByPhi(event).solve3())

search_by_phi = SearchByPhi(event)
solutions["search_by_phi"] = search_by_phi.solve()

save_to_file = True

print_event_2d_phi_r(event, tracks, filename='non_ghost_tracks', save_to_file=True)

# print_event_3d_phi(event, tracks, filename='cilindertry', save_to_file=False)
# plot_phi_by_module(event, solutions["search_by_phi"], filename='phi_by_module_search_by_phi', save_to_file=True)

# print_event_2d_phi(event, tracks, x=2, phix=0, phiy=1, filename="minibias event0 phi by z", save_to_file=True)
# print_event_2d_phi(event, tracks, x=0, phix=1, phiy=2, filename="minibias event0 theta by x", save_to_file=True)
# print_event_2d_phi(event, tracks, x=1, phix=2, phiy=0, filename="minibias event0 psi by y", save_to_file=True)
# print_event_2d_phi(event, tracks, x=1, phix=0, phiy=1, filename="minibias event0 phi by y", save_to_file=False)

# print_event_2d_phi_r(event, tracks, phix=0, phiy=1, filename="minibias event0 phi by r", save_to_file=True)
# print_event_2d_phi_r(event, tracks, phix=1, phiy=2, filename="minibias event0 theta by r", save_to_file=True)
# print_event_2d_phi_r(event, tracks, phix=2, phiy=0, filename="minibias event0 psi by r", save_to_file=True)

# print_event_2d_2phi(event, tracks, phix=[0,1], phiy=[1,2], filename="phi-theta", save_to_file=save_to_file)
# print_event_2d_2phi(event, tracks, phix=[1,2], phiy=[2,0], filename="theta-psi", save_to_file=save_to_file)
# print_event_2d_2phi(event, tracks, phix=[2,0], phiy=[0,1], filename="psi-phi", save_to_file=save_to_file)
# print_event_2d_2phi(event, tracks, phix=[0,1], phiy=[2,0], filename="phi-psi", save_to_file=save_to_file)

# print_event_2d_2phi_r(event, tracks, phix=0, phiy=1, x=2, filename="phi by angle(r, z)", save_to_file=save_to_file)
# print_event_2d_2phi_r(event, tracks, phix=1, phiy=2, x=0, filename="theta by angle(r, x)", save_to_file=save_to_file)
# print_event_2d_2phi_r(event, tracks, phix=2, phiy=0, x=1, filename="psi by angle(r, y)", save_to_file=save_to_file)

# print_event_3d_phi(event, tracks, x=0, y=1, phix=0, phiy=1, filename='minibias event0 phi by x-y', save_to_file=True)
# print_event_3d_phi(event, tracks, x=1, y=2, phix=1, phiy=2, filename='minibias event0 theta by y-z', save_to_file=True)
# print_event_3d_phi(event, tracks, x=2, y=0, phix=2, phiy=0, filename='minibias event0 psi by z-x', save_to_file=True)

# for o in range(-40, 101, 20):
#     print_event_3d_3phi(event, tracks, offset=o, filename=f'3d 3phi {o}', save_to_file=False)

# plot_projection(event, tracks, filename='projectionphi_phi', save_to_file=False)
