#!/usr/bin/python3

import event_model.event_model as em
import validator.validator_lite as vl
import json

# Solvers
from algorithms.graph_dfs import graph_dfs
from algorithms.track_following import track_following
solutions = {}

# Get an event
f = open("../events/minibias/velo_event_0.json")
json_data = json.loads(f.read())
event = em.event(json_data)
f.close()

for s_number in range(0, event.number_of_modules):
  event.modules[s_number].z = event.modules[s_number].hits()[0].z

# # Solve with the classic method
# track_following = track_following()
# solutions["classic"] = track_following.solve(event)
#
# # Solve with the DFS method
# dfs = graph_dfs(
#   allow_cross_track=False,
#   allowed_skip_modules=1,
#   max_slopes=(0.7, 0.7),
#   max_tolerance=(0.3, 0.3)
# )
# solutions["dfs"] = dfs.solve(event)

from visual.print_phi import print_event_2d_phi, print_event_3d_phi, print_event_2d_2phi
from data_analysis.events_analysis import tracks_from_data

tracks = tracks_from_data(json_data, only_reconstructible=True)
print_event_2d_phi(event, tracks, x=2, filename="minibias event0 phi by z", save_to_file=True)
# print_event_2d_2phi(event, tracks, filename="xz-xy")
# print_event_2d_phi(event, [], x=1, filename="empty_phi_y")
# print_event_2d_phi(event, [], filename="empty_phi_z")

# print_event_3d_phi(event, tracks, x=0, y=1, filename='empty3D_phi_xy')
# print_event_3d_phi(event, tracks, x=0, y=2, rotate=True, filename='empty3D_phi_xz')
# print_event_3d_phi(event, tracks, x=1, y=2, filename='empty3D_phi_yz')
# print_event_2d_phi(event, solutions["dfs"], filename="event_phi_dfs")
