#!/usr/bin/python3

import event_model.event_model as em
import validator.validator_lite as vl
import json

# Solvers
from algorithms.graph_dfs import graph_dfs
from algorithms.track_following import track_following
solutions = {}

# Get an event
f = open("../events/small_dataset/velo_event_0.json")
json_data = json.loads(f.read())
event = em.event(json_data)
f.close()

for s_number in range(0, event.number_of_modules):
  event.modules[s_number].z = event.modules[s_number].hits()[0].z

# Solve with the classic method
track_following = track_following()
solutions["classic"] = track_following.solve(event)

# Solve with the DFS method
dfs = graph_dfs(
  allow_cross_track=False,
  allowed_skip_modules=1,
  max_slopes=(0.7, 0.7),
  max_tolerance=(0.3, 0.3)
)
solutions["dfs"] = dfs.solve(event)

from visual.print_phi import print_event_2d_phi
print_event_2d_phi(event, solutions["classic"], filename="event_phi_classic")
print_event_2d_phi(event, solutions["dfs"], filename="event_phi_dfs")
