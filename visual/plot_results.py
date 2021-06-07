#!/usr/bin/python3

import event_model.event_model as em
import validator.validator_lite as vl
import json

# Solvers
from algorithms.clustering import Clustering
from algorithms.graph_dfs import graph_dfs
from algorithms.track_following import track_following
from data_analysis.events_analysis import tracks_from_data
from algorithms.search_by_phi import SearchByPhi

# Visuals
from visual.base import print_event_evaluation_3d

solutions = {}

# Get an event
f = open("../events/small_dataset/velo_event_0.json")
json_data = json.loads(f.read())
event = em.event(json_data)
f.close()

for s_number in range(0, event.number_of_modules):
    event.modules[s_number].z = event.modules[s_number].hits()[0].z

# Solve with the classic method
search_by_phi = SearchByPhi(event)
solutions["search_by_phi"] = search_by_phi.solve3()

tracks = tracks_from_data(json_data, only_reconstructible=True)
print_event_evaluation_3d(event, real_tracks=tracks, found_tracks=solutions["search_by_phi"],
                          save_to_file=True, filename='search_by_phi_visual_evaluation')

# # Solve with the DFS method
# dfs = graph_dfs(
#   allow_cross_track=False,
#   allowed_skip_modules=1,
#   max_slopes=(0.7, 0.7),
#   max_tolerance=(0.3, 0.3)
# )
# solutions["dfs"] = dfs.solve(event)

# from visual.base import print_event_2d
# print_event_2d(event, solutions["track_following"], filename="classic_solution_xz.png", save_to_file=True)
# print_event_2d(event, solutions["track_following"], y=1, filename="classic_solution_yz.png", save_to_file=True)
#
# print_event_2d(event, solutions["dfs"], filename="dfs_solution_xz.png", track_color=4, save_to_file=True)
# print_event_2d(event, solutions["dfs"], y=1, filename="dfs_solution_yz.png", track_color=4, save_to_file=True)
# solutions["clustering"] = tracks_from_data(json_data=json_data)
#
# from visual.base import print_clustered_2d
# bins = Clustering().get_angle_bins(event, x=2, y=0)
# print_clustered_2d(event, [], bins, filename="clustering_solution_xz", with_modules=False, save_to_file=True)
#
# bins = Clustering().get_angle_bins(event, x=1, y=2)
# print_clustered_2d(event, [], bins, y=1, filename="clustering_solution_yz", with_modules=False, save_to_file=True)
#
# bins = Clustering().get_angle_bins(event, x=0, y=1)
# print_clustered_2d(event, [], bins, x=0, y=1, filename="clustering_solution_xy", with_modules=False, save_to_file=True)
#
# bins = Clustering().get_space_bins(event, coordinate=2)
# print_clustered_2d(event, [], bins, x=2, y=1, filename="clustering_space_yz", with_modules=False, save_to_file=True)
#
# bins = Clustering().get_space_bins(event, coordinate=0)
# print_clustered_2d(event, [], bins, x=2, y=0, filename="clustering_space_xz", with_modules=False, save_to_file=True)
#
# bins = Clustering().get_space_bins(event, coordinate=1)
# print_clustered_2d(event, [], bins, x=2, y=1, filename="clustering_space_y_yz", with_modules=False, save_to_file=True)
