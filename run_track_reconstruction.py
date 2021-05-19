#!/usr/bin/python3

import json
import os
from event_model import event_model as em
from validator import validator_lite as vl

# Solvers
from algorithms.track_following import track_following
from algorithms.search_by_phi import SearchByPhi
from algorithms.search_by_constant import SearchByConstant
from algorithms.forward_search import ForwardSearch
from algorithms.template_matching import TemplateMatching
solutions = {
  "search_by_phi": []
}
validation_data = []

# Instantiate algorithm
track_following = track_following()

# Iterate all events
for (dirpath, dirnames, filenames) in os.walk("C:/Users/tjerk/Documents/GitHub/velopix_tracking-master/events/small_dataset"):
    for i, filename in enumerate(filenames):
        # Get an event
        f = open(os.path.realpath(os.path.join(dirpath, filename)))
        json_data = json.loads(f.read())
        event = em.event(json_data)
        f.close()

        # Do track reconstruction

        print("Reconstructing event %i..." % (i))

        # Append the solution and json_data
      
        solutions["search_by_phi"].append(ForwardSearch(event).solve())
        validation_data.append(json_data)

# Validate the solutions
for k, v in iter(sorted(solutions.items())):
    print("\nValidating tracks from %s:" % (k))
    vl.validate_print(validation_data, v)
    print()
