import math
import json
import events
from event_model import event_model as em

# Iterate all events
for (dirpath, dirnames, filenames) in os.walk("events"):
  for i, filename in enumerate(filenames):
    # Get an event
    f = open(os.path.realpath(os.path.join(dirpath, filename)))
    json_data = json.loads(f.read())
    event = em.event(json_data)
    f.close()

print(event.hits[0])

def filter_by_x(e):
  return e[0]


def filter_by_y(e):
  return e[1]


def filter_by_z(e):
  return e[2]


def sort_by_angle(e):
  origin=(0,0,0)
  return math.atan2(origin[0],e[0])


hits = [(3, 2, 1), (4, 1, 3)]

hits.sort(key=sort_by_angle)

print(hits)


