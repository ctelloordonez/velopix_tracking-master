import json
import math
import os

import numpy as np
import scipy.stats as stats
import matplotlib.pyplot as plt

import event_model.event_model as em
from validator.validator_lite import MCParticle


def tracks_from_data(json_data):
  tracks = []
  hits = []
  for hid, (x, y, z) in enumerate(zip(json_data["x"], json_data["y"], json_data["z"])):
    hits.append(em.hit(x, y, z, hid))
  description = json_data["montecarlo"]["description"]
  particles = json_data["montecarlo"]["particles"]
  for p in particles:
    d = {description[i]: p[i] for i in range(len(description))}
    trackhits = [hits[hit_number] for hit_number in d["hits"]]
    mcp = MCParticle(d.get("key", 0), d.get("pid", 0), d.get("p", 0), d.get("pt", 0), d.get("eta", 0),
                     d.get("phi", 0), d.get("charge", 0), trackhits)
    tracks.append(em.track(trackhits))
  return tracks


def noise_from_data(json_data):
  noise = 0
  hits_in_tracks = 0
  tracks = []
  hits = []
  for hid, (x, y, z) in enumerate(zip(json_data["x"], json_data["y"], json_data["z"])):
    hits.append(em.hit(x, y, z, hid))
  description = json_data["montecarlo"]["description"]
  particles = json_data["montecarlo"]["particles"]
  for p in particles:
    d = {description[i]: p[i] for i in range(len(description))}
    trackhits = [hits[hit_number] for hit_number in d["hits"]]
    hits_in_tracks += len(trackhits)

  noise = len(hits) - hits_in_tracks
  return noise


def origins_from_data(json_data, limit=False):
  origins = []
  hits = []
  for hid, (x, y, z) in enumerate(zip(json_data["x"], json_data["y"], json_data["z"])):
    hits.append(em.hit(x, y, z, hid))
  description = json_data["montecarlo"]["description"]
  particles = json_data["montecarlo"]["particles"]
  for p in particles:
    d = {description[i]: p[i] for i in range(len(description))}
    trackhits = [hits[hit_number] for hit_number in d["hits"]]
    if limit:
      if -15 <= trackhits[0].x <= 15 and -10 <= trackhits[0].y <= 10:
        origins.append(trackhits[0])
    else:
      origins.append(trackhits[0])

  return origins


def distribution_of_tracks():
  events_tracks = []
  for (dirpath, dirnames, filenames) in os.walk("../events/small_dataset"):
    for i, filename in enumerate(filenames):
      # Get an event
      f = open(os.path.realpath(os.path.join(dirpath, filename)))
      json_data = json.loads(f.read())
      event = em.event(json_data)
      f.close()

      events_tracks.append(tracks_from_data(json_data=json_data))
  total_tracks = 0
  for event in events_tracks:
    total_tracks += len(event)

  mu = total_tracks / len(events_tracks)

  aux = 0
  for event in events_tracks:
    aux += ((len(event) - mu) ** 2)

  variance = aux / len(events_tracks)

  sigma = math.sqrt(variance)

  print(f'Mean: {mu}')
  print(f'Variance: {variance}')
  print(f'Standard deviation: {sigma}')

  x = np.linspace(mu - 4*sigma, mu + 4*sigma, 100)
  plt.plot(x, stats.norm.pdf(x, mu, sigma), label='pdf')
  plt.title('Distribution of #tracks over 10 events')
  plt.text(-300, .0020, fr'$\mu={mu},\ \sigma={round(sigma, 4)}$')
  # plt.plot(x, stats.norm.cdf(x, mu, sigma), label='cdf')
  plt.legend()
  plt.grid(True)
  plt.show()


def distribution_of_noise():
  events_noise = []
  for (dirpath, dirnames, filenames) in os.walk("../events/minibias"):
    for i, filename in enumerate(filenames):
      # Get an event
      print(f'opening: {filename}')
      f = open(os.path.realpath(os.path.join(dirpath, filename)))
      json_data = json.loads(f.read())
      event = em.event(json_data)
      f.close()

      events_noise.append(noise_from_data(json_data=json_data))
      print(f'{filename} closed')

  total_noise = 0
  for event in events_noise:
    total_noise += event

  mu = total_noise / len(events_noise)

  aux = 0
  for event in events_noise:
    aux += ((event - mu) ** 2)

  variance = aux / len(events_noise)

  sigma = math.sqrt(variance)

  print(f'Mean: {mu}')
  print(f'Variance: {variance}')
  print(f'Standard deviation: {sigma}')

  x = np.linspace(mu - 4 * sigma, mu + 4 * sigma, 100)
  plt.plot(x, stats.norm.pdf(x, mu, sigma), label='pdf')
  plt.title('Distribution of #noise over 10 events')
  plt.text(-250, .0025, fr'$\mu={mu},\ \sigma={round(sigma, 4)}$')
  # plt.plot(x, stats.norm.cdf(x, mu, sigma), label='cdf')
  plt.legend()
  plt.grid(True)
  plt.show()


def noise_histogram():
  events_noise = []
  for (dirpath, dirnames, filenames) in os.walk("../events/minibias"):
    for i, filename in enumerate(filenames):
      # Get an event
      print(f'opening: {filename}')
      f = open(os.path.realpath(os.path.join(dirpath, filename)))
      json_data = json.loads(f.read())
      event = em.event(json_data)
      f.close()

      events_noise.append(noise_from_data(json_data=json_data))
      print(f'{filename} closed')

  x = np.array(events_noise)
  plt.hist(x, bins=100)
  plt.show()


def tracks_histogram():
  events_tracks = []
  for (dirpath, dirnames, filenames) in os.walk("../events/bsphiphi"):
    for i, filename in enumerate(filenames):
      # Get an event
      print(f'opening: {filename}')
      f = open(os.path.realpath(os.path.join(dirpath, filename)))
      json_data = json.loads(f.read())
      event = em.event(json_data)
      f.close()
      print(f'{filename} closed')

      events_tracks.append(len(tracks_from_data(json_data=json_data)))

  x = np.array(events_tracks)
  plt.hist(x, bins=get_bins(x))
  plt.show()


def tracks_by_noise():
  events_noise = []
  events_tracks = []
  for (dirpath, dirnames, filenames) in os.walk("../events/bsphiphi"):
    for i, filename in enumerate(filenames):
      # Get an event
      print(f'opening: {filename}')
      f = open(os.path.realpath(os.path.join(dirpath, filename)))
      json_data = json.loads(f.read())
      event = em.event(json_data)
      f.close()
      print(f'{filename} closed')

      events_noise.append(noise_from_data(json_data=json_data))
      events_tracks.append(len(tracks_from_data(json_data=json_data)))

  x = np.array(events_tracks)
  y = np.array(events_noise)

  plt.scatter(x, y, s=1, marker='o')

  m, b = np.polyfit(x, y, 1)
  plt.plot(x, m*x + b, c='y', linewidth=0.5)

  plt.xlabel('#tracks')
  plt.ylabel('#noise')
  plt.title('#tracks by #noise on bsphiphi')
  plt.grid(True)
  plt.show()


def track_origin_analysis():
  events_origins = []
  events_origins_limited = []
  for (dirpath, dirnames, filenames) in os.walk("../events/small_dataset"):
    for i, filename in enumerate(filenames):
      # Get an event
      f = open(os.path.realpath(os.path.join(dirpath, filename)))
      json_data = json.loads(f.read())
      event = em.event(json_data)
      f.close()

      events_origins.append(origins_from_data(json_data=json_data))
      events_origins_limited.append(origins_from_data(json_data=json_data, limit=True))

  total_origins = 0
  events_origin_means = []
  percentajes = []
  for e in range(len(events_origins)):
    percentajes.append(len(events_origins_limited[e])/len(events_origins[e]))

    event_sum = np.zeros(2)
    for origin in events_origins[e]:
      event_sum = event_sum + np.array([origin[0], origin[1]])
    event_mean = event_sum / len(events_origins[e])
    events_origin_means.append(event_mean)

  for p in percentajes:
    print(p)


def get_bins(x):
    q25, q75 = np.percentile(x, [.25, .75])
    bin_width = 2 * (q75 - q25) * len(x) ** (-1 / 3)
    bins = round((x.max() - x.min()) / bin_width)
    print("Freedman–Diaconis number of bins:", bins)
    return bins


tracks_histogram()
