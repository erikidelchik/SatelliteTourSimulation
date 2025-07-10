Generate an **animated satellite ground‑track** with camera fly‑along and dynamic overlay for playback in Google Earth.

---

## 🌍 Overview

`simulate_orbit.py` builds a KML file that

* plots the full ground track of a fictitious circular orbit;
* drops **20 evenly‑spaced target placemarks** along the first revolution;
* creates a **mask overlay** that fades in/out as the spacecraft approaches any target: <img width="1568" height="961" alt="image" src="https://github.com/user-attachments/assets/3c065bce-979a-4cfb-bc0c-95ff6d06bc41" />

* scripts a gx\:Tour so Google Earth automatically chases the satellite, tilts toward a nearby target, and blends the overlay for a cinematic feel.

The result is a self‑contained `orbit_tour_… .kml` ready to double‑click and play in Google Earth (desktop).

---

## ✨ Features

| Feature            | Details                                                                                                          |
| ------------------ | ---------------------------------------------------------------------------------------------------------------- |
| Parametric orbit   | Altitude, inclination, period, time‑step and total duration are command‑line flags.                              |
| Keplerian default  | If `--period` is omitted the script solves Kepler’s third law for a circular orbit at the chosen altitude.       |
| gx\:Track          | Full 3‑D track with per‑sample timestamps so you can scrub the time slider.                                      |
| Auto‑camera        | Smooth fly‑to key‑frames keep the satellite centred; heading and tilt adjust toward the current in‑range target. |
| Dynamic mask       | A PNG screen overlay whose alpha ramps 0→255 as the satellite comes within `range_km` of any target.             |
| Anti‑meridian safe | Longitudes and distance calculations are wrapped so nothing breaks at ±180°.                                     |

---

## 🚀 Quick start

```bash
# 1. Install dependencies (Python ≥ 3.9)
pip install numpy simplekml pillow

# 2. Run with defaults (500 km, 53° inc, 1.7 h tour)
python simulate_orbit.py

# 3. Or customise e.g. a 550 km × 30° orbit for 3 h
python simulate_orbit.py --alt 550 --inc 30 --hours 3

# The script prints
#   Saved mask.png
#   Tour-enabled KML saved to orbit_tour_evenly_spaced_20_targets550km_30deg.kml

```

## ▶️ How to run

1. Open the generated KML file in Google Earth Pro on desktop.

2. In the left‑hand Places panel, expand Tours (if it’s collapsed) and click Orbit fly‑along:
<img width="348" height="268" alt="image" src="https://github.com/user-attachments/assets/4b1047d3-88c5-46e8-8058-425e533afaed" />

3. Press the Play Tour button (usually at the bottom‑right of the 3‑D view) to watch the animation:
<img width="349" height="68" alt="image" src="https://github.com/user-attachments/assets/4dd859db-473a-4e97-9367-4909be821f95" />


## 🗄️ Project structure

```
.
├── simulate_orbit.py   # main script
├── mask.png            # auto‑generated donut‑hole overlay
├── orbit_tour_*.kml    # output KML tour(s)

```

---

## ⚙️ Command‑line reference

| Flag       | Default | Description                                      |
| ---------- | ------- | ------------------------------------------------ |
| `--alt`    | `250`   | Orbit altitude above mean sea level (km).        |
| `--inc`    | `0`     | Inclination (degrees).                           |
| `--period` | *calc*  | Orbital period (s). Leave blank to auto‑compute. |
| `--step`   | `10`    | Simulation time‑step (s).                        |
| `--hours`  | `2`     | Total simulated time (h).                        |

Example for a fast LEO demo:

```bash
python simulate_orbit.py --alt 400 --inc 97.7 --hours 1.3 --step 5
```

