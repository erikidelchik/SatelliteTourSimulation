"""
Generate a KML file that animates a circular satellite orbit in Google Earth.

Usage:
    run the file to get the KML file with default parameters
    or, run with parameters as such:
    python simulate_orbit.py --alt 550 --inc 53 --period 5580 --step 10 --hours 3


Dependencies:
    pip install numpy simplekml
"""

import argparse
from datetime import datetime, timedelta, timezone
import numpy as np
import simplekml
from math import radians, sin, degrees, cos, atan2, sqrt
import random
import sys
from PIL import Image, ImageDraw

# constants
R_EARTH = 6378.137  # (WGS-84 equatorial radius)
MU = 3.986004418e5  # Earth's GM

# create mask img
size = 1024  # square texture
radius = 0.35  # hole radius = 35% of width

img = Image.new("RGBA", (size, size), color=(0, 0, 0, 255))  # full black
draw = ImageDraw.Draw(img)
hole_r = int(size * radius / 2)
cx = cy = size // 2
draw.ellipse((cx - hole_r, cy - hole_r, cx + hole_r, cy + hole_r), fill=(0, 0, 0, 0))
img.save("mask.png")
print("Saved mask.png")


# Haversine distance (km)
def gc_distance(lat1, lon1, lat2, lon2):
    # wrap delta-lon so it is always inside –180…180°
    dlon = ((lon2 - lon1 + 540) % 360) - 180

    fi1, fi2, dg = map(radians, (lat1, lat2, dlon))
    a = sin((fi2 - fi1) / 2)**2 + cos(fi1) * cos(fi2) * sin(dg / 2)**2
    return 2 * R_EARTH * atan2(sqrt(a), sqrt(1 - a))


# helper: true bearing from point 1 → point 2 (deg, 0 = North)
def bearing(lat1, lon1, lat2, lon2):
    fi1, fi2 = np.radians([lat1, lat2])
    dg = np.radians(lon2 - lon1)
    x = np.sin(dg) * np.cos(fi2)
    y = np.cos(fi1) * np.sin(fi2) - np.sin(fi1) * np.cos(fi2) * np.cos(dg)
    brng = np.degrees(np.arctan2(x, y))
    return (brng + 360) % 360  # wrap into 0‥360


# coordinate generator
def eci_to_llh(ang_displacement):
    """
    Convert ang_displacement to geodetic latitude/longitude.
    latitude  = asin(sin(i) * sin θ)
    longitude = θ (for equatorial ascending node at 0°E) – Earth’s rotation is ignored
    """
    lat = np.arcsin(np.sin(inc_rad) * np.sin(ang_displacement))
    lon = (np.rad2deg(ang_displacement) + 540) % 360 - 180   # wrap to –180…180
    return np.rad2deg(lat), lon


def gc_destination(lat, lon, bearing_deg, dist_km):
    """Return lat,lon reached by starting at (lat,lon) and going
       'bearing_deg' ° for 'dist_km' km on a sphere (WGS-84 radius)."""
    R = 6371.0
    fi1, g1, ang_displacement = map(radians, (lat, lon, bearing_deg))
    d = dist_km / R
    fi2 = sin(fi1) * cos(d) + cos(fi1) * sin(d) * cos(ang_displacement)
    fi2 = atan2(fi2, sqrt(1 - fi2 ** 2))  # lat
    g2 = g1 + atan2(sin(ang_displacement) * sin(d) * cos(fi1),
                    cos(d) - sin(fi1) * sin(fi2))  # lon
    return degrees(fi2), (degrees(g2) + 540) % 360 - 180  # wrap lon


def generate_evenly_spaced_targets(track_latlon,
                                   n_points=20,
                                   lateral_km=0):
    """
    Pick n_points along the given track_latlon so that each target
    lies roughly 'total_track_length / n_points' km apart.

    1. Compute the cumulative distance along track_latlon (in km).
    2. Let spacing = (total_length / n_points).
    3. Select the first sample whose cumulative_dist >= k * spacing, for k = 1..n_points.
    4. (Optionally) offset each selected point laterally by up to lateral_km in a random direction.
       If lateral_km == 0, no offset is applied and the target lies exactly on the computed ground track.

    Returns a list of (lat, lon) for each of the n_points targets.
    """

    # list of cumulative distances along track_latlon
    cumdist = [0.0]  # cumdist[i] = distance from track_latlon[0] to track_latlon[i]
    for i in range(1, len(track_latlon)):
        latA, lonA = track_latlon[i - 1]
        latB, lonB = track_latlon[i]
        d = gc_distance(latA, lonA, latB, lonB)
        cumdist.append(cumdist[-1] + d)

    total_length = cumdist[-1]
    if total_length <= 0:
        raise ValueError("Track is too short (zero total length).")

    # compute how far apart each target should be
    spacing_km = total_length / n_points

    out = []
    next_threshold = spacing_km  # we’ll look for the first index where cumdist >= this
    idx = 0

    # for each of the n_points, find the appropriate index
    for _ in range(n_points):
        # Advance idx until cumdist[idx] >= next_threshold (or reach end)
        while idx < len(cumdist) and cumdist[idx] < next_threshold:
            idx += 1

        # If we ran off the end, clamp to the last sample
        if idx >= len(cumdist):
            idx = len(cumdist) - 1

        base_lat, base_lon = track_latlon[idx]

        if lateral_km > 0:
            # If lateral_km > 0, randomly offset that base point by up to lateral_km
            bearing_rand = random.uniform(0, 360)
            tgt_lat, tgt_lon = gc_destination(base_lat, base_lon, bearing_rand, lateral_km)
        else:
            # No offset: target is exactly on the ground track
            tgt_lat, tgt_lon = base_lat, base_lon

        out.append((tgt_lat, tgt_lon))

        # Increase the threshold for the next target
        next_threshold += spacing_km

    return out



# argument parsing
parser = argparse.ArgumentParser(description="Create animated satellite track KML.")
parser.add_argument("--alt", type=float, default=250, help="Altitude above mean sea level in km")
parser.add_argument("--inc", type=float, default=0, help="Orbital inclination in degrees")
parser.add_argument("--period", type=float, default=None, help="Orbital period in seconds "
                                                               "(leave blank for Keplerian)")
parser.add_argument("--step", type=float, default=10, help="Time step in seconds")
parser.add_argument("--hours", type=float, default=2, help="Total simulated hours")
if 'ipykernel' in sys.modules:
    args = parser.parse_args(args=['--alt', '500', '--inc', '53', '--hours', '1.7'])
else:
    args = parser.parse_args()

h = args.alt
inc_rad = np.deg2rad(args.inc)
r = R_EARTH + h  # orbital radius (km)

if args.period is None:
    args.period = 2 * np.pi * np.sqrt(r ** 3 / MU)  # seconds

ang_velocity = 2 * np.pi / args.period  # angular velocity of the satellite (radians per sec)

# time vector
start_utc = datetime.now(timezone.utc)
t_span = np.arange(0, args.hours * 3600 + args.step, args.step)

# build KML document
kml = simplekml.Kml()

mask = kml.newscreenoverlay(name="mask")
mask.icon.href = "mask.png"  # local file → Google Earth loads it
mask.color = '00000000'            # begin with transparent mask

# centre the overlay on the screen
mask.overlayxy = simplekml.OverlayXY(
    x=0.5, y=0.5,
    xunits=simplekml.Units.fraction, yunits=simplekml.Units.fraction)

mask.screenxy = simplekml.ScreenXY(
    x=0.5, y=0.5,
    xunits=simplekml.Units.fraction, yunits=simplekml.Units.fraction)

# start at full width/height (1 = whole screen, fraction units)
mask.size = simplekml.Size(
    x=1.0, y=1.0,
    xunits=simplekml.Units.fraction, yunits=simplekml.Units.fraction)
MASK_ID = mask.id  # remember for animated updates

track = kml.newgxtrack(name=f"Satellite @ {h} km, inc {args.inc}°")


# Builds the tour
# list that contains the time when_list[i] at which the satellite has been in coord_list[i]
when_list = []

# list that contains coordinates of where the satellite has been
coord_list = []

# for every time stemp in time t we create a point (this is what makes the orbit)
for t in t_span:
    ang_displacement = ang_velocity * t
    # geographic coordinates used to pinpoint locations on earth
    lat, lon = eci_to_llh(ang_displacement)
    alt_m = h * 1000
    when_list.append((start_utc + timedelta(seconds=float(t))).isoformat())
    coord_list.append((lon, lat, alt_m))

track.newwhen(when_list)
track.newgxcoord(coord_list)

track.altitudemode = simplekml.AltitudeMode.absolute
track.extrude = 0
track.style.iconstyle.icon.href = "http://maps.google.com/mapfiles/kml/shapes/satellite.png"

tour = kml.newgxtour(name="Orbit fly-along")
playlist = tour.newgxplaylist()

track_samples = [eci_to_llh(ang_velocity * t) for t in t_span]

one_rev = int(args.period // args.step) + 1        # samples in 0 … 1 period
orbit_samples = track_samples[:one_rev]

# generate list of the targets
targets = generate_evenly_spaced_targets(orbit_samples,
                                         n_points=20,
                                         lateral_km=0)

targets_named = [(lat, lon, f"Target #{i + 1}") for i, (lat, lon) in enumerate(targets)]


# Give each target a stable id so AnimatedUpdate can address it
target_ids = {}
for lat, lon, name in targets_named:
    pm = kml.newpoint(name=name, coords=[(lon, lat)])
    pm.style.iconstyle.icon.href = (
        "http://maps.google.com/mapfiles/kml/paddle/red-circle.png")
    pm.style.iconstyle.scale = 1.0
    target_ids[name] = pm.id

step_cam = 10  # seconds between camera key‐frames
idx_step = int(step_cam / args.step)  # every idx_step we create new flyto
range_km = 700  # trigger radius, in km

watch_target = None  # hold the target (lat, lon, name) when in range_km, otherwise will remain None


def alpha_hex(a):  # 0‥255 → 'AA' two-digit hex
    return f"{max(0, min(255, int(a))):02x}"


for i in range(idx_step, len(t_span) - idx_step, idx_step):
    t_now = t_span[i]
    theta_now = ang_velocity * t_now  # theta = Current lap angle in radians
    sat_lat, sat_lon = eci_to_llh(theta_now)  # The geographical position of the satellite on the Earth's surface

    # find a target within range.
    # scan through all targets each time. if multiple are inside radius,
    # pick the first one we see. Otherwise, watch_target = None.
    watch_target = None
    for tgt_lat, tgt_lon, tgt_name in targets_named:
        if t_now > 180 and gc_distance(sat_lat, sat_lon, tgt_lat, tgt_lon) <= range_km:
            watch_target = (tgt_lat, tgt_lon, tgt_name)
            break

    # Build the next FlyTo
    fly = playlist.newgxflyto(
        gxduration=step_cam)  # key-frame are the time points where the camera state is explicitly defined
    # satellite coordinates
    fly.camera.latitude = sat_lat
    fly.camera.longitude = sat_lon
    # camera height in meters above sea level
    fly.camera.altitude = h * 1000
    # the altitude is measured above the sea level
    fly.camera.altitudemode = simplekml.AltitudeMode.absolute
    # creates an effect of continuous movement
    fly.gxflytomode = simplekml.GxFlyToMode.smooth

    # nearest target distance
    dist_km, *_ = min(
        ((gc_distance(sat_lat, sat_lon, tl, tn),)
         for tl, tn, _ in targets_named),
        key=lambda x: x[0])

    # alpha ramps 0 → 255 inside `range_km`
    if dist_km >= range_km:
        alpha = 0
    else:
        alpha = 255 * (1 - dist_km / range_km)  # linear fade

    # colour string for Google Earth: AABBGGRR, black = 000000

    colour = f"{alpha_hex(alpha)}000000"

    upd = simplekml.Update()
    if t_now > 180:
        upd.change = (f'<ScreenOverlay targetId="{MASK_ID}">'
                      f'  <color>{colour}</color>'
                      f'</ScreenOverlay>')

    # let GE interpolate opacity over the whole key-frame
    playlist.newgxanimatedupdate(gxduration=step_cam, update=upd)

    if watch_target is not None:
        tgt_lat, tgt_lon, _ = watch_target

        # heading point directly from satellite to target (calculated with bearing())
        fly.camera.heading = bearing(sat_lat, sat_lon, tgt_lat, tgt_lon)

        # tilt: look “down” at target.
        # ground distance:
        ground_km = gc_distance(sat_lat, sat_lon, tgt_lat, tgt_lon)
        tilt_rad = atan2(ground_km, h)  # relative inclination between length and width
        fly.camera.tilt = degrees(tilt_rad)

    else:
        # no target in range, normal chase‐cam orientation
        next_lat, next_lon = eci_to_llh(ang_velocity * t_span[i + idx_step])
        fly.camera.heading = bearing(sat_lat, sat_lon, next_lat, next_lon)  # fly and look at the next point
        fly.camera.tilt = 35  # a normal state where the camera looks down

outfile = f"orbit_tour_evenly_spaced_20_targets{h:.0f}km_{args.inc:.0f}deg.kml"
kml.save(outfile)
print(f"Tour-enabled KML saved to {outfile}")
