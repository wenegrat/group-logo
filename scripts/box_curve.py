"""Darkness curve for the soft dark box behind the web banner, synced to the trails animation.

Reads the rendered trail frames, measures the amount of trail in each (sum of alpha), smooths it
over 5 frames and shapes it: the box follows the square root of the trail amount relative to its
peak, both while the trails build up and while they fade out, so it starts darkening as soon as
the first streaks appear and lingers slightly as they disappear. Writes a JSON list (one value per
frame, 0..1) to paste into the page's script (see web/banner.html).

Usage: python scripts/box_curve.py <frames dir> <out.json> [--stop-day 20] [--frames-per-day 8]
"""
import argparse, glob, json, os
import numpy as np
from PIL import Image

p = argparse.ArgumentParser()
p.add_argument("frames"); p.add_argument("out")
p.add_argument("--stop-day", type=float, default=20.0, help="day after which particles are no longer re-seeded")
p.add_argument("--frames-per-day", type=int, default=8)
a = p.parse_args()

files = sorted(glob.glob(os.path.join(a.frames, "frame_*.png")))
amt = np.array([np.array(Image.open(f).convert("RGBA"))[..., 3].astype(float).sum() for f in files])
amt = np.convolve(np.pad(amt / amt.max(), 2, mode="edge"), np.ones(5) / 5, mode="valid")
ipk = int(np.argmax(amt[: int(a.stop_day * a.frames_per_day)]))       # peak before re-seeding stops
curve = np.sqrt(np.clip(amt / amt[ipk], 0, 1))
json.dump([round(float(v), 3) for v in curve], open(a.out, "w"))
print(f"{len(curve)} frames, peak at day {ipk / a.frames_per_day:.2f}, saved {a.out}")
