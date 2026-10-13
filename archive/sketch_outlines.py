"""Hand-drawn ("xkcd"-style) vorticity outlines around the logo, using matplotlib's path.sketch
(the line-jitter mechanism behind plt.xkcd). Reads data from scripts/export_outline_data.jl and
writes transparent PNG frames: sketched outlines (at LINE_ALPHA opacity) with the letters on top.

Usage: python scripts/sketch_outlines.py <data dir> <out dir> [--levels 0.2] [--lw 1.2]
           [--sketch scale,length,randomness] [--alpha 0.5] [--days 12] [--seed 0]
           [--smooth cells] [--sign both|pos|neg] [--minlen km] [--colors cyclone,anticyclone]
           [--level-alphas a1,a2,...]
"""
import argparse, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from PIL import Image

p = argparse.ArgumentParser()
p.add_argument("data"); p.add_argument("out")
p.add_argument("--levels", default="0.2"); p.add_argument("--lw", type=float, default=1.2)
p.add_argument("--sketch", default="1,100,2"); p.add_argument("--alpha", type=float, default=0.5)
p.add_argument("--days", default=None, help="comma-separated days to render (default: all frames)")
p.add_argument("--seed", type=int, default=0)
p.add_argument("--smooth", type=float, default=0.0, help="Gaussian smoothing of the field, in grid cells (simplifies shapes)")
p.add_argument("--sign", default="both", choices=["both", "pos", "neg"], help="contour cyclones (pos), anticyclones (neg), or both")
p.add_argument("--minlen", type=float, default=0.0, help="drop contour loops shorter than this (km)")
p.add_argument("--colors", default="white,white", help="line colors for cyclones,anticyclones (names or hex)")
p.add_argument("--level-alphas", default=None, help="opacity per |level| (same order as --levels), multiplied by --alpha")
a = p.parse_args()

meta = open(os.path.join(a.data, "meta.txt")).read().split("\n")
nx, ny, nf, x0, x1, y0, y1, X0, X1, Y0, Y1, W, H = meta[0].split()
nx, ny, nf, W, H = int(nx), int(ny), int(nf), int(W), int(H)
x0, x1, y0, y1, X0, X1, Y0, Y1 = map(float, (x0, x1, y0, y1, X0, X1, Y0, Y1))
times = np.array([float(t) for t in meta[1].split()])
Z = np.fromfile(os.path.join(a.data, "zeta.f32"), dtype=np.float32).reshape(nf, ny, nx)   # Julia column-major
xs, ys = np.linspace(x0, x1, nx), np.linspace(y0, y1, ny)
letters = Image.open(os.path.join(a.data, "letters.png")).convert("RGBA")
signs = {"both": (-1, 1), "pos": (1,), "neg": (-1,)}[a.sign]
lev = sorted({s * float(v) for v in a.levels.split(",") for s in signs})

def smooth(f, sigma):
    """Gaussian smoothing on the doubly periodic grid (FFT)."""
    if sigma <= 0:
        return f
    ky = np.fft.fftfreq(f.shape[0])[:, None]; kx = np.fft.rfftfreq(f.shape[1])[None, :]
    G = np.exp(-2 * (np.pi * sigma) ** 2 * (kx ** 2 + ky ** 2))
    return np.fft.irfft2(np.fft.rfft2(f) * G, s=f.shape)
sketch = tuple(float(v) for v in a.sketch.split(","))

frames = range(nf) if a.days is None else [int(np.argmin(abs(times - float(d)))) for d in a.days.split(",")]
os.makedirs(a.out, exist_ok=True)
plt.rcParams["path.sketch"] = sketch
dpi = 100
for k, n in enumerate(frames, start=1):
    np.random.seed(a.seed + n)
    fig = plt.figure(figsize=(W / dpi, H / dpi), dpi=dpi)
    ax = fig.add_axes([0, 0, 1, 1]); ax.set_axis_off()
    ax.set_xlim(X0, X1); ax.set_ylim(Y0, Y1)
    F = smooth(Z[n].astype(np.float64), 2 * a.smooth)
    cpos, cneg = [c if not c.lstrip("#").isalnum() or len(c.lstrip("#")) != 6 else "#" + c.lstrip("#") for c in a.colors.split(",")]
    mags = [float(v) for v in a.levels.split(",")]
    lalpha = [1.0] * len(mags) if a.level_alphas is None else [float(v) for v in a.level_alphas.split(",")]
    def rgba(l):
        c = matplotlib.colors.to_rgb(cneg if l < 0 else cpos)
        return (*c, lalpha[min(range(len(mags)), key=lambda i: abs(mags[i] - abs(l)))])
    cs = ax.contour(xs, ys, F, levels=lev, colors=[rgba(l) for l in lev], linewidths=a.lw,
                    negative_linestyles="solid")
    if a.minlen > 0:          # redraw only contour pieces longer than minlen
        segs = [(l, seg) for l, level_segs in zip(cs.levels, cs.allsegs) for seg in level_segs]
        cs.remove()
        for l, seg in segs:
            if np.sum(np.hypot(*np.diff(seg, axis=0).T)) >= a.minlen:
                ax.plot(seg[:, 0], seg[:, 1], color=cneg if l < 0 else cpos, lw=a.lw, solid_capstyle="round")
    fig.savefig(os.path.join(a.out, "_lines.png"), dpi=dpi, transparent=True)
    plt.close(fig)
    lines = Image.open(os.path.join(a.out, "_lines.png")).convert("RGBA").resize((W, H))
    arr = np.array(lines).astype(np.float32); arr[..., 3] *= a.alpha
    lines = Image.fromarray(arr.clip(0, 255).astype(np.uint8))
    frame = Image.alpha_composite(lines, letters)
    name = f"frame_{k:03d}.png" if a.days is None else f"sketch_day{times[n]:05.2f}.png"
    frame.save(os.path.join(a.out, name))
os.remove(os.path.join(a.out, "_lines.png"))
print(f"Saved {len(frames)} frame(s) ({W} x {H}) to {a.out}")
