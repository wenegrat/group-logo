"""Animated particle trails (nullschool-style) in the simulated surface flow around the logo.

Particles are advected (RK2) by the saved velocities, bilinear in space and linear in time between
saved frames; each is drawn with a fading tail whose brightness scales with the local speed, so calm
water stays clear. Particles that hit land or exceed their lifetime are re-seeded at random in water.
Writes transparent PNG frames (trails with the letters on top), one per saved output time.

Usage: python scripts/particle_trails.py <data dir> <out dir> [--n 8000] [--sub 6] [--trail 36]
           [--vmax 0.6] [--lw 0.7] [--alpha 0.8] [--color white] [--life 60,120] [--seed 0]
           [--stills d1,d2] [--taper 0..1] [--sketch scale,length,randomness] [--vmin m/s]
           [--colorby zeta|b] [--cmap color1,color2] [--crange value]
           [--wobble px] [--fade steps] [--soft] [--edgefade fx,fy] [--fpo frames] [--stopseed day] [--taperseed day] [--rampup days]
"""
import argparse, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from PIL import Image

p = argparse.ArgumentParser()
p.add_argument("data"); p.add_argument("out")
p.add_argument("--n", type=int, default=8000, help="number of particles")
p.add_argument("--sub", type=int, default=6, help="integration substeps per saved frame")
p.add_argument("--trail", type=int, default=36, help="tail length in substeps")
p.add_argument("--vmax", type=float, default=0.6, help="speed (m/s) at which tails are fully bright")
p.add_argument("--lw", type=float, default=0.7); p.add_argument("--alpha", type=float, default=0.8)
p.add_argument("--color", default="white")
p.add_argument("--life", default="60,120", help="lifetime range in substeps")
p.add_argument("--seed", type=int, default=0)
p.add_argument("--stills", default=None, help="comma-separated days: only save these frames")
p.add_argument("--colorby", default=None, choices=[None, "zeta", "b"], help="color strokes by vorticity sign or buoyancy anomaly")
p.add_argument("--cmap", default="#9FE8E0,#FFC2A8", help="two colors: (cyclonic,anticyclonic) or (cold,warm)")
p.add_argument("--crange", type=float, default=None, help="field value at full color (default: 0.5 for zeta, 98th pct for b)")
p.add_argument("--vmin", type=float, default=0.0, help="don't draw particles whose recent mean speed (m/s) is below this")
p.add_argument("--taper", type=float, default=0.0, help="0: constant width; 1: tails thin to zero (width ∝ position along tail)")
p.add_argument("--sketch", default=None, help="matplotlib path.sketch scale,length,randomness (redrawn randomly each frame)")
p.add_argument("--wobble", type=float, default=0.0, help="fixed hand-drawn wobble amplitude (pixels), attached to each particle")
p.add_argument("--wobblefreq", type=float, default=1.0, help="wobble frequency along the track (lower = gentler bends)")
p.add_argument("--fade", type=int, default=0, help="steps to fade particles in at birth and out at death (incl. hitting land)")
p.add_argument("--edgefade", default="0,0", help="fade tracks near the domain edges: width as fraction of domain (x,y)")
p.add_argument("--stopseed", type=float, default=None, help="day after which particles are no longer re-seeded; all tracks fade out by the last frame (for looping)")
p.add_argument("--noletters", action="store_true", help="don't draw the letters into the frames (overlay them separately, e.g. in HTML)")
p.add_argument("--rampup", type=float, default=0.0, help="days over which particles are introduced at the start (random birth times), mirroring the end")
p.add_argument("--endfade", action="store_true", help="with --stopseed: also fade all tracks out between stopseed and the end")
p.add_argument("--taperseed", type=float, default=None, help="day from which re-seeding tapers to zero at the last frame; every particle is given a death time before the end, so tracks thin out gradually and the last frame is empty (for looping)")
p.add_argument("--fpo", type=int, default=1, help="frames saved per output interval (must divide --sub); >1 for smoother motion")
p.add_argument("--soft", action="store_true", help="fade particles gradually around --vmin instead of a hard cutoff")
a = p.parse_args()

meta = open(os.path.join(a.data, "meta.txt")).read().split("\n")
Nx, Ny, nf, dx, X0, X1, Y0, Y1, W, H = meta[0].split()
Nx, Ny, nf, W, H = int(Nx), int(Ny), int(nf), int(W), int(H)
dx, X0, X1, Y0, Y1 = map(float, (dx, X0, X1, Y0, Y1))
times = np.array([float(t) for t in meta[1].split()])
U = np.fromfile(os.path.join(a.data, "u.f32"), dtype=np.float32).reshape(nf, Ny, Nx)
V = np.fromfile(os.path.join(a.data, "v.f32"), dtype=np.float32).reshape(nf, Ny, Nx)
land = np.fromfile(os.path.join(a.data, "land.u8"), dtype=np.uint8).reshape(Ny, Nx).astype(bool)
letters = Image.open(os.path.join(a.data, "letters.png")).convert("RGBA")
C = None
if a.colorby:
    C = np.fromfile(os.path.join(a.data, ("zeta" if a.colorby == "zeta" else "b") + ".f32"), dtype=np.float32).reshape(nf, Ny, Nx)
    crange = a.crange or (0.5 if a.colorby == "zeta" else float(np.percentile(C[-1][~np.fromfile(os.path.join(a.data, "land.u8"), dtype=np.uint8).reshape(Ny, Nx).astype(bool)], 98)))
c_lo, c_hi = [np.array(matplotlib.colors.to_rgb(c)) for c in a.cmap.split(",")]
pale = np.array([1.0, 1.0, 1.0])
Lx, Ly = X1 - X0, Y1 - Y0
dt = (times[1] - times[0]) * 86400 / a.sub                      # seconds per substep
rng = np.random.default_rng(a.seed)
water = np.argwhere(~land)

def seed(k):
    c = water[rng.integers(len(water), size=k)]
    return np.column_stack([X0 + (c[:, 1] + rng.random(k)) * dx, Y0 + (c[:, 0] + rng.random(k)) * dx])

def interp(F, x, y):
    """Bilinear interpolation of a cell-centered periodic field at positions (km)."""
    fx = (x - X0) / dx - 0.5; fy = (y - Y0) / dx - 0.5
    i0 = np.floor(fx).astype(int); j0 = np.floor(fy).astype(int); ax = fx - i0; ay = fy - j0
    i0 %= Nx; j0 %= Ny; i1 = (i0 + 1) % Nx; j1 = (j0 + 1) % Ny
    return ((1-ax)*(1-ay)*F[j0, i0] + ax*(1-ay)*F[j0, i1] + (1-ax)*ay*F[j1, i0] + ax*ay*F[j1, i1])

def vel(x, y, n, s):
    w = s / a.sub; n1 = min(n + 1, nf - 1)
    u = (1-w) * interp(U[n], x, y) + w * interp(U[n1], x, y)
    v = (1-w) * interp(V[n], x, y) + w * interp(V[n1], x, y)
    return u, v

def on_land(x, y):
    i = (np.floor((x - X0) / dx).astype(int)) % Nx; j = (np.floor((y - Y0) / dx).astype(int)) % Ny
    return land[j, i]

lo, hi = (int(v) for v in a.life.split(","))
P = seed(a.n); age = rng.integers(0, hi, size=a.n); life = rng.integers(lo, hi, size=a.n)
dying = np.zeros(a.n, bool)                                     # hit land: frozen in place while fading out
inactive = np.zeros(a.n, bool)                                  # retired after --stopseed (never re-seeded)
cur_t = times[0]
tapering = False
steps_per_day = a.sub / (times[1] - times[0])
# Ramp-up: each particle is born at a random step during the first --rampup days (invisible until then)
birth = rng.integers(0, max(int(a.rampup * steps_per_day), 1), size=a.n) if a.rampup > 0 else np.zeros(a.n, int)
pending = birth > 0
step_count = 0
phase = rng.random((a.n, 2)) * 2 * np.pi
km_per_px = (X1 - X0) / W
def wobble_offset(age_now):                                     # smooth, fixed function of particle age
    k = a.wobblefreq
    return a.wobble * km_per_px * (np.sin(0.35 * k * age_now + phase[:, 0]) + 0.5 * np.sin(0.9 * k * age_now + phase[:, 1])) / 1.5
hist = np.repeat(P[None], a.trail + 1, axis=0)                  # [trail+1, n, 2], newest last
spd = np.zeros((a.trail + 1, a.n))
fld = np.zeros((a.trail + 1, a.n))                               # coloring field along each tail
off = np.zeros((a.trail + 1, a.n))                               # wobble offset recorded with each tail point
stills = None if a.stills is None else {int(np.argmin(abs(times - float(d)))) for d in a.stills.split(",")}
os.makedirs(a.out, exist_ok=True)
fade = np.linspace(0, 1, a.trail + 1)[1:]                       # per segment, old → new
dpi = 100

def draw(n, jf):
    # Tail points displaced perpendicular to the tail by the fixed per-point wobble offsets
    pts = hist.copy()
    if a.wobble > 0:
        d = np.gradient(hist, axis=0)
        d[..., 0] = np.where(np.abs(d[..., 0]) > Lx / 2, 0, d[..., 0]); d[..., 1] = np.where(np.abs(d[..., 1]) > Ly / 2, 0, d[..., 1])
        norm = np.hypot(d[..., 0], d[..., 1]) + 1e-12
        pts[..., 0] += -d[..., 1] / norm * off; pts[..., 1] += d[..., 0] / norm * off
    # Segments between consecutive tail points, skipping periodic wrap-arounds
    a0, a1 = pts[:-1], pts[1:]
    ok = (np.abs(a1[..., 0] - a0[..., 0]) < Lx / 2) & (np.abs(a1[..., 1] - a0[..., 1]) < Ly / 2)
    bright = np.clip(spd[1:] / a.vmax, 0, 1) * fade[:, None]
    ok &= bright > 0.02
    # Per-particle opacity: fade in at birth, out at death; soft or hard speed cutoff
    pa = np.ones(a.n)
    if a.fade > 0:
        sm = lambda x: np.clip(x, 0, 1) ** 2 * (3 - 2 * np.clip(x, 0, 1))
        pa *= sm(age / a.fade) * sm((life - age) / a.fade)
    if a.vmin > 0:
        mean_spd = spd.mean(axis=0)
        if a.soft:
            x = np.clip((mean_spd - 0.5 * a.vmin) / a.vmin, 0, 1); pa *= x * x * (3 - 2 * x)
        else:
            pa *= (mean_spd >= a.vmin)
    if a.taperseed is not None:
        pa = pa * ~inactive
    if a.stopseed is not None:   # retired particles vanish (no re-seeding after --stopseed)
        pa = pa * ~inactive
        if a.endfade:            # optionally force everything to fade out between stopseed and the end
            x = np.clip((times[-1] - cur_t) / max(times[-1] - a.stopseed, 1e-9), 0, 1)
            pa = pa * (x * x * (3 - 2 * x))
    pa = pa * ~pending
    bright = bright * pa[None, :]
    ex, ey = (float(v) for v in a.edgefade.split(","))
    if ex > 0 or ey > 0:         # fade by distance of each segment's midpoint to the nearest domain edge
        mx = 0.5 * (a0[..., 0] + a1[..., 0]); my = 0.5 * (a0[..., 1] + a1[..., 1])
        sm = lambda x: np.clip(x, 0, 1) ** 2 * (3 - 2 * np.clip(x, 0, 1))
        fx = sm(np.minimum(mx - X0, X1 - mx) / (ex * Lx)) if ex > 0 else 1
        fy = sm(np.minimum(my - Y0, Y1 - my) / (ey * Ly)) if ey > 0 else 1
        bright = bright * fx * fy
    ok &= bright > 0.02
    segs = np.stack([a0[ok], a1[ok]], axis=1)
    rgba = np.zeros((ok.sum(), 4)); rgba[:, 3] = bright[ok]
    if C is None:
        rgba[:, :3] = matplotlib.colors.to_rgb(a.color)
    else:
        f = fld[1:][ok]
        if a.colorby == "zeta":     # sign → hue, magnitude → saturation (pale near zero)
            t = np.clip(np.abs(f) / crange, 0, 1)[:, None]
            base = np.where((f > 0)[:, None], c_lo, c_hi)
            rgba[:, :3] = pale * (1 - t) + base * t
        else:                       # cold → warm
            t = np.clip(f / crange, 0, 1)[:, None]
            rgba[:, :3] = c_lo * (1 - t) + c_hi * t
    fig = plt.figure(figsize=(W / dpi, H / dpi), dpi=dpi)
    ax = fig.add_axes([0, 0, 1, 1]); ax.set_axis_off(); ax.set_xlim(X0, X1); ax.set_ylim(Y0, Y1)
    widths = a.lw * (1 - a.taper + a.taper * np.broadcast_to(fade[:, None], ok.shape)[ok])
    if a.sketch:
        plt.rcParams["path.sketch"] = tuple(float(v) for v in a.sketch.split(","))
    ax.add_collection(LineCollection(segs, colors=rgba, linewidths=widths, capstyle="round"))
    tmp = os.path.join(a.out, "_lines.png")
    fig.savefig(tmp, dpi=dpi, transparent=True); plt.close(fig)
    lines = Image.open(tmp).convert("RGBA").resize((W, H))
    arr = np.array(lines).astype(np.float32); arr[..., 3] *= a.alpha
    frame = Image.fromarray(arr.clip(0, 255).astype(np.uint8))
    if not a.noletters:
        frame = Image.alpha_composite(frame, letters)
    name = f"frame_{n * a.fpo + jf:04d}.png" if stills is None else f"trails_day{times[n]:05.2f}.png"
    frame.save(os.path.join(a.out, name))

for n in range(nf):
    for s in range(a.sub):
        u1, v1 = vel(P[:, 0], P[:, 1], n, s)
        xm = P[:, 0] + 0.5 * dt * u1 / 1e3; ym = P[:, 1] + 0.5 * dt * v1 / 1e3
        u2, v2 = vel(xm, ym, n, s + 0.5)
        Pold = P.copy()
        P[:, 0] = X0 + (P[:, 0] + dt * u2 / 1e3 - X0) % Lx
        P[:, 1] = Y0 + (P[:, 1] + dt * v2 / 1e3 - Y0) % Ly
        P[dying] = Pold[dying]; u2 = np.where(dying, 0, u2); v2 = np.where(dying, 0, v2)
        age += 1
        if a.fade > 0:                     # hitting land: freeze and fade out instead of vanishing
            hit = on_land(P[:, 0], P[:, 1]) & ~dying
            P[hit] = Pold[hit]; dying |= hit; life[hit] = np.minimum(life[hit], age[hit] + a.fade)
            dead = age > life
        else:
            dead = on_land(P[:, 0], P[:, 1]) | (age > life)
        cur_t = times[n] + (s + 1) / a.sub * (times[1] - times[0])
        step_count += 1
        born = pending & (birth <= step_count)
        if born.any():           # bring in particles scheduled for this step, freshly seeded
            k = born.sum(); P[born] = seed(k); age[born] = 0; life[born] = rng.integers(lo, hi, size=k)
            hist[:, born] = P[born][None]; spd[:, born] = 0; fld[:, born] = 0; off[:, born] = 0
            pending &= ~born
        if a.stopseed is not None and cur_t >= a.stopseed:
            inactive |= dead; dead = np.zeros_like(dead)
        if a.taperseed is not None and cur_t >= a.taperseed:
            left = int((times[-1] - cur_t) * steps_per_day) - a.fade - a.sub   # steps left before fade-out must start (1-frame margin)
            if not tapering:     # give every live particle a random death time before the end
                tapering = True
                cap = age + rng.integers(0, max(left, 1), size=a.n) + a.fade
                life = np.minimum(life, cap)
            # re-seed only with probability decreasing linearly to zero at the end
            pkeep = np.clip((times[-1] - cur_t) / (times[-1] - a.taperseed), 0, 1)
            if left <= 0:        # too late for a new particle to fade in and out before the end
                pkeep = 0.0
            retire = dead & (rng.random(a.n) > pkeep)
            inactive |= retire; dead = dead & ~retire
        if dead.any():
            k = dead.sum(); P[dead] = seed(k); age[dead] = 0; life[dead] = rng.integers(lo, hi, size=k)
            if tapering:         # newly seeded particles must also be gone before the end
                left = int((times[-1] - cur_t) * steps_per_day) - a.fade - a.sub
                life[dead] = np.minimum(life[dead], rng.integers(0, max(left, 1), size=k) + a.fade)
            dying[dead] = False; phase[dead] = rng.random((k, 2)) * 2 * np.pi
            hist[:, dead] = P[dead][None]; spd[:, dead] = 0; fld[:, dead] = 0; off[:, dead] = 0
        hist = np.roll(hist, -1, axis=0); hist[-1] = P
        spd = np.roll(spd, -1, axis=0); spd[-1] = np.hypot(u2, v2)
        # wobble scaled down for near-stationary particles, whose perpendicular direction is ill-defined
        off = np.roll(off, -1, axis=0); off[-1] = wobble_offset(age) * np.clip(np.hypot(u2, v2) / 0.1, 0, 1) if a.wobble > 0 else 0
        if C is not None:
            w = s / a.sub; n1 = min(n + 1, nf - 1)
            fld = np.roll(fld, -1, axis=0)
            fld[-1] = (1 - w) * interp(C[n], P[:, 0], P[:, 1]) + w * interp(C[n1], P[:, 0], P[:, 1])
        step = a.sub // a.fpo
        if (s + 1) % step == 0:
            if stills is None:
                draw(n, (s + 1) // step - 1)
            elif n in stills and s + 1 == a.sub:
                draw(n, 0)
if os.path.exists(os.path.join(a.out, "_lines.png")):
    os.remove(os.path.join(a.out, "_lines.png"))
print(f"Done: {nf} frames advected, {'all' if stills is None else len(stills)} saved to {a.out}")
