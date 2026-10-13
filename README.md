# group-logo

An animated logo for the UMD Ocean Dynamics group, made with [Oceananigans.jl](https://github.com/CliMA/Oceananigans.jl).

The words "OCEAN DYNAMICS" are a set of full-depth islands in a rotating, stratified ocean. The
islands' coastlines heat the water next to them, and the warm coastal currents this sets up go
baroclinically unstable and shed eddies that gradually fill the domain. The web banner shows the
surface flow as animated particle trails around the lettering.

## Physics

**Model.** `HydrostaticFreeSurfaceModel`: hydrostatic Boussinesq equations on an f-plane with a
buoyancy tracer, implicit free surface (preconditioned conjugate-gradient solver on the immersed
grid), WENO advection for momentum (vector-invariant) and buoyancy, and no explicit viscosity or
diffusivity. The domain is doubly periodic with a flat bottom and free-slip walls everywhere.

**Initial state.** At rest, with uniform stratification, b = N² z (plus small random noise).

**Forcing.** A constant buoyancy flux F_b (warming) through every side wall between ocean and land,
at all depths. Because the heating is uniform with depth, the coastal buoyancy gradient it builds
is depth-independent. Once the free surface adjusts, it is in thermal-wind balance with a uniform
vertical shear and zero depth-averaged flow. Locally this is an Eady problem whose horizontal
buoyancy gradient grows until the rim currents go unstable. Nothing removes the heat.

**Parameters** (`scripts/heated_logo.jl`):

| | |
|---|---|
| Coriolis parameter | f = 10⁻⁴ s⁻¹ |
| Depth | H = 1000 m, 8 levels |
| Deformation radius | L_d = NH/f = 10 km, set equal to half the stroke width of the letters |
| Stratification | N² = (L_d f / H)² = 10⁻⁶ s⁻² |
| Wall buoyancy flux | F_b = Δb·L_d / t = 5×10⁻⁴ m s⁻² × 10 km / 2 days ≈ 2.9×10⁻⁵ m² s⁻³ |
| Grid | 5 cells per L_d (Δx = 2 km), 694 × 231 × 8 (1388 × 462 km, about 3:1) |
| Open water around the text | 17.5 L_d east/west, 19.5 L_d north/south |
| Run length | 35 days, output every 3 hours |

The wall flux is an idealized forcing: expressed as a heat flux through a thin wall into a wide
ocean, it is far larger than anything in nature.

**What to expect.**
- At the surface, warm water carries anticyclonic flow; at the bottom, the same warm water carries
  cyclonic flow.
- Because cyclones can become much stronger than anticyclones (ζ can far exceed f, but not fall
  much below −f), the surface vorticity shows broad, weak anticyclonic patches. The bottom
  vorticity shows small, intense cyclone cores.

**Known limitations.**
- Water sealed inside letters (the counters of O, D, A) is heated with no way to escape and warms
  without bound. It is dynamically isolated; `ZERO_LAKES=1` hides it in snapshots.
- Long straight, grid-aligned fronts appeared late in earlier runs. They are probably numerical and
  haven't been investigated.
- Velocities are saved every 3 hours and interpolated linearly in time for the particle trails, so
  long trails around translating eddies show step-like displacements; the 24-hour trails used for
  the banner keep this mostly invisible.

## How the logo becomes land

`src/Logo.jl` renders text with Makie/FreeType into a bitmap and samples it onto the model grid,
scaled so that a median stroke of the letters is `2 × cells_per_Ld` cells wide. The current layout
(`:text`) is "OCEAN DYNAMICS" in Montserrat ExtraBold (included, SIL Open Font License). An earlier
layout (`:emblem`: the Terrapin emblem followed by the text in a bubble font, with a drop shadow
that is also land) is kept for reproducing earlier runs.

`src/logo_config.jl` defines the layout once, and both the simulation and the plotting scripts
use it, so the drawn logo always sits exactly on the model's land. The figures draw the logo from
a 4×-finer version of the same mask, so outlines are smooth, while the model sees the
grid-resolution mask.

## Requirements

- Julia 1.12 and the packages in `Project.toml` (Oceananigans 0.113, CairoMakie, JLD2,
  FreeTypeAbstraction). `Manifest.toml` pins the exact versions used:
  `julia --project=. -e 'using Pkg; Pkg.instantiate()'`.
- For the web banner: Python 3 with the packages in `requirements.txt` (numpy, matplotlib,
  Pillow), and `ffmpeg` (with libvpx-vp9; on macOS the HEVC encoder `hevc_videotoolbox`) and
  `cwebp` from libwebp.
- **Only for the earlier `:emblem` layout**, two inputs that are not included in this repository:
  - The **emblem image**, a PNG with a transparent background. Put it at `assets/emblem.png` or set
    `EMBLEM_PATH`. The Terrapin emblem is a registered trademark of the University of Maryland
    and is not redistributed here.
  - The **Bubble Vibes** font (`BubbleVibes-Regular.otf`, from 1001fonts.com; licensed for personal
    use only). Put it at `fonts/bubblevibes/BubbleVibes-Regular.otf` or set `LOGO_FONT_PATH`.

## Running

```bash
# Check the layout (land mask and rendered logo) without simulating
julia --project=. scripts/preview_logo.jl 5             # → output/logo_preview.png

# Optional: verify the sign and magnitude of the wall flux on a small test case
julia --project=. scripts/test_wall_flux.jl

# Run the simulation (about 1 hour on 8 CPU threads; add `gpu` to use a GPU)
julia -t 8 --project=. scripts/heated_logo.jl           # → output/heated_text_694x231x8_{surface,bottom}.jld2

# Four-panel animation and snapshot figures
julia --project=. scripts/animate.jl heated_text_694x231x8
julia --project=. scripts/snapshot.jl heated_text_694x231x8 [frame]
```

The simulation refuses to overwrite existing output unless you add `overwrite` to the command.
`scripts/heated_logo.jl emblem` reruns the earlier emblem layout.

## Web banner

The banner (see `web/banner.html`) has three layers, back to front: a soft dark box whose darkness
follows the animation, a transparent video of white particle trails in the surface flow, and the
lettering as a separate sharp image. The particles are introduced gradually over the first 10
days, no longer replaced after day 20, and gone by day 28, so the video loops seamlessly (about
22.5 s at 10 frames per second).

```bash
# 1. Surface velocities (and land mask) for days 0–28
julia --project=. scripts/export_velocity.jl heated_text_694x231x8 surface 0 28 1100

# 2. Particle-trail frames (transparent PNGs, no lettering)
python scripts/particle_trails.py output/velocity_data_surface output/web_frames/trails \
    --n 5000 --trail 48 --life 96,384 --taper 1 --fade 24 --soft --vmin 0.05 --lw 0.9 --alpha 0.6 \
    --color white --edgefade 0.08,0.2 --rampup 10 --stopseed 20 --noletters

# 3. Videos with transparency: VP9 for most browsers (≈1.8 MB), HEVC for Safari (≈2.0 MB)
ffmpeg -framerate 10 -i output/web_frames/trails/frame_%04d.png -c:v libvpx-vp9 -pix_fmt yuva420p \
    -b:v 0 -crf 45 -row-mt 1 -deadline good -cpu-used 1 -g 240 output/banner/trails.webm
ffmpeg -framerate 10 -i output/web_frames/trails/frame_%04d.png -vf format=bgra \
    -c:v hevc_videotoolbox -alpha_quality 0.5 -q:v 40 -tag:v hvc1 -movflags +faststart output/banner/trails_hevc.mov

# 4. Lettering image (2× resolution) and the darkness curve for the background box
julia --project=. scripts/export_letters.jl heated_text_694x231x8 2200 output/banner/letters.png
cwebp -lossless output/banner/letters.png -o output/banner/letters.webp
python scripts/box_curve.py output/web_frames/trails output/banner/box_curve.json
```

`particle_trails.py --help` lists the other options (colouring by vorticity or buoyancy, hand-drawn
wobble, more frames per output interval, …).

**Plotting options** (environment variables for the Julia figures; colours are applied only when
plotting, so they can be changed without rerunning):
- `LOGO_SCHEME`: text layout `white` (default) or `yellow`; emblem layout `red` (default),
  `white`, `red_rim`, `two_tone`
- `ZERO_LAKES=1` (snapshots): draw trapped water as unheated

## Repository layout

```
scripts/
  heated_logo.jl      the simulation
  animate.jl          four-panel animation (surface/bottom vorticity, surface buoyancy, speed)
  snapshot.jl         single-frame figure with zooms
  preview_logo.jl     land mask + rendered logo, no simulation
  export_velocity.jl  surface velocities, vorticity, buoyancy and land mask for the web banner
  particle_trails.py  particle-trail frames for the web banner
  export_letters.jl   lettering image aligned with the banner
  box_curve.py        darkness curve for the banner's background box
  test_wall_flux.jl   check of the immersed-boundary wall flux
src/
  logo_config.jl      the logo layouts shared by simulation and plots
  Logo.jl             text/emblem rasterization and layout tools
  letter_styles.jl    decorative letter rendering (flat, bubble shading, outlines, shadows)
  emblem_overlay.jl   draws the styled logo on top of model output
  plotting.jl         small plotting helpers
  TextMask.jl         simple text-to-mask tools (used by archived scripts)
web/banner.html       HTML/CSS/JS for the banner in the website header
fonts/                Montserrat (used), plus Cinzel, Playfair Display and DM Serif Display from
                      the lettering study (all SIL Open Font License); fonts/bubblevibes/ is where
                      the Bubble Vibes font goes for the emblem layout (not included)
assets/emblem.png     where the emblem image goes for the emblem layout (not included)
archive/              earlier experiments: shallow-water letters, flow past letters, Eady
                      channels, lettering and colour studies, contour/outline and banner-style
                      studies, earlier web-banner scripts. Some use free fonts that are no longer
                      included (Google Fonts: Fredoka, Baloo 2, Titan One, Lilita One, Chewy,
                      Luckiest Guy, Sniglet, Yellowtail, Kaushan Script, Damion); download them
                      into fonts/<name>/ to rerun those scripts.
```

## License

The code is released under the MIT License (see `LICENSE`). The included fonts are under the SIL
Open Font License (see each `fonts/<name>/OFL.txt`). The emblem artwork and the Bubble Vibes font
are not part of this repository and have their own terms.
Rendered animations that show the emblem should not be redistributed without permission from the
University of Maryland.
