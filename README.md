# group-logo

An animated logo for the UMD Ocean Dynamics group, made with [Oceananigans.jl](https://github.com/CliMA/Oceananigans.jl).

The logo (the Terrapin emblem followed by "OCEAN DYNAMICS" in bubble letters) is a set of
full-depth islands in a rotating, stratified ocean. The islands' coastlines heat the water next to
them, and the warm coastal currents this sets up go baroclinically unstable and shed eddies that
gradually fill the domain.

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
| Grid | 7.5 cells per L_d (Δx ≈ 1.33 km), 700 × 170 × 8 |
| Open water around the logo | 9.4 L_d east/west, 7.7 L_d north/south |
| Run length | 25 days, output every 6 hours |

The wall flux is an idealized forcing: expressed as a heat flux through a thin wall into a wide
ocean, it is far larger than anything in nature.

**What to expect.**
- At the surface, warm water carries anticyclonic flow; at the bottom, the same warm water carries
  cyclonic flow.
- Because cyclones can become much stronger than anticyclones (ζ can far exceed f, but not fall
  much below −f), the surface vorticity shows broad, weak anticyclonic patches. The bottom
  vorticity shows small, intense cyclone cores.

**Known limitations.**
- Water trapped inside the logo (sealed-off pockets in the emblem and letters) is heated with no
  way to escape and warms without bound. It is dynamically isolated; `ZERO_LAKES=1` hides it in
  snapshots.
- Long straight, grid-aligned fronts appear late in the runs. They are probably numerical and
  haven't been investigated.
- Fine details of the emblem (head, fist, feet) are only a few grid cells across as land. The
  artwork drawn on top is far more detailed than what the model sees.

## How the logo becomes land

`src/Logo.jl` renders text with Makie/FreeType into a bitmap, and loads the emblem PNG as a
silhouette (alpha ≥ 0.5, with enclosed gaps filled). It places them side by side, adds the drop
shadow under the text, and samples the result onto the model grid. The scale is set so that a
median stroke of the letters is `2 × cells_per_Ld` cells wide.

`src/logo_config.jl` defines the layout once, and both the simulation and the plotting scripts
use it, so the drawn logo always sits exactly on the model's land. The figures draw the logo from
a 4×-finer version of the same mask, so outlines are smooth, while the model sees the
grid-resolution mask.

## Requirements

- Julia 1.12 and the packages in `Project.toml` (Oceananigans 0.113, CairoMakie, JLD2,
  FreeTypeAbstraction). `Manifest.toml` pins the exact versions used:
  `julia --project=. -e 'using Pkg; Pkg.instantiate()'`.
- **Two inputs that are not included in this repository:**
  - The **emblem image**, a PNG with a transparent background. Put it at `assets/emblem.png` or set
    `EMBLEM_PATH`. The Terrapin emblem is a registered trademark of the University of Maryland
    and is not redistributed here.
  - The **Bubble Vibes** font (`BubbleVibes-Regular.otf`, from 1001fonts.com; licensed for personal
    use only). Put it at `fonts/bubblevibes/BubbleVibes-Regular.otf` or set `LOGO_FONT_PATH`.

## Running

```bash
# Check the layout (land mask and rendered logo) without simulating
julia --project=. scripts/preview_logo.jl 7.5          # → output/logo_preview.png

# Optional: verify the sign and magnitude of the wall flux on a small test case
julia --project=. scripts/test_wall_flux.jl

# Run the simulation (about 1 hour on 8 CPU threads; add `gpu` to use a GPU)
julia -t 8 --project=. scripts/heated_logo.jl          # → output/heated_emblem_700x170x8_{surface,bottom}.jld2

# Animation and snapshot figures
julia --project=. scripts/animate.jl heated_emblem_700x170x8
julia --project=. scripts/snapshot.jl heated_emblem_700x170x8 [frame]
```

The simulation refuses to overwrite existing output unless you add `overwrite` to the command.

**Plotting options** (environment variables; colours are applied only when plotting, so they can
be changed without rerunning):
- `LOGO_SCHEME`: `red` (default), `white`, `red_rim`, `two_tone`
- `ZERO_LAKES=1` (snapshots): draw trapped water as unheated

## Repository layout

```
scripts/
  heated_logo.jl      the simulation
  animate.jl          four-panel animation (surface/bottom vorticity, surface buoyancy, speed)
  snapshot.jl         single-frame figure with zooms
  preview_logo.jl     land mask + rendered logo, no simulation
  test_wall_flux.jl   check of the immersed-boundary wall flux
src/
  logo_config.jl      the logo layout shared by simulation and plots
  Logo.jl             text/emblem rasterization and layout tools
  letter_styles.jl    decorative letter rendering (bubble shading, outlines, shadows)
  emblem_overlay.jl   draws the styled logo on top of model output
  plotting.jl         small plotting helpers
  TextMask.jl         simple text-to-mask tools (used by archived scripts)
fonts/bubblevibes/    where the Bubble Vibes font goes (not included)
assets/emblem.png     where the emblem image goes (not included)
archive/              earlier experiments: shallow-water letters, flow past letters, Eady
                      channels, lettering and colour studies. Some of these use other free
                      fonts (Google Fonts: Fredoka, Baloo 2, Titan One, Lilita One, Chewy,
                      Luckiest Guy, Sniglet, Yellowtail, Kaushan Script, Damion) that are no
                      longer included; download them into fonts/<name>/ to rerun those scripts.
```

## License

The code is released under the MIT License (see `LICENSE`). This does not cover the emblem
artwork or the Bubble Vibes font, which are not part of this repository and have their own terms.
Rendered animations that show the emblem should not be redistributed without permission from the
University of Maryland.
