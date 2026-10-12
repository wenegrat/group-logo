# The logo as full-depth islands in a baroclinically unstable Eady channel.
#
# Uniform stratification N² and uniform cross-channel buoyancy gradient M² = ∂b/∂y, with the
# initial flow in thermal-wind balance, f ∂u/∂z = -M², and zero depth-averaged velocity.
# Periodic in x, walls in y. The letters are land at every depth (GridFittedBottom with
# bottom height at the surface). Implicit free surface (PCG solver on the immersed grid).
#
# Scales are set by Ri = N² f² / M⁴ and by Ld = N H / f relative to the letter stroke width.
# The grid has `cells_per_Ld` cells per deformation radius.
#
# Logos:
#   vibes (default): "UMD" over "OCEAN DYNAMICS" in Bubble Vibes, Ld = stroke / 2
#   text:            "UMD Ocean Dynamics" in Heros Bold on one line, Ld = stroke
#
# Usage: julia -t 8 --project=. scripts/eady_islands.jl [vibes | text] [gpu] [overwrite]

include(joinpath(@__DIR__, "..", "src", "TextMask.jl"))
include(joinpath(@__DIR__, "..", "src", "Logo.jl"))
using .TextMask
using Oceananigans
using Oceananigans.Units
using Printf, Random

Random.seed!(1234)

arch = "gpu" in ARGS ? GPU() : CPU()

# ── Parameters ───────────────────────────────────────────────────────────────
logo         = "text" in ARGS ? :text : :vibes
cells_per_Ld = 10                   # 1 km grid spacing (was 5 cells per Ld, 2 km)
Nz           = 8

f   = 1e-4                          # s⁻¹
H   = 1000meters
Ld  = 10kilometers
Ri  = 4

N²  = (Ld * f / H)^2                # Ld = N H / f
M²  = sqrt(N²) * f / sqrt(Ri)       # Ri = N² f² / M⁴

Δx = Ld / cells_per_Ld

# Land mask on cell centers. The domain is 896 km × 308 km, leaving ≈ 9.4 Ld of open water east and
# west of the logo and ≈ 7.7 Ld between the logo and the north/south walls.
if logo == :vibes
    vibes = FTFont(fontfile("bubblevibes/BubbleVibes-Regular.otf"))
    stroke_per_Ld = 2                                   # Ld = stroke / 2
    domain_in_Ld = (89.6, 30.8)
    land = logo_mask(two_line(vibes); cells_per_stroke = stroke_per_Ld * cells_per_Ld,
                     grid = round.(Int, domain_in_Ld .* cells_per_Ld))
    Nx, Ny = size(land)
else
    text            = "UMD Ocean Dynamics"
    stroke_fraction = 0.0143                            # measured stroke width / text width for this font
    Nx, Ny          = 512, 128
    text_width      = Ld / stroke_fraction              # Ld = stroke
    xs = range(-(Nx - 1) / 2, (Nx - 1) / 2, length = Nx) .* Δx
    ys = range(-(Ny - 1) / 2, (Ny - 1) / 2, length = Ny) .* Δx
    land = text_mask(text, collect(xs), collect(ys); width = text_width) .> 0.5
end
Lx, Ly = Nx * Δx, Ny * Δx

noise      = 1e-2 * M² * Ld         # random perturbation on b
stop_time  = 10days
save_every = 6hours

name = "eady_$(logo)_$(Nx)x$(Ny)x$(Nz)"

# Refuse to clobber earlier output unless explicitly asked
outputs = [joinpath("output", name * suffix) for suffix in ("_surface.jld2", "_bottom.jld2")]
if any(isfile, outputs) && !("overwrite" in ARGS)
    error("Output files for $name already exist; rerun with `overwrite` to replace them: $(filter(isfile, outputs))")
end

@info @sprintf("%s: Nx × Ny × Nz = %d × %d × %d, Lx = %.0f km, Ly = %.0f km, Δx = %.1f km, land fraction = %.3f",
               logo, Nx, Ny, Nz, Lx / 1e3, Ly / 1e3, Δx / 1e3, sum(land) / length(land))
@info @sprintf("N² = %.2e s⁻², M² = %.2e s⁻², shear M²/f = %.2e s⁻¹, ΔU over depth = %.2f m/s",
               N², M², M² / f, M² / f * H)
@info @sprintf("Eady growth rate ≈ 0.31 f/√Ri = %.2e s⁻¹ (e-folding %.2f days)",
               0.31f / sqrt(Ri), sqrt(Ri) / (0.31f) / day)

# ── Grid with letters as full-depth islands ──────────────────────────────────
underlying_grid = RectilinearGrid(arch; size = (Nx, Ny, Nz),
                                  x = (-Lx/2, Lx/2), y = (-Ly/2, Ly/2), z = (-H, 0),
                                  topology = (Periodic, Bounded, Bounded), halo = (7, 7, 4))

bottom_height = ifelse.(land, 0.0, -H)

grid = ImmersedBoundaryGrid(underlying_grid, GridFittedBottom(bottom_height))

# ── Model ────────────────────────────────────────────────────────────────────
model = HydrostaticFreeSurfaceModel(grid;
                                    coriolis = FPlane(; f),
                                    buoyancy = BuoyancyTracer(),
                                    tracers = :b,
                                    free_surface = ImplicitFreeSurface(),
                                    momentum_advection = WENOVectorInvariant(),
                                    tracer_advection = WENO())

# Thermal wind with zero depth average: u = -(M²/f)(z + H/2)
bᵢ(x, y, z) = N² * z + M² * y + noise * randn()
uᵢ(x, y, z) = -(M² / f) * (z + H / 2)
set!(model, b = bᵢ, u = uᵢ)

# ── Simulation ───────────────────────────────────────────────────────────────
simulation = Simulation(model; Δt = 10minutes, stop_time)
conjure_time_step_wizard!(simulation, IterationInterval(20), cfl = 0.3, max_Δt = 30minutes)

u, v, w = model.velocities
b = model.tracers.b
ζ = Field(∂x(v) - ∂y(u))

mkpath("output")
simulation.output_writers[:surface] = JLD2Writer(model, (; ζ, b, u, v);
                                                 filename = joinpath("output", name * "_surface.jld2"),
                                                 indices = (:, :, Nz),
                                                 schedule = TimeInterval(save_every),
                                                 overwrite_files = true)

simulation.output_writers[:bottom] = JLD2Writer(model, (; ζ, b);
                                                filename = joinpath("output", name * "_bottom.jld2"),
                                                indices = (:, :, 1),
                                                schedule = TimeInterval(save_every),
                                                overwrite_files = true)

wall = Ref(time_ns())
function progress(sim)
    elapsed = (time_ns() - wall[]) / 1e9
    @printf("t = %s  iter = %6d  Δt = %s  max|u| = %.3f  max|v| = %.3f  max|w| = %.2e  wall = %.1f s\n",
            prettytime(sim), iteration(sim), prettytime(sim.Δt),
            maximum(abs, u), maximum(abs, v), maximum(abs, w), elapsed)
    flush(stdout)
    wall[] = time_ns()
end
add_callback!(simulation, progress, TimeInterval(1day))

Oceananigans.Diagnostics.erroring_NaNChecker!(simulation)
run!(simulation)
