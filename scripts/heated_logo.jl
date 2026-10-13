# The logo ("OCEAN DYNAMICS"; see src/logo_config.jl) as full-depth islands that heat the ocean
# through their side walls.
#
# The ocean starts at rest with uniform stratification N². A uniform buoyancy flux Fb (positive =
# warming) is imposed on every lateral face between ocean and land, at all depths. The coastal
# buoyancy gradient this builds is depth-independent, so after the free surface adjusts it is in
# thermal-wind balance with a uniform shear and zero depth-mean flow: locally an Eady problem whose
# M² grows with time until the rim currents go baroclinically unstable. Nothing removes the heat.
#
# Doubly periodic, hydrostatic, implicit free surface. Ld = N H / f = stroke / 2.
#
# Usage: julia -t 8 --project=. scripts/heated_logo.jl [emblem] [gpu] [overwrite]
#   default layout: text (Montserrat); `emblem`: Terrapin emblem + Bubble Vibes (earlier runs)

include(joinpath(@__DIR__, "..", "src", "logo_config.jl"))
using Oceananigans
using Oceananigans.Units
using Printf, Random

Random.seed!(1234)

arch = "gpu" in ARGS ? GPU() : CPU()

# ── Parameters ───────────────────────────────────────────────────────────────
layout       = "emblem" in ARGS ? :emblem : :text
cells_per_Ld = layout == :text ? 5 : 7.5   # text: Δx = 2 km; emblem (earlier run): ≈ 1.33 km
Nz           = 8

f   = 1e-4                          # s⁻¹
H   = 1000meters
Ld  = 10kilometers
N²  = (Ld * f / H)^2                # Ld = N H / f

# Wall buoyancy flux, chosen so the coastal buoyancy contrast reaches Δb★ over a width ~Ld in t★
Δb★, t★ = 5e-4, 2days
Fb      = Δb★ * Ld / t★             # m² s⁻³

noise      = 5e-6                   # random perturbation on b (m s⁻²)
stop_time  = layout == :text ? 35days : 25days
save_every = layout == :text ? 3hours : 6hours

Δx = Ld / cells_per_Ld

# Land mask (logo layout defined in src/logo_config.jl)
land = logo_layout(cells_per_Ld; layout).land
Nx, Ny = size(land)
Lx, Ly = Nx * Δx, Ny * Δx

name = "heated_$(layout)_$(Nx)x$(Ny)x$(Nz)"

# Refuse to clobber earlier output unless explicitly asked
outputs = [joinpath("output", name * suffix) for suffix in ("_surface.jld2", "_bottom.jld2")]
if any(isfile, outputs) && !("overwrite" in ARGS)
    error("Output files for $name already exist; rerun with `overwrite` to replace them: $(filter(isfile, outputs))")
end

@info @sprintf("Nx × Ny × Nz = %d × %d × %d, Lx = %.0f km, Ly = %.0f km, Δx = %.1f km, land fraction = %.3f",
               Nx, Ny, Nz, Lx / 1e3, Ly / 1e3, Δx / 1e3, sum(land) / length(land))
@info @sprintf("N² = %.2e s⁻², wall buoyancy flux Fb = %.2e m² s⁻³ (≈ %.0f W m⁻² for α = 2e-4 K⁻¹)",
               N², Fb, Fb * 1025 * 3990 / (9.81 * 2e-4))
flush(stderr)

# ── Grid with the logo as full-depth islands ──────────────────────────────────
underlying_grid = RectilinearGrid(arch; size = (Nx, Ny, Nz),
                                  x = (-Lx/2, Lx/2), y = (-Ly/2, Ly/2), z = (-H, 0),
                                  topology = (Periodic, Periodic, Bounded), halo = (7, 7, 4))

grid = ImmersedBoundaryGrid(underlying_grid, GridFittedBottom(ifelse.(land, 0.0, -H)))

# ── Model ────────────────────────────────────────────────────────────────────
# Positive flux on immersed lateral faces adds buoyancy to the ocean (checked in scripts/test_wall_flux.jl)
wall_flux = FluxBoundaryCondition(Fb)
b_bcs = FieldBoundaryConditions(immersed = ImmersedBoundaryCondition(west = wall_flux, east = wall_flux,
                                                                     south = wall_flux, north = wall_flux))

model = HydrostaticFreeSurfaceModel(grid;
                                    coriolis = FPlane(; f),
                                    buoyancy = BuoyancyTracer(),
                                    tracers = :b,
                                    boundary_conditions = (; b = b_bcs),
                                    free_surface = ImplicitFreeSurface(),
                                    momentum_advection = WENOVectorInvariant(),
                                    tracer_advection = WENO())

bᵢ(x, y, z) = N² * z + noise * randn()
set!(model, b = bᵢ)

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
    @printf("t = %s  iter = %6d  Δt = %s  max|u| = %.3f  max|v| = %.3f  max|w| = %.2e  max b′ = %.2e  wall = %.1f s\n",
            prettytime(sim), iteration(sim), prettytime(sim.Δt),
            maximum(abs, u), maximum(abs, v), maximum(abs, w), maximum(b) - N² * znodes(grid, Center())[Nz], elapsed)
    flush(stdout)
    wall[] = time_ns()
end
add_callback!(simulation, progress, TimeInterval(1day))

Oceananigans.Diagnostics.erroring_NaNChecker!(simulation)
run!(simulation)
