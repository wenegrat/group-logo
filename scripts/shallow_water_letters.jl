# Option B: "UMD Ocean Dynamics" as thickness anomalies in a 1.5-layer (reduced-gravity)
# shallow-water model on an f-plane. Each letter is a geostrophically balanced lens whose
# edge jet goes unstable and sheds eddies.
#
# Nondimensionalized with f = g′ = H = 1, so lengths are in deformation radii Ld = √(g′H)/f
# and times in 1/f (one inertial period = 2π).
#
# Usage: julia --project=. scripts/shallow_water_letters.jl [gpu]

include(joinpath(@__DIR__, "..", "src", "TextMask.jl"))
using .TextMask
using Oceananigans
using Oceananigans.Models: ShallowWaterModel
using Oceananigans.Units: Time
using Printf, Random

Random.seed!(1234)

arch = "gpu" in ARGS ? GPU() : CPU()

# ── Parameters ───────────────────────────────────────────────────────────────
text        = "UMD Ocean Dynamics"
Lx, Ly      = 256, 64          # domain size (Ld)
text_width  = 200              # text spans this width (Ld); stroke ≈ 2.4 Ld
cells_per_Ld = arch isa GPU ? 16 : 4
Nx, Ny      = Lx * cells_per_Ld, Ly * cells_per_Ld

f, g, H     = 1, 1, 1
Δh          = 0.3              # letter thickness anomaly (fraction of H); >0 → anticyclonic lenses
σ_smooth    = 0.5              # Gaussian smoothing of letter edges (Ld)
noise       = 1e-3             # random perturbation amplitude on h

Δt          = 0.2 / cells_per_Ld   # gravity-wave CFL ≈ 0.25 at √(gH) ≈ 1.3
stop_time   = 300              # ≈ 48 inertial periods
save_every  = 1

name = @sprintf("sw_letters_%dx%d", Nx, Ny)

# ── Grid, model ──────────────────────────────────────────────────────────────
grid = RectilinearGrid(arch; size = (Nx, Ny), x = (-Lx/2, Lx/2), y = (-Ly/2, Ly/2),
                       topology = (Periodic, Periodic, Flat), halo = (5, 5))

model = ShallowWaterModel(grid; gravitational_acceleration = g, coriolis = FPlane(; f),
                          momentum_advection = WENO(order = 5), tracers = :c)

# ── Initial condition: balanced letters ──────────────────────────────────────
xc, yc = xnodes(grid, Center()), ynodes(grid, Center())
M = text_mask(text, Array(xc), Array(yc); width = text_width)
M = gaussian_smooth(M, σ_smooth * cells_per_Ld)

hᵢ = H .+ Δh .* M .+ noise .* randn(Nx, Ny)
set!(model, h = hᵢ, c = M)

# Geostrophic transport: uh = -(g/f) h ∂y h,  vh = (g/f) h ∂x h
h = model.solution.h
uhᵢ = Field(@at((Face, Center, Center), h * ∂y(h)))
vhᵢ = Field(@at((Center, Face, Center), h * ∂x(h)))
compute!(uhᵢ); compute!(vhᵢ)
interior(uhᵢ) .*= -g / f
interior(vhᵢ) .*=  g / f
set!(model, uh = uhᵢ, vh = vhᵢ)

# ── Diagnostics and output ───────────────────────────────────────────────────
uh, vh, h = model.solution
u, v = uh / h, vh / h
ζ = Field(∂x(v) - ∂y(u))          # relative vorticity (units of f)

simulation = Simulation(model; Δt, stop_time)

mkpath("output")
simulation.output_writers[:fields] = JLD2Writer(model, (; h, ζ, c = model.tracers.c);
                                                filename = joinpath("output", name * ".jld2"),
                                                schedule = TimeInterval(save_every),
                                                overwrite_files = true)

wall = Ref(time_ns())
function progress(sim)
    elapsed = (time_ns() - wall[]) / 1e9
    @printf("t = %6.1f  iter = %6d  max|u| = %.3f  h ∈ [%.3f, %.3f]  wall = %.1f s\n",
            time(sim), iteration(sim), maximum(abs, u), minimum(h), maximum(h), elapsed)
    flush(stdout)
    wall[] = time_ns()
end
add_callback!(simulation, progress, TimeInterval(10))

Oceananigans.Diagnostics.erroring_NaNChecker!(simulation)
run!(simulation)
