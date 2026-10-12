# Option 1: "UMD Ocean Dynamics" as a solid immersed boundary in a 2D (x-y) nonhydrostatic
# flow. A sponge band restores the flow to uniform upward velocity U, so each letter sheds
# its own wake upward and the text stays readable underneath. The domain is doubly periodic;
# wakes leaving the top re-enter at the bottom and are erased by the sponge.
#
# Nondimensional: U = 1, lengths in the same units as `text_width`.
#
# Usage: julia --project=. scripts/flow_past_letters.jl [gpu]

include(joinpath(@__DIR__, "..", "src", "TextMask.jl"))
using .TextMask
using Oceananigans
using Oceananigans.Solvers: ConjugateGradientPoissonSolver
using Printf, Random

Random.seed!(1234)

arch = "gpu" in ARGS ? GPU() : CPU()

# ── Parameters ───────────────────────────────────────────────────────────────
text        = "UMD Ocean Dynamics"
Lx, Ly      = 128, 48
text_width  = 100              # stroke width ≈ 1.2
text_y      = -Ly / 4          # text center (wakes develop above it)
Nx          = arch isa GPU ? 4096 : 1024
Ny          = round(Int, Nx * Ly / Lx)

U           = 1
ν           = 2e-3             # Re ≈ U * stroke / ν ≈ 600
sponge_y    = -Ly / 2 + 3      # center of the restoring band
sponge_w    = 1.5              # Gaussian half-width of the band
sponge_rate = 2                # restoring rate (U / length)

stop_time   = 200              # ≈ 4 flow-through times
save_every  = 0.5

name = @sprintf("flow_letters_%dx%d", Nx, Ny)

# ── Grid with the text as an immersed boundary ───────────────────────────────
underlying_grid = RectilinearGrid(arch; size = (Nx, Ny), x = (-Lx/2, Lx/2), y = (-Ly/2, Ly/2),
                                  topology = (Periodic, Periodic, Flat), halo = (5, 5))

xc, yc = xnodes(underlying_grid, Center()), ynodes(underlying_grid, Center())
M = text_mask(text, Array(xc), Array(yc); width = text_width, center = (0, text_y))
solid = reshape(M .> 0.5, Nx, Ny, 1)

grid = ImmersedBoundaryGrid(underlying_grid, GridFittedBoundary(solid))

# ── Model ────────────────────────────────────────────────────────────────────
sponge(x, y, args...) = exp(-(y - sponge_y)^2 / (2sponge_w^2))
u_forcing = Relaxation(rate = sponge_rate, mask = sponge, target = 0)
v_forcing = Relaxation(rate = sponge_rate, mask = sponge, target = U)

no_slip = ValueBoundaryCondition(0)
boundary_conditions = (u = FieldBoundaryConditions(immersed = no_slip),
                       v = FieldBoundaryConditions(immersed = no_slip))

pressure_solver = ConjugateGradientPoissonSolver(grid; maxiter = 100)

model = NonhydrostaticModel(grid; advection = WENO(order = 5),
                            closure = ScalarDiffusivity(; ν),
                            forcing = (u = u_forcing, v = v_forcing),
                            boundary_conditions, pressure_solver)

vᵢ(x, y) = U * (1 + 1e-2 * randn())
uᵢ(x, y) = 1e-2 * U * randn()
set!(model, u = uᵢ, v = vᵢ)

# ── Output ───────────────────────────────────────────────────────────────────
u, v, w = model.velocities
ω = Field(∂x(v) - ∂y(u))

simulation = Simulation(model; Δt = 0.1 * (Lx / Nx) / U, stop_time)
wizard = TimeStepWizard(cfl = 0.7, max_change = 1.1)
simulation.callbacks[:wizard] = Callback(wizard, IterationInterval(10))

mkpath("output")
simulation.output_writers[:fields] = JLD2Writer(model, (; ω);
                                                filename = joinpath("output", name * ".jld2"),
                                                schedule = TimeInterval(save_every),
                                                overwrite_files = true)

wall = Ref(time_ns())
function progress(sim)
    elapsed = (time_ns() - wall[]) / 1e9
    @printf("t = %6.1f  iter = %6d  Δt = %.2e  max|u| = %.2f  max|v| = %.2f  wall = %.1f s\n",
            time(sim), iteration(sim), sim.Δt, maximum(abs, u), maximum(abs, v), elapsed)
    flush(stdout)
    wall[] = time_ns()
end
add_callback!(simulation, progress, TimeInterval(5))

Oceananigans.Diagnostics.erroring_NaNChecker!(simulation)
run!(simulation)
