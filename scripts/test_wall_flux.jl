# Check the sign and magnitude of a buoyancy flux imposed on the lateral faces of a full-depth
# immersed island (GridFittedBottom with bottom at the surface), hydrostatic model at rest.
using Oceananigans, Oceananigans.Units, Printf

Nx, Ny, Nz = 32, 32, 4
Δ, H = 2kilometers, 1000meters
grid = RectilinearGrid(size = (Nx, Ny, Nz), x = (0, Nx * Δ), y = (0, Ny * Δ), z = (-H, 0),
                       topology = (Periodic, Periodic, Bounded), halo = (7, 7, 4))
land = [12 ≤ i ≤ 20 && 10 ≤ j ≤ 22 for i in 1:Nx, j in 1:Ny]          # 9 × 13 cell island
grid = ImmersedBoundaryGrid(grid, GridFittedBottom(ifelse.(land, 0.0, -H)))

Fb = 3e-5                                                             # m² s⁻³ (buoyancy flux per wall area)
for sgn in (+1, -1)
    q = FluxBoundaryCondition(sgn * Fb)
    b_bcs = FieldBoundaryConditions(immersed = ImmersedBoundaryCondition(west = q, east = q, south = q, north = q))
    model = HydrostaticFreeSurfaceModel(grid; buoyancy = BuoyancyTracer(), tracers = :b,
                                        free_surface = ImplicitFreeSurface(), coriolis = FPlane(f = 1e-4),
                                        boundary_conditions = (; b = b_bcs))
    B() = sum(interior(compute!(Field(Integral(model.tracers.b)))))
    B₀ = B()
    T = 1hour
    simulation = Simulation(model; Δt = 1minute, stop_time = T, verbose = false)
    run!(simulation)
    perimeter = 2 * (9 + 13) * Δ                                       # wall length in m
    expected = Fb * perimeter * H * T                                   # m³ · m s⁻² of buoyancy content
    @printf("FluxBoundaryCondition(%+.0e): ΔB = %+.4e, expected magnitude %.4e, ratio %+.3f\n",
            sgn * Fb, B() - B₀, expected, (B() - B₀) / expected)
end
