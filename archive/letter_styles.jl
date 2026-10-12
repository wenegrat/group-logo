# Preview decorative letter styles over a saved surface-vorticity frame (plotting only).
# Usage: julia --project=. scripts/letter_styles.jl [name] [frame]
using Oceananigans, CairoMakie, Printf
include(joinpath(@__DIR__, "..", "src", "Logo.jl"))
include(joinpath(@__DIR__, "..", "src", "letter_styles.jl"))

name = length(ARGS) ≥ 1 ? ARGS[1] : "eady_vibes_448x128x8"
ζs = FieldTimeSeries("output/$(name)_surface.jld2", "ζ"; backend = OnDisk())
n = length(ARGS) ≥ 2 ? parse(Int, ARGS[2]) : length(ζs.times)
grid = ζs.grid; ug = grid.underlying_grid
Nx, Ny = size(grid, 1), size(grid, 2)
Δx = xspacings(ug, Center())[1] / 1e3
xf = xnodes(ug, Face())[1:Nx] ./ 1e3; yf = ynodes(ug, Face())[1:Ny] ./ 1e3
ζ = interior(ζs[n], 1:Nx, 1:Ny, 1) ./ 1e-4

# Fine mask (4× the grid resolution) of the same layout, aligned with the simulation mask
land = Array(grid.immersed_boundary.bottom_height[1:Nx, 1:Ny, 1]) .== 0
r = 4
vibes = FTFont(fontfile("bubblevibes/BubbleVibes-Regular.otf"))
cells_per_stroke = 10                                  # as in the runs (Ld = stroke / 2, 5 cells per Ld)
Mf = logo_mask(two_line(vibes); cells_per_stroke = r * cells_per_stroke, grid = (r * Nx, r * Ny))
xfine = range(xf[1] + Δx / (2r), step = Δx / r, length = r * Nx)
yfine = range(yf[1] + Δx / (2r), step = Δx / r, length = r * Ny)
@info "fine mask land cells / (r² × grid land cells) = $(round(sum(Mf) / (r^2 * sum(land)), digits = 3))"

styles = [(:plain, "current: white + black outline"), (:retro, "retro pastel + drop shadow"),
          (:candy, "candy colors, glossy"), (:bubblegum, "bubblegum pink, glossy"),
          (:sunset, "sunset gradient + drop shadow"), (:aqua, "aqua gradient, glossy + sparkles")]

fig = Figure(size = (1800, 1100))
for (k, (style, label)) in enumerate(styles)
    ax = Axis(fig[(k + 1) ÷ 2, 2 - k % 2]; aspect = DataAspect(), title = label)
    hidedecorations!(ax); hidespines!(ax)
    heatmap!(ax, xf, yf, ζ; colormap = :balance, colorrange = (-2, 2))
    draw_styled_letters!(ax, collect(xfine), collect(yfine), Mf; style, stroke_px = r * cells_per_stroke)
end
resize_to_layout!(fig)
save("output/letter_styles.png", fig)
@info "Saved output/letter_styles.png"
