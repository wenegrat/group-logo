# Preview the retro and sunset letter styles over different background fields (plotting only).
# Usage: julia --project=. scripts/letter_backgrounds.jl [name] [frame]
using Oceananigans, CairoMakie, Printf
include(joinpath(@__DIR__, "..", "src", "Logo.jl"))
include(joinpath(@__DIR__, "..", "src", "letter_styles.jl"))

name = length(ARGS) ≥ 1 ? ARGS[1] : "eady_vibes_448x128x8"
fs = "output/$(name)_surface.jld2"
ζs = FieldTimeSeries(fs, "ζ"; backend = OnDisk()); bs = FieldTimeSeries(fs, "b"; backend = OnDisk())
us = FieldTimeSeries(fs, "u"; backend = OnDisk()); vs = FieldTimeSeries(fs, "v"; backend = OnDisk())
n = length(ARGS) ≥ 2 ? parse(Int, ARGS[2]) : length(ζs.times)
grid = ζs.grid; ug = grid.underlying_grid
Nx, Ny = size(grid, 1), size(grid, 2)
Δx = xspacings(ug, Center())[1] / 1e3
xf = xnodes(ug, Face())[1:Nx] ./ 1e3; yf = ynodes(ug, Face())[1:Ny] ./ 1e3
xc = xnodes(ug, Center()) ./ 1e3; yc = ynodes(ug, Center()) ./ 1e3

ζ = interior(ζs[n], 1:Nx, 1:Ny, 1) ./ 1e-4
b = interior(bs[n], 1:Nx, 1:Ny, 1) .* 1e3
ui = interior(us[n], 1:Nx, 1:Ny, 1); vi = interior(vs[n], 1:Nx, 1:Ny, 1)
uc = (ui .+ circshift(ui, (-1, 0))) ./ 2
vc = size(interior(vs[n]), 2) == Ny ? (vi .+ circshift(vi, (0, -1))) ./ 2 : (vi .+ interior(vs[n], 1:Nx, 2:Ny+1, 1)) ./ 2
speed = sqrt.(uc .^ 2 .+ vc .^ 2)

r = 4; cells_per_stroke = 10
vibes = FTFont(fontfile("bubblevibes/BubbleVibes-Regular.otf"))
Mf = logo_mask(two_line(vibes); cells_per_stroke = r * cells_per_stroke, grid = (r * Nx, r * Ny))
xfine = collect(range(xf[1] + Δx / (2r), step = Δx / r, length = r * Nx))
yfine = collect(range(yf[1] + Δx / (2r), step = Δx / r, length = r * Ny))

navy = cgrad([RGBf(0.02, 0.06, 0.18), RGBf(0.10, 0.30, 0.55), RGBf(0.75, 0.88, 0.95)])
backgrounds = [
    ("vorticity ζ/f, red/blue",         xf, yf, ζ,     :balance,     (-2, 2)),
    ("vorticity |ζ/f|, navy",           xf, yf, abs.(ζ), navy,       (0, 2)),
    ("surface buoyancy, deep blues",    xc, yc, b,     Reverse(:deep), extrema(b)),
    ("surface speed, ice",              xc, yc, speed, :ice,         (0, 1.4)),
]

fig = Figure(size = (1800, 1500))
for (row, (label, x, y, A, colormap, colorrange)) in enumerate(backgrounds), (col, style) in enumerate((:retro, :sunset))
    ax = Axis(fig[row, col]; aspect = DataAspect(), title = "$(style) on $(label)")
    hidedecorations!(ax); hidespines!(ax)
    heatmap!(ax, x, y, A; colormap, colorrange)
    draw_styled_letters!(ax, xfine, yfine, Mf; style, stroke_px = r * cells_per_stroke)
end
resize_to_layout!(fig)
save("output/letter_backgrounds.png", fig)
@info "Saved output/letter_backgrounds.png"
