# Variants of the retro and sunset letter styles over a navy |ζ/f| background (plotting only).
# Usage: julia --project=. scripts/letter_style_variants.jl [name] [frame]
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

r = 4; cells_per_stroke = 10
vibes = FTFont(fontfile("bubblevibes/BubbleVibes-Regular.otf"))
Mf = logo_mask(two_line(vibes); cells_per_stroke = r * cells_per_stroke, grid = (r * Nx, r * Ny))
xfine = collect(range(xf[1] + Δx / (2r), step = Δx / r, length = r * Nx))
yfine = collect(range(yf[1] + Δx / (2r), step = Δx / r, length = r * Ny))
navy = cgrad([RGBf(0.02, 0.06, 0.18), RGBf(0.10, 0.30, 0.55), RGBf(0.75, 0.88, 0.95)])

styles = [(:retro, "retro (as before): puffy shading"), (:retro_flat, "retro, flatter print-style shading"),
          (:sunset, "sunset (as before): one gradient over the whole logo"), (:sunset_lines, "sunset, gradient within each line"),
          (:sunset_gloss, "sunset, per-line gradient + gloss"), (:terps, "UMD colors: red UMD, gold OCEAN DYNAMICS")]
fig = Figure(size = (1800, 1100))
for (k, (style, label)) in enumerate(styles)
    ax = Axis(fig[(k + 1) ÷ 2, 2 - k % 2]; aspect = DataAspect(), title = label)
    hidedecorations!(ax); hidespines!(ax)
    heatmap!(ax, xf, yf, abs.(ζ); colormap = navy, colorrange = (0, 2))
    draw_styled_letters!(ax, xfine, yfine, Mf; style, stroke_px = r * cells_per_stroke)
end
resize_to_layout!(fig)
save("output/letter_style_variants.png", fig)
@info "Saved output/letter_style_variants.png"
