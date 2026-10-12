# Preview logo size vs. domain for the Bubble Vibes two-line layout (no simulation).
# Grid spacing is fixed at Δx = Ld / 5 = 2 km.
#
# Usage: julia --project=. scripts/logo_size_options.jl
using CairoMakie, FreeTypeAbstraction, Printf
include(joinpath(@__DIR__, "..", "src", "Logo.jl"))
include(joinpath(@__DIR__, "..", "src", "plotting.jl"))

cells_per_Ld = 5
lines = two_line(FTFont(fontfile("bubblevibes/BubbleVibes-Regular.otf")))
B = layout_bitmap(lines)
stroke_px = minimum(median_stroke(layout_bitmap([l])) for l in lines)

# Place the logo at `cells_per_stroke` in a grid: either with margins (cells) or into a fixed Nx × Ny
function place(cells_per_stroke; margins = nothing, grid = nothing)
    scale = stroke_px / cells_per_stroke
    isnothing(grid) && return rasterize(B, scale; margins)
    nx, ny = floor.(Int, size(B) ./ scale)
    return rasterize(B, scale; margins = ((grid[1] - nx) ÷ 2, (grid[2] - ny) ÷ 2))
end

options = [
    ("Current: Ld = ½ stroke, margins 16 Ld (x) / 9 Ld (y)",        10, (margins = (80, 45),)),
    ("Same 544 × 192 grid, logo ×1.2 (Ld ≈ 0.42 stroke)",            12, (grid = (544, 192),)),
    ("Same 544 × 192 grid, logo ×1.35 (Ld ≈ 0.37 stroke)",           13.5, (grid = (544, 192),)),
    ("Ld = ½ stroke, margins 8 Ld / 4.5 Ld (smaller domain)",        10, (margins = (40, 22),)),
    ("Ld = ½ stroke, margins 5 Ld / 3 Ld (smaller domain)",          10, (margins = (25, 15),)),
    ("Ld = ⅓ stroke, margins 8 Ld / 4.5 Ld",                         15, (margins = (40, 22),)),
]

ocean = RGBf(0.05, 0.12, 0.25)
fig = Figure(size = (1800, 1500))
for (k, (label, cps, kw)) in enumerate(options)
    M = place(cps; kw...)
    Nx, Ny = size(M)
    cost = Nx * Ny / (544 * 192)
    I = findall(M)
    fx = (maximum(i[1] for i in I) - minimum(i[1] for i in I) + 1) / Nx
    fy = (maximum(i[2] for i in I) - minimum(i[2] for i in I) + 1) / Ny
    ax = Axis(fig[(k + 1) ÷ 2, 2 - k % 2]; aspect = DataAspect(), backgroundcolor = ocean,
              title = @sprintf("%s\ngrid %d × %d = %.0f × %.0f km, %.2f× cost of last run\nlogo fills %.0f%% of width, %.0f%% of height",
                               label, Nx, Ny, 2Nx, 2Ny, cost, 100fx, 100fy))
    hidedecorations!(ax); hidespines!(ax)
    draw_letters!(ax, 1:Nx, 1:Ny, M; outline = 1.5)
    limits!(ax, 0.5, Nx + 0.5, 0.5, Ny + 0.5)
    @info @sprintf("%-60s %4d × %4d  cost %.2f×  fill %.0f%% × %.0f%%", label, Nx, Ny, cost, 100fx, 100fy)
end
resize_to_layout!(fig)
save("output/logo_size_options.png", fig)
@info "Saved output/logo_size_options.png"
