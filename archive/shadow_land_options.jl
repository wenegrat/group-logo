# Preview a drop shadow made part of the land (immersed boundary), on the 448 × 154 Eady grid.
# Left: the land mask the model would see (grid resolution). Right: the same layout drawn retro-style.
# Usage: julia --project=. scripts/shadow_land_options.jl
using CairoMakie, Printf
include(joinpath(@__DIR__, "..", "src", "Logo.jl"))
include(joinpath(@__DIR__, "..", "src", "plotting.jl"))
include(joinpath(@__DIR__, "..", "src", "letter_styles.jl"))

vibes = FTFont(fontfile("bubblevibes/BubbleVibes-Regular.otf"))
lines = two_line(vibes)
grid, cps, r = (448, 154), 10, 4
navy = RGBf(0.05, 0.12, 0.25)

options = [("no shadow", nothing), ("shadow 0.25 stroke", (0.25, -0.25)), ("shadow 0.35 stroke", (0.35, -0.35))]
fig = Figure(size = (1800, 1300))
for (row, (label, sh)) in enumerate(options)
    land = logo_mask(lines; cells_per_stroke = cps, grid, shadow = sh)
    lakes = enclosed_water(land; periodic = (true, false))
    stats = @sprintf("land %.1f%% of domain, %d trapped-water cells", 100sum(land) / length(land), sum(lakes))
    @info "$label: $stats"

    ax = Axis(fig[row, 1]; aspect = DataAspect(), backgroundcolor = navy, title = "$label: model land mask\n$stats")
    hidedecorations!(ax); hidespines!(ax)
    heatmap!(ax, 1:grid[1], 1:grid[2], ifelse.(land, 1.0, NaN); colormap = [:white], colorrange = (0, 1))
    limits!(ax, 0.5, grid[1] + 0.5, 0.5, grid[2] + 0.5)

    finegrid = r .* grid
    Lf = logo_mask(lines; cells_per_stroke = r * cps, grid = finegrid, shadow = sh, letters_only = !isnothing(sh))
    Sf = isnothing(sh) ? falses(finegrid) : logo_mask(lines; cells_per_stroke = r * cps, grid = finegrid, shadow = sh)
    xs = collect(range(0.5 + 0.5 / r, step = 1 / r, length = finegrid[1]))
    ys = collect(range(0.5 + 0.5 / r, step = 1 / r, length = finegrid[2]))
    ax = Axis(fig[row, 2]; aspect = DataAspect(), backgroundcolor = navy, title = "$label: drawn (retro)")
    hidedecorations!(ax); hidespines!(ax)
    draw_styled_letters!(ax, xs, ys, Lf; style = :retro, stroke_px = r * cps, shadow = isnothing(sh) ? falses(finegrid) : Sf)
    limits!(ax, 0.5, grid[1] + 0.5, 0.5, grid[2] + 0.5)
end
resize_to_layout!(fig)
save("output/shadow_land_options.png", fig)
@info "Saved output/shadow_land_options.png"
