# Animate a run: surface vorticity, surface buoyancy (anomaly for heated runs), bottom vorticity,
# surface speed. Emblem runs get the styled logo; other runs show land as white with a black outline.
#
# Usage: julia --project=. scripts/animate.jl [name] [dir]
#   LOGO_SCHEME = red (default) | white | red_rim | two_tone   (colors of OCEAN DYNAMICS)
using Oceananigans, CairoMakie, Printf
include(joinpath(@__DIR__, "..", "src", "plotting.jl"))
include(joinpath(@__DIR__, "..", "src", "logo_config.jl"))
include(joinpath(@__DIR__, "..", "src", "letter_styles.jl"))
include(joinpath(@__DIR__, "..", "src", "emblem_overlay.jl"))

name = length(ARGS) ≥ 1 ? ARGS[1] : "heated_emblem_467x112x8"
dir  = length(ARGS) ≥ 2 ? ARGS[2] : "output"
f = 1e-4

fs, fb = joinpath(dir, "$(name)_surface.jld2"), joinpath(dir, "$(name)_bottom.jld2")
ζs = FieldTimeSeries(fs, "ζ"; backend = OnDisk())
bs = FieldTimeSeries(fs, "b"; backend = OnDisk())
us = FieldTimeSeries(fs, "u"; backend = OnDisk())
vs = FieldTimeSeries(fs, "v"; backend = OnDisk())
ζb = FieldTimeSeries(fb, "ζ"; backend = OnDisk())
times = ζs.times[1:min(length(ζs.times), length(ζb.times))]

grid = ζs.grid; ug = grid.underlying_grid
Nx, Ny = size(grid, 1), size(grid, 2)
xc = xnodes(ug, Center()) ./ 1e3; yc = ynodes(ug, Center()) ./ 1e3
xf = xnodes(ug, Face())[1:Nx] ./ 1e3; yf = ynodes(ug, Face())[1:Ny] ./ 1e3
land = Array(grid.immersed_boundary.bottom_height[1:Nx, 1:Ny, 1]) .== 0

# Runs with the emblem layout get the styled logo (scheme from ENV LOGO_SCHEME: red, white, red_rim, two_tone)
emblem_run = occursin("emblem", name)
if emblem_run
    Δx_km = xf[2] - xf[1]
    ED = emblem_overlay_data(; cells_per_Ld = 10 / Δx_km, Δx_km, x₀_km = xf[1], y₀_km = yf[1])
    scheme = Symbol(get(ENV, "LOGO_SCHEME", "red"))
end

function speed(n)
    ui = interior(us[n], 1:Nx, 1:Ny, 1); vi = interior(vs[n], 1:Nx, 1:Ny, 1)
    uc = (ui .+ circshift(ui, (-1, 0))) ./ 2
    vc = size(interior(vs[n]), 2) == Ny ? (vi .+ circshift(vi, (0, -1))) ./ 2 :   # periodic y
                                             (vi .+ interior(vs[n], 1:Nx, 2:Ny+1, 1)) ./ 2
    return sqrt.(uc .^ 2 .+ vc .^ 2)
end

n = Observable(1)
ζsurf = @lift interior(ζs[$n], 1:Nx, 1:Ny, 1) ./ f
heated = startswith(name, "heated")
b₀ = interior(bs[1], 1:Nx, 1:Ny, 1)
bsurf = heated ? (@lift (interior(bs[$n], 1:Nx, 1:Ny, 1) .- b₀) .* 1e3) :   # anomaly from the initial state
                 (@lift interior(bs[$n], 1:Nx, 1:Ny, 1) .* 1e3)
ζbot  = @lift interior(ζb[$n], 1:Nx, 1:Ny, 1) ./ f
spd   = @lift speed($n)
blims = if heated   # 0 to the open-water 99th percentile of the final anomaly
    bend = sort(((interior(bs[length(times)], 1:Nx, 1:Ny, 1) .- b₀) .* 1e3)[.!land])
    (0, bend[ceil(Int, 0.99 * length(bend))])
else
    extrema(b₀[.!land]) .* 1e3
end
ζlim  = 2.0

aspect = (xf[end] - xf[1]) / (yf[end] - yf[1])
W = 1800
fig = Figure(size = (W, round(Int, 4 * W / aspect + 4 * 35 + 70)))
title = @lift @sprintf("%s, t = %.1f days", heated ? "Heated logo islands" : "Eady channel with logo islands", times[$n] / 86400)
Label(fig[0, 1], title, fontsize = 24, tellwidth = false)

panels = ((xf, yf, ζsurf, "surface vorticity", :balance, (-ζlim, ζlim)),
          (xc, yc, bsurf, heated ? "surface buoyancy anomaly" : "surface buoyancy", :thermal, blims),
          (xf, yf, ζbot,  "bottom vorticity", :balance, (-ζlim, ζlim)),
          (xc, yc, spd,   "surface speed", :viridis, (0, 1.4)))

for (row, (x, y, A, label, colormap, colorrange)) in enumerate(panels)
    ax = Axis(fig[row, 1]; title = label, aspect = DataAspect())
    hidedecorations!(ax)
    hm = heatmap!(ax, x, y, A; colormap, colorrange)
    emblem_run ? draw_emblem_logo!(ax, ED; scheme) : draw_letters!(ax, xc, yc, land)
end

out = joinpath("output", "$(name)_panels.mp4")
record(fig, out, eachindex(times); framerate = 12) do i
    n[] = i
end
@info "Saved $out ($(length(times)) frames, t = 0–$(round(times[end] / 86400, digits = 1)) days)"
