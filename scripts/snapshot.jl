# Snapshot figure of one saved frame: surface vorticity and buoyancy, bottom vorticity, surface
# speed, and zooms. Emblem runs get the styled logo.
#
# Usage: julia --project=. scripts/snapshot.jl [name] [frame index; default = last] [dir]
#   LOGO_SCHEME = red (default) | white | red_rim | two_tone
#   ZERO_LAKES=1  show water trapped inside the logo as unheated (plotting only)
using Oceananigans, CairoMakie, Printf
include(joinpath(@__DIR__, "..", "src", "plotting.jl"))
include(joinpath(@__DIR__, "..", "src", "logo_config.jl"))
include(joinpath(@__DIR__, "..", "src", "letter_styles.jl"))
include(joinpath(@__DIR__, "..", "src", "emblem_overlay.jl"))

name = length(ARGS) ≥ 1 ? ARGS[1] : "heated_emblem_467x112x8"
dir = length(ARGS) ≥ 3 ? ARGS[3] : "output"
fs, fb = joinpath(dir, "$(name)_surface.jld2"), joinpath(dir, "$(name)_bottom.jld2")
f = 1e-4

ζs = FieldTimeSeries(fs, "ζ"; backend = OnDisk())
us = FieldTimeSeries(fs, "u"; backend = OnDisk())
vs = FieldTimeSeries(fs, "v"; backend = OnDisk())
ζb = FieldTimeSeries(fb, "ζ"; backend = OnDisk())
bs = FieldTimeSeries(fs, "b"; backend = OnDisk())

n = length(ARGS) ≥ 2 ? parse(Int, ARGS[2]) : length(ζs.times)
t = ζs.times[n] / 86400

grid = ζs.grid; ug = grid.underlying_grid
Nx, Ny = size(grid, 1), size(grid, 2)
xc = xnodes(ug, Center()) ./ 1e3; yc = ynodes(ug, Center()) ./ 1e3
xf = xnodes(ug, Face())[1:Nx] ./ 1e3; yf = ynodes(ug, Face())[1:Ny] ./ 1e3

land = Array(grid.immersed_boundary.bottom_height[1:Nx, 1:Ny, 1]) .== 0

# Runs built from a logo layout (src/logo_config.jl) get the styled logo.
# LOGO_SCHEME: text layout white (default) | yellow; emblem layout red (default) | white | red_rim | two_tone
layout = layout_from_name(name)
emblem_run = !isnothing(layout)
if emblem_run
    Δx_km = xf[2] - xf[1]
    ED = emblem_overlay_data(; cells_per_Ld = 10 / Δx_km, Δx_km, x₀_km = xf[1], y₀_km = yf[1], layout)
    scheme = Symbol(get(ENV, "LOGO_SCHEME", string(default_scheme(layout))))
end

ui = interior(us[n], 1:Nx, 1:Ny, 1); vi = interior(vs[n], 1:Nx, 1:Ny, 1)
uc = (ui .+ circshift(ui, (-1, 0))) ./ 2                       # to cell centers (periodic x)
vc = size(interior(vs[n]), 2) == Ny ? (vi .+ circshift(vi, (0, -1))) ./ 2 :   # periodic y
                                             (vi .+ interior(vs[n], 1:Nx, 2:Ny+1, 1)) ./ 2
speed = sqrt.(uc .^ 2 .+ vc .^ 2)
ζsurf = interior(ζs[n], 1:Nx, 1:Ny, 1) ./ f
ζbot  = interior(ζb[n], 1:Nx, 1:Ny, 1) ./ f
heated = startswith(name, "heated")
if heated   # buoyancy anomaly relative to the initial state; color range from the open-ocean 99th percentile
    bsurf = (interior(bs[n], 1:Nx, 1:Ny, 1) .- interior(bs[1], 1:Nx, 1:Ny, 1)) .* 1e3
    # ZERO_LAKES=1: show water trapped inside letters as unheated (b′ = 0); plotting only
    zero_lakes = get(ENV, "ZERO_LAKES", "0") == "1"
    zero_lakes && (bsurf[enclosed_water(land; periodic = (true, size(interior(vs[n]), 2) == Ny))] .= 0)
    bvals = sort(bsurf[.!land]); blims = (0, bvals[ceil(Int, 0.99 * length(bvals))])
else
    bsurf = interior(bs[n], 1:Nx, 1:Ny, 1) .* 1e3
    blims = extrema(interior(bs[1], 1:Nx, 1:Ny, 1)[.!land]) .* 1e3   # fixed to the initial range
end

fig = Figure(size = (1800, 1850))
Label(fig[0, 1:2], @sprintf("%s, t = %.2f days", heated ? "Heated logo islands" : "Eady channel with logo islands", t), fontsize = 24)

function panel!(pos, x, y, A, title; colormap, colorrange, xlims = nothing, ylims = nothing)
    ax = Axis(pos[1, 1]; title, aspect = DataAspect())
    hidedecorations!(ax)
    hm = heatmap!(ax, x, y, A; colormap, colorrange)
    emblem_run ? draw_emblem_logo!(ax, ED; scheme) : draw_letters!(ax, xc, yc, land)
    isnothing(xlims) || xlims!(ax, xlims...)
    isnothing(ylims) || ylims!(ax, ylims...)
    return ax
end

ζlim = 2.0
panel!(fig[1, 1:2], xf, yf, ζsurf, "surface vorticity"; colormap = :balance, colorrange = (-ζlim, ζlim))
panel!(fig[2, 1:2], xc, yc, bsurf, heated ? "surface buoyancy anomaly" : "surface buoyancy";
       colormap = :thermal, colorrange = blims)
panel!(fig[3, 1:2], xf, yf, ζbot,  "bottom vorticity";  colormap = :balance, colorrange = (-ζlim, ζlim))
panel!(fig[4, 1:2], xc, yc, speed, "surface speed"; colormap = :viridis, colorrange = (0, 1.4))

# Zoom on the top-center of the logo ("UMD" in the two-line layout), padded by 15 km
I = findall(land); xl = extrema(xc[i[1]] for i in I); yl = extrema(yc[i[2]] for i in I)
xm, w = (xl[1] + xl[2]) / 2, xl[2] - xl[1]
zoom = ((xm - 0.2w, xm + 0.2w), (yl[1] + 0.4 * (yl[2] - yl[1]), yl[2] + 15))
if emblem_run   # emblem plus the start of the text, padded by 25 km
    bx, by, bw, bh = ED.box
    zoom = ((bx - 25, bx + 3.2bw), (yl[1] - 25, yl[2] + 25))
end
panel!(fig[5, 1], xf, yf, ζsurf, "zoom: surface vorticity"; colormap = :balance,
       colorrange = (-ζlim, ζlim), xlims = zoom[1], ylims = zoom[2])
panel!(fig[5, 2], xc, yc, speed, "zoom: surface speed"; colormap = :viridis,
       colorrange = (0, 1.4), xlims = zoom[1], ylims = zoom[2])

out = @sprintf("output/%s_snapshot_t%05.2fd%s.png", name, t, (heated && zero_lakes) ? "_nolakes" : "")
save(out, fig)
@info "Saved $out"
