# Render frames for a web header: one panel of surface vorticity with the styled logo, at a fixed
# pixel width. Modes:
#   white    white filaments, alpha ∝ |ζ|/f (transparent where the flow is calm)
#   redblue  red/blue vorticity, alpha ∝ |ζ|/f
#   opaque   red/blue vorticity on an opaque background
# Frames go to output/web_frames/<mode>/frame_NNN.png (RGBA). In the transparent modes the field fades
# out toward the edges (TAPER_X, TAPER_Y: fraction of width/height, defaults 0.15 and 0.35).
#
# Usage: julia --project=. scripts/web_frames.jl [name] [mode] [width px]
using Oceananigans, CairoMakie, Printf
using CairoMakie.Colors: RGBA, RGB
include(joinpath(@__DIR__, "..", "src", "plotting.jl"))
include(joinpath(@__DIR__, "..", "src", "logo_config.jl"))
include(joinpath(@__DIR__, "..", "src", "letter_styles.jl"))
include(joinpath(@__DIR__, "..", "src", "emblem_overlay.jl"))

name = length(ARGS) ≥ 1 ? ARGS[1] : "heated_emblem_700x170x8"
mode = length(ARGS) ≥ 2 ? Symbol(ARGS[2]) : :white
W    = length(ARGS) ≥ 3 ? parse(Int, ARGS[3]) : 1400

ζs = FieldTimeSeries(joinpath("output", "$(name)_surface.jld2"), "ζ"; backend = OnDisk())
grid = ζs.grid; ug = grid.underlying_grid
Nx, Ny = size(grid, 1), size(grid, 2)
xf = xnodes(ug, Face())[1:Nx] ./ 1e3; yf = ynodes(ug, Face())[1:Ny] ./ 1e3
Δx_km = xf[2] - xf[1]
land = Array(grid.immersed_boundary.bottom_height[1:Nx, 1:Ny, 1]) .== 0
ED = emblem_overlay_data(; cells_per_Ld = 10 / Δx_km, Δx_km, x₀_km = xf[1], y₀_km = yf[1])

ζmax = 2.0                                       # ζ/f at full color
alpha(z) = clamp(abs(z) / 1.2, 0, 1)^1.2         # opacity ramp: transparent near 0, opaque by |ζ| ≈ 1.2 f
cmap = cgrad(:balance)
function colorize(ζ)
    map(ζ) do z
        isnan(z) && return RGBA(0, 0, 0, 0)
        c = cmap[clamp((z / ζmax + 1) / 2, 0, 1)]
        mode == :white   ? RGBA(1, 1, 1, 0.85alpha(z)) :
        mode == :redblue ? RGBA(c.r, c.g, c.b, alpha(z)) :
                           RGBA(c.r, c.g, c.b, 1)
    end
end
# Edge taper (transparent modes): the field fades to transparent over the outer `taper` fraction of
# the width and height, so the banner blends into the page instead of ending in a hard rectangle.
taper = (parse(Float64, get(ENV, "TAPER_X", "0.15")), parse(Float64, get(ENV, "TAPER_Y", "0.35")))
smooth(s) = (s = clamp(s, 0, 1); s^2 * (3 - 2s))                    # smoothstep
ramp(k, N, w) = w ≤ 0 ? 1.0 : smooth(min(k - 0.5, N - k + 0.5) / (w * N))
edge = [ramp(i, Nx, taper[1]) * ramp(j, Ny, taper[2]) for i in 1:Nx, j in 1:Ny]

function frame(n)
    ζ = interior(ζs[n], 1:Nx, 1:Ny, 1) ./ 1e-4
    ζ[land] .= NaN
    C = colorize(ζ)
    mode == :opaque && return C
    return [RGBA(c.r, c.g, c.b, c.alpha * e) for (c, e) in zip(C, edge)]
end

H = round(Int, W * Ny / Nx)
fig = Figure(size = (W, H), backgroundcolor = mode == :opaque ? :white : :transparent, figure_padding = 0)
ax = Axis(fig[1, 1]; backgroundcolor = :transparent)
hidedecorations!(ax); hidespines!(ax)
img = Observable(frame(1))
image!(ax, (xf[1], xf[end] + Δx_km), (yf[1], yf[end] + Δx_km), img; interpolate = true)
draw_emblem_logo!(ax, ED; scheme = Symbol(get(ENV, "LOGO_SCHEME", "red")))
limits!(ax, xf[1], xf[end] + Δx_km, yf[1], yf[end] + Δx_km)

dir = joinpath("output", "web_frames", string(mode)); mkpath(dir)
for n in eachindex(ζs.times)
    img[] = frame(n)
    save(joinpath(dir, @sprintf("frame_%03d.png", n)), fig; px_per_unit = 1)
end
@info "Saved $(length(ζs.times)) frames ($W × $H) to $dir"
