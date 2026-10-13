# Render frames for a web header: one panel of the flow with the styled logo, at a fixed pixel width.
# Modes (see `style` below): white, redblue, aqua, signed, signedD, bfoam, opaque.
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
layout = layout_from_name(name)
ED = emblem_overlay_data(; cells_per_Ld = 10 / Δx_km, Δx_km, x₀_km = xf[1], y₀_km = yf[1], layout)

# Edge taper (transparent modes): the field fades to transparent over the outer `taper` fraction of
# the width and height, so the banner blends into the page instead of ending in a hard rectangle.
taper = (parse(Float64, get(ENV, "TAPER_X", "0.15")), parse(Float64, get(ENV, "TAPER_Y", "0.35")))
smooth(s) = (s = clamp(s, 0, 1); s^2 * (3 - 2s))                    # smoothstep
ramp(k, N, w) = w ≤ 0 ? 1.0 : smooth(min(k - 0.5, N - k + 0.5) / (w * N))
edge = [ramp(i, Nx, taper[1]) * ramp(j, Ny, taper[2]) for i in 1:Nx, j in 1:Ny]
coastfade = smooth.((distance_to(land) .- 1) ./ 4)               # fades out the band hugging the coast

aqua, foam, peach = RGB(0.62, 0.95, 0.92), RGB(0.85, 1.0, 0.95), RGB(1.0, 0.82, 0.62)
deep_aqua, coral, pale = RGB(0.10, 0.75, 0.78), RGB(0.98, 0.55, 0.40), RGB(0.92, 0.97, 0.96)
lerp(a, b, t) = RGB(a.r + (b.r - a.r) * t, a.g + (b.g - a.g) * t, a.b + (b.b - a.b) * t)
# Signed vorticity, nearly opaque except near ζ = 0 (banner_styles.jl, BANNER_SET=signed, option D)
signed_opaque(z; z₀ = 0.05, z₁ = 0.3, amax = 0.9) =
    (lerp(pale, z > 0 ? deep_aqua : coral, clamp(abs(z) / 1.5, 0, 1)), amax * smooth((abs(z) - z₀) / (z₁ - z₀)))
cmap = cgrad(:balance)

# Per mode: (field, value → (color, alpha), fade the coastal band?)
#   white      white, alpha ∝ |ζ|/f                      (original)
#   redblue    red/blue vorticity, alpha ∝ |ζ|/f
#   aqua       aqua tint, alpha ∝ |ζ|/f, coast faded      (banner_styles option 3)
#   signed     cyclones aqua, anticyclones peach           (option 5)
#   signedD    signed ζ, opaque except near 0: cyclones deep aqua, anticyclones coral
#   bfoam      buoyancy anomaly as a sea-foam glow         (option 8)
#   opaque     red/blue vorticity on an opaque background
ramp_ζ(z; s = 1.2, p = 1.2) = clamp(abs(z) / s, 0, 1)^p
# bfoam opacity: maximum, and ramp exponent (smaller = more opaque at intermediate warming)
bfoam_max = parse(Float64, get(ENV, "BFOAM_MAX", "0.9"))
bfoam_p   = parse(Float64, get(ENV, "BFOAM_P", "1.0"))
style = Dict(
    :white   => (:ζ, z -> (RGB(1, 1, 1), 0.85ramp_ζ(z)), false),
    :redblue => (:ζ, z -> (cmap[clamp((z / 2 + 1) / 2, 0, 1)], ramp_ζ(z)), false),
    :aqua    => (:ζ, z -> (aqua, 0.9ramp_ζ(z)), true),
    :signed  => (:ζ, z -> (z > 0 ? aqua : peach, 0.9ramp_ζ(z)), true),
    :signedD => (:ζ, signed_opaque, true),
    :bfoam   => (:b, x -> (foam, bfoam_max * ramp_ζ(x; s = 1.0, p = bfoam_p)), false),
    :opaque  => (:ζ, z -> (cmap[clamp((z / 2 + 1) / 2, 0, 1)], 1.0), false))
fieldname, colorfun, coast = style[mode]

if fieldname == :b      # buoyancy anomaly, scaled by the 98th percentile of open water at the last frame
    bs = FieldTimeSeries(joinpath("output", "$(name)_surface.jld2"), "b"; backend = OnDisk())
    b₀ = interior(bs[1], 1:Nx, 1:Ny, 1)
    blast = sort((interior(bs[length(bs.times)], 1:Nx, 1:Ny, 1) .- b₀)[.!land .& .!enclosed_water(land; periodic = (true, true))])
    bscale = blast[ceil(Int, 0.98length(blast))]
end
field(n) = fieldname == :ζ ? interior(ζs[n], 1:Nx, 1:Ny, 1) ./ 1e-4 :
                             (interior(bs[n], 1:Nx, 1:Ny, 1) .- b₀) ./ bscale

# Water sealed inside letters (e.g. the counters of O, D, A) is heated with no escape; leave it out
lakes = enclosed_water(land; periodic = (true, true))

function frame(n)
    F = field(n)
    fieldname == :b && (F = max.(F, 0))                  # glow only where the water has warmed
    # Transparent pixels keep the overlay's color (not black), so smooth interpolation at the edges of
    # the glow doesn't blend toward black and leave dark smudges.
    c₀, _ = colorfun(0.0)
    O = fill(RGBA(c₀.r, c₀.g, c₀.b, 0.0), Nx, Ny)
    for j in 1:Ny, i in 1:Nx
        (land[i, j] || lakes[i, j]) && continue
        c, a = colorfun(F[i, j])
        mode == :opaque || (a *= edge[i, j] * (coast ? coastfade[i, j] : 1.0))
        O[i, j] = RGBA(c.r, c.g, c.b, a)
    end
    return O
end

H = round(Int, W * Ny / Nx)
fig = Figure(size = (W, H), backgroundcolor = mode == :opaque ? :white : :transparent, figure_padding = 0)
ax = Axis(fig[1, 1]; backgroundcolor = :transparent)
hidedecorations!(ax); hidespines!(ax)
img = Observable(frame(1))
image!(ax, (xf[1], xf[end] + Δx_km), (yf[1], yf[end] + Δx_km), img; interpolate = true)
draw_emblem_logo!(ax, ED; scheme = Symbol(get(ENV, "LOGO_SCHEME", string(default_scheme(layout)))))
limits!(ax, xf[1], xf[end] + Δx_km, yf[1], yf[end] + Δx_km)

dir = joinpath("output", "web_frames", string(mode)); mkpath(dir)
for n in eachindex(ζs.times)
    img[] = frame(n)
    save(joinpath(dir, @sprintf("frame_%03d.png", n)), fig; px_per_unit = 1)
end
@info "Saved $(length(ζs.times)) frames ($W × $H) to $dir"
