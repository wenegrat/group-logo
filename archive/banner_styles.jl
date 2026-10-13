# Compare ways of drawing the flow over the website's header photo (plotting only): one saved frame
# rendered with several fields, tints and opacity ramps, each composited over the photo with the logo.
#
# Usage: julia --project=. scripts/banner_styles.jl [name] [day] [background image]
#   BANNER_SET = default | signed | bright
using Oceananigans, CairoMakie, Printf, Statistics
using CairoMakie.Colors: RGBA, RGB
include(joinpath(@__DIR__, "..", "src", "plotting.jl"))
include(joinpath(@__DIR__, "..", "src", "logo_config.jl"))
include(joinpath(@__DIR__, "..", "src", "letter_styles.jl"))
include(joinpath(@__DIR__, "..", "src", "emblem_overlay.jl"))

name = length(ARGS) ≥ 1 ? ARGS[1] : "heated_emblem_700x170x8"
day  = length(ARGS) ≥ 2 ? parse(Float64, ARGS[2]) : 15.0
bg_path = length(ARGS) ≥ 3 ? ARGS[3] : joinpath(REPO, "output", "website_mock", "site", "img", "Landsat.jpg")

fs = joinpath("output", "$(name)_surface.jld2")
ζs = FieldTimeSeries(fs, "ζ"; backend = OnDisk()); bs = FieldTimeSeries(fs, "b"; backend = OnDisk())
n = argmin(abs.(ζs.times ./ 86400 .- day))
grid = ζs.grid; ug = grid.underlying_grid
Nx, Ny = size(grid, 1), size(grid, 2)
xf = xnodes(ug, Face())[1:Nx] ./ 1e3; yf = ynodes(ug, Face())[1:Ny] ./ 1e3
Δx_km = xf[2] - xf[1]
land = Array(grid.immersed_boundary.bottom_height[1:Nx, 1:Ny, 1]) .== 0
layout = layout_from_name(name)
ED = emblem_overlay_data(; cells_per_Ld = 10 / Δx_km, Δx_km, x₀_km = xf[1], y₀_km = yf[1], layout)

ζ = interior(ζs[n], 1:Nx, 1:Ny, 1) ./ 1e-4
b′ = interior(bs[n], 1:Nx, 1:Ny, 1) .- interior(bs[1], 1:Nx, 1:Ny, 1)
lakes = enclosed_water(land; periodic = (true, true))
bopen = sort(b′[.!land .& .!lakes]); bscale = bopen[ceil(Int, 0.98length(bopen))]

# Shared masks: edge taper and (optional) suppression of the band hugging the coast
smooth(s) = (s = clamp(s, 0, 1); s^2 * (3 - 2s))
ramp(k, N, w) = smooth(min(k - 0.5, N - k + 0.5) / (w * N))
edge = [ramp(i, Nx, 0.15) * ramp(j, Ny, 0.35) for i in 1:Nx, j in 1:Ny]
dcoast = distance_to(land)                                       # cells to the nearest land
coastfade = smooth.((dcoast .- 1) ./ 4)                          # 0 at the coast → 1 beyond ~5 cells

aqua, foam, peach, navy, amber = RGB(0.62, 0.95, 0.92), RGB(0.85, 1.0, 0.95), RGB(1.0, 0.82, 0.62), RGB(0.03, 0.10, 0.22), RGB(1.0, 0.70, 0.25)
unsigned(c; s = 1.2, p = 1.2, amax = 0.85) = z -> (c, amax * clamp(abs(z) / s, 0, 1)^p)
signed(cpos, cneg; s = 1.2, p = 1.2, amax = 0.85) = z -> (z > 0 ? cpos : cneg, amax * clamp(abs(z) / s, 0, 1)^p)

# Signed vorticity, nearly opaque except near ζ = 0: alpha rises from 0 at |ζ| = z₀ to `amax` at
# |ζ| = z₁; color goes from a pale tint at small |ζ| to the full color at |ζ| ≥ 1.5 f.
lerp(a, b, t) = RGB(a.r + (b.r - a.r) * t, a.g + (b.g - a.g) * t, a.b + (b.b - a.b) * t)
function signed_opaque(cpos, cneg; z₀ = 0.05, z₁ = 0.3, amax = 0.9, pale = RGB(0.92, 0.97, 0.96), zfull = 1.5)
    z -> begin
        t = clamp(abs(z) / zfull, 0, 1)
        (lerp(pale, z > 0 ? cpos : cneg, t), amax * smooth((abs(z) - z₀) / (z₁ - z₀)))
    end
end
deep_aqua, coral = RGB(0.10, 0.75, 0.78), RGB(0.98, 0.55, 0.40)

bright_aqua, bright_coral = RGB(0.15, 0.90, 0.92), RGB(1.0, 0.50, 0.35)
which_set = get(ENV, "BANNER_SET", "default")

# (label, field, color/alpha function, suppress coastal band?)
options = which_set == "signed" ? [
    ("A. signed, transparent only for |ζ| < 0.05 f (opaque by 0.3 f)", ζ, signed_opaque(aqua, peach), true),
    ("B. same, transparent for |ζ| < 0.1 f",                           ζ, signed_opaque(aqua, peach; z₀ = 0.1, z₁ = 0.4), true),
    ("C. same, max opacity 70%",                                       ζ, signed_opaque(aqua, peach; amax = 0.7), true),
    ("D. stronger colors (deep aqua / coral), |ζ| < 0.05 f clear",    ζ, signed_opaque(deep_aqua, coral), true),
] : which_set == "bright" ? [
    ("D (current)",                                                    ζ, signed_opaque(deep_aqua, coral), true),
    ("E. D at full opacity, full color by |ζ| = 1 f",                  ζ, signed_opaque(deep_aqua, coral; amax = 1.0, zfull = 1.0), true),
    ("F. brighter hues (vivid aqua / coral), full color by 1 f",       ζ, signed_opaque(bright_aqua, bright_coral; amax = 1.0, zfull = 1.0), true),
    ("G. F with full color by 0.6 f (most saturated)",                 ζ, signed_opaque(bright_aqua, bright_coral; amax = 1.0, zfull = 0.6), true),
] : [
    ("1. current: white, |ζ| ramp, max 85%",                         ζ,  unsigned(RGB(1, 1, 1)), false),
    ("2. white, steep ramp (strong filaments only), coast faded",    ζ,  unsigned(RGB(1, 1, 1); s = 1.5, p = 2.5, amax = 1.0), true),
    ("3. aqua tint, |ζ| ramp, coast faded",                          ζ,  unsigned(aqua; amax = 0.9), true),
    ("4. sea-foam, soft (max 55%), coast faded",                     ζ,  unsigned(foam; s = 1.0, p = 1.0, amax = 0.55), true),
    ("5. signed ζ: cyclones aqua, anticyclones peach",               ζ,  signed(aqua, peach; amax = 0.9), true),
    ("6. signed ζ: cyclones white, anticyclones navy",               ζ,  signed(RGB(1, 1, 1), navy; amax = 0.85), true),
    ("7. buoyancy anomaly: amber glow",                              b′ ./ bscale, unsigned(amber; s = 1.0, p = 1.5, amax = 0.75), false),
    ("8. buoyancy anomaly: sea-foam glow",                           b′ ./ bscale, unsigned(foam; s = 1.0, p = 1.5, amax = 0.7), false),
]

function overlay(field, f, coast)
    c₀, _ = f(0.0)
    O = fill(RGBA(c₀.r, c₀.g, c₀.b, 0.0), Nx, Ny)
    for j in 1:Ny, i in 1:Nx
        (land[i, j] || lakes[i, j]) && continue
        c, a = f(field[i, j])
        a *= edge[i, j] * (coast ? coastfade[i, j] : 1.0)
        O[i, j] = RGBA(c.r, c.g, c.b, a)
    end
    return O
end

bg = CairoMakie.FileIO.load(bg_path); bg = permutedims(reverse(bg, dims = 1))
bw = size(bg, 1); bh = min(size(bg, 2), round(Int, bw * Ny / Nx)); j₀ = (size(bg, 2) - bh) ÷ 2
X = (xf[1], xf[end] + Δx_km); Y = (yf[1], yf[end] + Δx_km)

fig = Figure(size = (1800, cld(length(options), 2) * round(Int, 880 * Ny / Nx) + 260))
for (k, (label, field, f, coast)) in enumerate(options)
    ax = Axis(fig[(k + 1) ÷ 2, 2 - k % 2]; aspect = DataAspect(), title = label)
    hidedecorations!(ax); hidespines!(ax)
    image!(ax, X, Y, bg[:, j₀+1:j₀+bh])
    image!(ax, X, Y, overlay(field, f, coast); interpolate = true)
    draw_emblem_logo!(ax, ED; scheme = Symbol(get(ENV, "LOGO_SCHEME", string(default_scheme(layout)))))
    limits!(ax, X..., Y...)
end
resize_to_layout!(fig)
out = @sprintf("output/banner_styles%s_t%05.2fd.png", which_set == "default" ? "" : "_" * which_set, ζs.times[n] / 86400)
save(out, fig)
@info "Saved $out"
