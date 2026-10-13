# A static banner (no animation): one frame of signed surface vorticity, smoothly upsampled, colored
# and optionally shaded, composited over a background photo with the text logo. Writes several
# variants to output/banner_still/ and a comparison sheet.
#
# Usage: julia --project=. scripts/banner_still.jl [name] [day] [width px] [background image]
#   SHOW_LAKES=1  color the water sealed inside letters instead of leaving it clear
#   NO_BG=1       no background photo; transparent PNGs (plus a version cropped to the content)
#   FADE=0        no edge vignette
#   ONLY=<key>    render one variant only (e.g. 3_shaded_coast); skips the comparison sheet
#   POS_COLOR, NEG_COLOR, PALE_COLOR (hex), TAG: colors for ζ > 0, ζ < 0, weak |ζ|; file-name tag
#   ONLY=custom with SHADE, Z0, Z1, AMAX: embossing strength, transparent range of |ζ|/f, max opacity
#   SMOOTH_Z=σ, FEATHER=σ (grid cells): smooth the vorticity / soften shape edges
#   HALO=d0,w (grid cells): fade out between d0 and d0+w from the letters
#   CONTINUOUS=1  color over the full |ζ| range with a pale → color → deeper ramp (less saturated)
using Oceananigans, CairoMakie, Printf, Statistics
using CairoMakie.Colors: RGBA, RGB
include(joinpath(@__DIR__, "..", "src", "plotting.jl"))
include(joinpath(@__DIR__, "..", "src", "logo_config.jl"))
include(joinpath(@__DIR__, "..", "src", "letter_styles.jl"))
include(joinpath(@__DIR__, "..", "src", "emblem_overlay.jl"))
include(joinpath(@__DIR__, "..", "src", "TextMask.jl"))
using .TextMask: gaussian_smooth

name = length(ARGS) ≥ 1 ? ARGS[1] : "heated_text_694x231x8"
day  = length(ARGS) ≥ 2 ? parse(Float64, ARGS[2]) : 10.0
W    = length(ARGS) ≥ 3 ? parse(Int, ARGS[3]) : 2200
bg_path = length(ARGS) ≥ 4 ? ARGS[4] : joinpath(REPO, "output", "website_mock", "site", "img", "Landsat.jpg")

ζs = FieldTimeSeries(joinpath("output", "$(name)_surface.jld2"), "ζ"; backend = OnDisk())
n = argmin(abs.(ζs.times ./ 86400 .- day))
grid = ζs.grid; ug = grid.underlying_grid
Nx, Ny = size(grid, 1), size(grid, 2)
xf = xnodes(ug, Face())[1:Nx] ./ 1e3; yf = ynodes(ug, Face())[1:Ny] ./ 1e3
Δx_km = xf[2] - xf[1]
land = Array(grid.immersed_boundary.bottom_height[1:Nx, 1:Ny, 1]) .== 0
lakes = enclosed_water(land; periodic = (true, true))
layout = layout_from_name(name)
r = 4
ED = emblem_overlay_data(; cells_per_Ld = 10 / Δx_km, Δx_km, x₀_km = xf[1], y₀_km = yf[1], layout, r)

# ── Fields on the fine (r×) grid ─────────────────────────────────────────────
ζ = interior(ζs[n], 1:Nx, 1:Ny, 1) ./ 1e-4
ζ[land] .= 0                                         # land is masked later; zero avoids spikes at coasts

# Bilinear upsampling by r (periodic), cell-centered
function upsample(A, r)
    nx, ny = size(A)
    B = zeros(r * nx, r * ny)
    for J in 1:r*ny, I in 1:r*nx
        x = (I - 0.5) / r + 0.5; y = (J - 0.5) / r + 0.5
        i₀ = floor(Int, x); j₀ = floor(Int, y); a = x - i₀; b = y - j₀
        i₁, i₂ = mod1(i₀, nx), mod1(i₀ + 1, nx); j₁, j₂ = mod1(j₀, ny), mod1(j₀ + 1, ny)
        B[I, J] = (1-a)*(1-b)*A[i₁, j₁] + a*(1-b)*A[i₂, j₁] + (1-a)*b*A[i₁, j₂] + a*b*A[i₂, j₂]
    end
    return B
end
Z = upsample(ζ, r)
# SMOOTH_Z=σ (grid cells): Gaussian-smooth the vorticity before coloring, removing small filaments
smooth_z = parse(Float64, get(ENV, "SMOOTH_Z", "0"))
smooth_z > 0 && (Z = gaussian_smooth(Z, smooth_z * r))
fine_land = ED.letters .| ED.shadow                  # land on the fine grid (text layout: letters only)
fine_lakes = upsample(Float64.(lakes), r) .> 0.5
Nfx, Nfy = size(Z)

# SHOW_LAKES=1: color the water sealed inside letters (O, D, A counters) instead of leaving it clear
show_lakes = get(ENV, "SHOW_LAKES", "0") == "1"
smooth(s) = (s = clamp(s, 0, 1); s^2 * (3 - 2s))
rampf(k, N, w) = smooth(min(k - 0.5, N - k + 0.5) / (w * N))
# FADE=0: no edge vignette (the colored patch simply ends where the vorticity drops to near zero)
fade = get(ENV, "FADE", "1") != "0"
edge = fade ? [rampf(i, Nfx, 0.12) * rampf(j, Nfy, 0.30) for i in 1:Nfx, j in 1:Nfy] : ones(Nfx, Nfy)
coastfade = smooth.((distance_to(fine_land) ./ r .- 1) ./ 4)
show_lakes && (coastfade[fine_lakes] .= 1.0)          # counters are narrow; don't fade them away

# Hillshade of |ζ| (smoothed) as a height field: light from the upper left
H = gaussian_smooth(abs.(Z), 1.5r)
function hillshade(H; az = deg2rad(135), alt = deg2rad(40), exag = 6.0)
    nx, ny = size(H); S = similar(H)
    for j in 1:ny, i in 1:nx
        dx = (H[min(i+1, nx), j] - H[max(i-1, 1), j]) / 2 * exag
        dy = (H[i, min(j+1, ny)] - H[i, max(j-1, 1)]) / 2 * exag
        nrm = sqrt(dx^2 + dy^2 + 1)
        S[i, j] = (sin(alt) + cos(alt) * (-dx * cos(az) - dy * sin(az))) / nrm
    end
    return S ./ sin(alt)                               # 1 = flat
end
S = hillshade(H)

lerp(a, b, t) = RGB(a.r + (b.r - a.r) * t, a.g + (b.g - a.g) * t, a.b + (b.b - a.b) * t)
# Colors (hex): POS_COLOR for cyclones (ζ > 0), NEG_COLOR for anticyclones, PALE_COLOR at weak |ζ|
hexcolor(h) = (h = lstrip(h, '#'); RGB(parse(Int, h[1:2], base = 16) / 255, parse(Int, h[3:4], base = 16) / 255, parse(Int, h[5:6], base = 16) / 255))
bright_aqua  = hexcolor(get(ENV, "POS_COLOR", "26E6EB"))
bright_coral = hexcolor(get(ENV, "NEG_COLOR", "FF8059"))
pale         = hexcolor(get(ENV, "PALE_COLOR", "EBF7F5"))
tag = get(ENV, "TAG", "")                                   # appended to output file names

# CONTINUOUS=1: spread the color over the actual range of |ζ| (full color at its 98th percentile in
# open water) with a three-stop ramp per sign, pale → color → deeper shade, instead of saturating at 1 f
continuous = get(ENV, "CONTINUOUS", "0") == "1"
zq = let w = abs.(Z[.!fine_land .& .!fine_lakes]); w = sort(w[w .> 0.05]); w[ceil(Int, 0.98length(w))] end
continuous && @info @sprintf("98th percentile of |ζ|/f in open water: %.2f", zq)
deeper(c) = RGB(0.55c.r, 0.55c.g, 0.55c.b)
function ramp_color(z, zfull)
    base = z > 0 ? bright_aqua : bright_coral
    t = clamp(abs(z) / zfull, 0, 1)
    continuous || return lerp(pale, base, t)
    t^0.8 < 0.6 ? lerp(pale, base, t^0.8 / 0.6) : lerp(base, deeper(base), (t^0.8 - 0.6) / 0.4)
end

function colorize(; shade = 0.0, coast = true, z₀ = 0.05, z₁ = 0.3, zfull = continuous ? zq : 1.0, amax = 1.0)
    O = Matrix{RGBA{Float64}}(undef, Nfx, Nfy)
    for j in 1:Nfy, i in 1:Nfx
        z = Z[i, j]
        c = ramp_color(z, zfull)
        a = amax * smooth((abs(z) - z₀) / (z₁ - z₀)) * edge[i, j] * (coast ? coastfade[i, j] : 1.0)
        (fine_land[i, j] || (fine_lakes[i, j] && !show_lakes)) && (a = 0.0)
        k = 1 + shade * (S[i, j] - 1)
        O[i, j] = RGBA(clamp(c.r * k, 0, 1), clamp(c.g * k, 0, 1), clamp(c.b * k, 0, 1), a)
    end
    # FEATHER=σ (grid cells): blur the opacity so shapes have soft edges
    feather = parse(Float64, get(ENV, "FEATHER", "0"))
    if feather > 0
        A = gaussian_smooth([o.alpha for o in O], feather * r)
        landish = fine_land .| (fine_lakes .& !show_lakes)
        O = [landish[i, j] ? RGBA(O[i, j].r, O[i, j].g, O[i, j].b, 0.0) : RGBA(O[i, j].r, O[i, j].g, O[i, j].b, A[i, j])
             for i in 1:Nfx, j in 1:Nfy]
    end
    # HALO=d0,w (grid cells): fade out between d0 and d0 + w from the letters
    halo = get(ENV, "HALO", "")
    if !isempty(halo)
        d₀, w = parse.(Float64, split(halo, ","))
        dl = distance_to(ED.letters) ./ r
        O = [RGBA(O[i, j].r, O[i, j].g, O[i, j].b, O[i, j].alpha * (1 - smooth((dl[i, j] - d₀) / w))) for i in 1:Nfx, j in 1:Nfy]
    end
    return O
end

# Soft dark glow around the letters (decorative, for legibility)
dout = distance_to(ED.letters) ./ r                       # in grid cells (2 km)
glow(strength) = [RGBA(0.02, 0.06, 0.12, strength * exp(-dout[i, j] / 2.5) * !ED.letters[i, j]) for i in 1:Nfx, j in 1:Nfy]

variant_specs = [
    ("1_smooth",          "1. smooth upsampled colors (style F)",                () -> colorize(), () -> nothing),
    ("2_shaded",          "2. + embossed shading",                               () -> colorize(shade = 0.6), () -> nothing),
    ("3_shaded_coast",    "3. + vorticity kept along the coastline",            () -> colorize(shade = 0.6, coast = false), () -> nothing),
    # SHADE (embossing strength, default 0.6), Z0/Z1 (transparent below |ζ| = Z0 f, opaque above Z1 f)
    ("custom",            "custom shade / transparency",
     () -> colorize(shade = parse(Float64, get(ENV, "SHADE", "0.6")), coast = false,
                    z₀ = parse(Float64, get(ENV, "Z0", "0.05")), z₁ = parse(Float64, get(ENV, "Z1", "0.3")),
                    amax = parse(Float64, get(ENV, "AMAX", "1.0"))), () -> nothing),
    ("4_shaded_glow",     "4. shaded, coastline kept, soft dark glow behind text", () -> colorize(shade = 0.6, coast = false), () -> glow(0.55)),
]
only₀ = get(ENV, "ONLY", "")
isempty(only₀) || filter!(v -> v[1] == only₀, variant_specs)
variants = [(k, l, o(), g()) for (k, l, o, g) in variant_specs]

bg = CairoMakie.FileIO.load(bg_path); bg = permutedims(reverse(bg, dims = 1))
bw = size(bg, 1); bh = min(size(bg, 2), round(Int, bw * Ny / Nx)); j₀ = (size(bg, 2) - bh) ÷ 2
bgcrop = bg[:, j₀+1:j₀+bh]
X = (xf[1], xf[end] + Δx_km); Y = (yf[1], yf[end] + Δx_km)
Hpx = round(Int, W * Ny / Nx)

# NO_BG=1: leave out the photo and save with a transparent background (to overlay on a web page)
no_bg = get(ENV, "NO_BG", "0") == "1"

function draw!(ax, overlay, gl)
    no_bg || image!(ax, X, Y, bgcrop)
    image!(ax, X, Y, overlay; interpolate = true)
    isnothing(gl) || image!(ax, X, Y, gl; interpolate = true)
    fill, outc, _ = TEXT_SCHEMES[:white]
    draw_styled_letters!(ax, ED.xs, ED.ys, ED.letters; style = :flat, stroke_px = ED.stroke_px,
                         fills = (fill, fill), outline_color = outc, outline_width = 2.0)
    hidedecorations!(ax); hidespines!(ax); limits!(ax, X..., Y...)
end

dir = joinpath("output", "banner_still"); mkpath(dir)
for (key, label, ov, gl) in variants
    fig = Figure(size = (W, Hpx), figure_padding = 0, backgroundcolor = no_bg ? :transparent : :white)
    draw!(Axis(fig[1, 1]; backgroundcolor = :transparent), ov, gl)
    path = joinpath(dir, @sprintf("banner_t%05.2fd_%s%s%s%s%s.png", ζs.times[n] / 86400, key, show_lakes ? "_lakes" : "",
                                  fade ? "" : "_nofade", no_bg ? "_transparent" : "", isempty(tag) ? "" : "_" * tag))
    save(path, fig; px_per_unit = 1)
    if no_bg   # crop transparent margins (small pad) so the banner sits tightly in the page layout
        img = CairoMakie.FileIO.load(path)
        A = [Float64(c.alpha) > 0.02 for c in img]
        rows = findall(vec(any(A, dims = 2))); cols = findall(vec(any(A, dims = 1)))
        pad = 8
        # symmetric about the image center, where the text is centered, so the words stay centered on the page
        hr = max(size(img, 1) ÷ 2 - first(rows), last(rows) - size(img, 1) ÷ 2) + pad
        hc = max(size(img, 2) ÷ 2 - first(cols), last(cols) - size(img, 2) ÷ 2) + pad
        r₁, r₂ = max(size(img, 1) ÷ 2 - hr, 1), min(size(img, 1) ÷ 2 + hr, size(img, 1))
        c₁, c₂ = max(size(img, 2) ÷ 2 - hc, 1), min(size(img, 2) ÷ 2 + hc, size(img, 2))
        CairoMakie.FileIO.save(replace(path, ".png" => "_cropped.png"), img[r₁:r₂, c₁:c₂])
    end
end

isempty(only₀) || (@info "Saved $(only₀)"; exit())
sheet = Figure(size = (1800, 2 * round(Int, 880 * Ny / Nx) + 140))
for (k, (key, label, ov, gl)) in enumerate(variants)
    ax = Axis(sheet[(k + 1) ÷ 2, 2 - k % 2]; aspect = DataAspect(), title = label)
    draw!(ax, ov, gl)
end
resize_to_layout!(sheet)
out = joinpath(dir, @sprintf("comparison_t%05.2fd%s.png", ζs.times[n] / 86400, show_lakes ? "_lakes" : ""))
save(out, sheet)
@info "Saved $(length(variants)) banners ($W × $Hpx) and $out"
