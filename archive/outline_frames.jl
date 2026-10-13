# Frames for an animated web banner of white vorticity outlines (contours of |ζ| at the given
# levels, no fill) around the text logo, on a transparent background. One PNG per saved output
# time between the given days. Assemble with ffmpeg (VP9-alpha WebM / HEVC-alpha).
#
# Usage: julia --project=. scripts/outline_frames.jl [name] [levels in f, comma-separated] [first day] [last day] [width px] [linewidth]
#   LINE_ALPHA = opacity of the outlines (default 1); LEVEL = surface (default) | bottom
using Oceananigans, CairoMakie, Printf
using CairoMakie.Colors: RGB
include(joinpath(@__DIR__, "..", "src", "plotting.jl"))
include(joinpath(@__DIR__, "..", "src", "logo_config.jl"))
include(joinpath(@__DIR__, "..", "src", "letter_styles.jl"))
include(joinpath(@__DIR__, "..", "src", "emblem_overlay.jl"))

name   = length(ARGS) ≥ 1 ? ARGS[1] : "heated_text_694x231x8"
levels = length(ARGS) ≥ 2 ? parse.(Float64, split(ARGS[2], ",")) : [0.4]
day₁   = length(ARGS) ≥ 3 ? parse(Float64, ARGS[3]) : 0.0
day₂   = length(ARGS) ≥ 4 ? parse(Float64, ARGS[4]) : 16.0
W      = length(ARGS) ≥ 5 ? parse(Int, ARGS[5]) : 1100
lw     = length(ARGS) ≥ 6 ? parse(Float64, ARGS[6]) : 1.5

level_name = get(ENV, "LEVEL", "surface")                  # surface | bottom
ζs = FieldTimeSeries(joinpath("output", "$(name)_$(level_name).jld2"), "ζ"; backend = OnDisk())
frames = findall(t -> day₁ ≤ t / 86400 ≤ day₂ + 1e-6, ζs.times)
grid = ζs.grid; ug = grid.underlying_grid
Nx, Ny = size(grid, 1), size(grid, 2)
xf = xnodes(ug, Face())[1:Nx] ./ 1e3; yf = ynodes(ug, Face())[1:Ny] ./ 1e3
Δx_km = xf[2] - xf[1]
land = Array(grid.immersed_boundary.bottom_height[1:Nx, 1:Ny, 1]) .== 0
lakes = enclosed_water(land; periodic = (true, true))
r = 2                                                      # contour on a 2× upsampled field
ED = emblem_overlay_data(; cells_per_Ld = 10 / Δx_km, Δx_km, x₀_km = xf[1], y₀_km = yf[1], layout = layout_from_name(name), r)

function upsample(A, r)                                    # bilinear, periodic, cell-centered
    nx, ny = size(A); B = zeros(r * nx, r * ny)
    for J in 1:r*ny, I in 1:r*nx
        x = (I - 0.5) / r + 0.5; y = (J - 0.5) / r + 0.5
        i₀ = floor(Int, x); j₀ = floor(Int, y); a = x - i₀; b = y - j₀
        i₁, i₂ = mod1(i₀, nx), mod1(i₀ + 1, nx); j₁, j₂ = mod1(j₀, ny), mod1(j₀ + 1, ny)
        B[I, J] = (1-a)*(1-b)*A[i₁, j₁] + a*(1-b)*A[i₂, j₁] + (1-a)*b*A[i₁, j₂] + a*b*A[i₂, j₂]
    end
    B
end
field(n) = (ζ = interior(ζs[n], 1:Nx, 1:Ny, 1) ./ 1e-4; ζ[land .| lakes] .= 0; upsample(ζ, r))

X = (xf[1], xf[end] + Δx_km); Y = (yf[1], yf[end] + Δx_km)
Hpx = round(Int, W * Ny / Nx)
fig = Figure(size = (W, Hpx), figure_padding = 0, backgroundcolor = :transparent)
ax = Axis(fig[1, 1]; backgroundcolor = :transparent)
hidedecorations!(ax); hidespines!(ax); limits!(ax, X..., Y...)
Z = Observable(field(frames[1]))
# Lines are drawn opaque white; LINE_ALPHA is applied afterwards to the saved frames, away from the
# letters (drawing semi-transparent strokes directly leaves dark anti-aliasing fringes)
line_alpha = parse(Float64, get(ENV, "LINE_ALPHA", "1.0"))
contour!(ax, ED.xs, ED.ys, Z; levels = sort(vcat(-levels, levels)), color = :white, linewidth = lw)
letter_fill, outc, _ = TEXT_SCHEMES[:white]
draw_styled_letters!(ax, ED.xs, ED.ys, ED.letters; style = :flat, stroke_px = ED.stroke_px,
                     fills = (letter_fill, letter_fill), outline_color = outc, outline_width = 1.2)

tag = "outline_" * (level_name == "surface" ? "" : level_name * "_") * join(replace.(string.(levels), "." => ""), "_")
dir = joinpath("output", "web_frames", tag); mkpath(dir); foreach(f -> rm(joinpath(dir, f)), readdir(dir))
# Letter region on the canvas (letters plus a 2-pixel margin for their outline), for the alpha step
near_letters = distance_to(ED.letters) .≤ 2r
function letter_pixel(row, col, H, Wc)
    nfx, nfy = size(near_letters)
    near_letters[clamp(ceil(Int, col / Wc * nfx), 1, nfx), clamp(ceil(Int, (H - row + 1) / H * nfy), 1, nfy)]
end
for (k, n) in enumerate(frames)
    Z[] = field(n)
    path = joinpath(dir, @sprintf("frame_%03d.png", k))
    save(path, fig; px_per_unit = 1)
    if line_alpha < 1
        img = CairoMakie.FileIO.load(path); H, Wc = size(img)
        for row in 1:H, col in 1:Wc
            c = img[row, col]
            (c.alpha > 0 && !letter_pixel(row, col, H, Wc)) &&
                (img[row, col] = typeof(c)(c.r, c.g, c.b, Float64(c.alpha) * line_alpha))
        end
        CairoMakie.FileIO.save(path, img)
    end
end
@info @sprintf("Saved %d frames (%d × %d), days %.2f–%.2f, to %s", length(frames), W, Hpx,
               ζs.times[frames[1]] / 86400, ζs.times[frames[end]] / 86400, dir)
