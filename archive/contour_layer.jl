# Transparent banner layer of white vorticity contour lines (no fill), on the same canvas as
# scripts/banner_still.jl (so it can be composited with its layers by scripts/layer_banners.jl).
#
# Usage: julia --project=. scripts/contour_layer.jl [name] [day] [levels, comma-separated, in f] [linewidth px] [tag]
using Oceananigans, CairoMakie, Printf
include(joinpath(@__DIR__, "..", "src", "plotting.jl"))
include(joinpath(@__DIR__, "..", "src", "logo_config.jl"))
include(joinpath(@__DIR__, "..", "src", "letter_styles.jl"))
include(joinpath(@__DIR__, "..", "src", "emblem_overlay.jl"))

name   = length(ARGS) ≥ 1 ? ARGS[1] : "heated_text_694x231x8"
day    = length(ARGS) ≥ 2 ? parse(Float64, ARGS[2]) : 17.875
levels = length(ARGS) ≥ 3 ? parse.(Float64, split(ARGS[3], ",")) : [0.6]
lw     = length(ARGS) ≥ 4 ? parse(Float64, ARGS[4]) : 2.0
tag    = length(ARGS) ≥ 5 ? ARGS[5] : "contours"
W = 2200

ζs = FieldTimeSeries(joinpath("output", "$(name)_surface.jld2"), "ζ"; backend = OnDisk())
n = argmin(abs.(ζs.times ./ 86400 .- day))
grid = ζs.grid; ug = grid.underlying_grid
Nx, Ny = size(grid, 1), size(grid, 2)
xf = xnodes(ug, Face())[1:Nx] ./ 1e3; yf = ynodes(ug, Face())[1:Ny] ./ 1e3
Δx_km = xf[2] - xf[1]
land = Array(grid.immersed_boundary.bottom_height[1:Nx, 1:Ny, 1]) .== 0
lakes = enclosed_water(land; periodic = (true, true))
r = 4
ED = emblem_overlay_data(; cells_per_Ld = 10 / Δx_km, Δx_km, x₀_km = xf[1], y₀_km = yf[1], layout = layout_from_name(name), r)

ζ = interior(ζs[n], 1:Nx, 1:Ny, 1) ./ 1e-4
ζ[land .| lakes] .= 0
function upsample(A, r)                                   # bilinear, periodic, cell-centered
    nx, ny = size(A); B = zeros(r * nx, r * ny)
    for J in 1:r*ny, I in 1:r*nx
        x = (I - 0.5) / r + 0.5; y = (J - 0.5) / r + 0.5
        i₀ = floor(Int, x); j₀ = floor(Int, y); a = x - i₀; b = y - j₀
        i₁, i₂ = mod1(i₀, nx), mod1(i₀ + 1, nx); j₁, j₂ = mod1(j₀, ny), mod1(j₀ + 1, ny)
        B[I, J] = (1-a)*(1-b)*A[i₁, j₁] + a*(1-b)*A[i₂, j₁] + (1-a)*b*A[i₁, j₂] + a*b*A[i₂, j₂]
    end
    B
end
Z = upsample(ζ, r)
mask = ED.letters .| (upsample(Float64.(lakes), r) .> 0.5)
Z[mask] .= 0                                             # Makie's contour doesn't accept NaN; erased below

X = (xf[1], xf[end] + Δx_km); Y = (yf[1], yf[end] + Δx_km)
Hpx = round(Int, W * Ny / Nx)
fig = Figure(size = (W, Hpx), figure_padding = 0, backgroundcolor = :transparent)
ax = Axis(fig[1, 1]; backgroundcolor = :transparent)
hidedecorations!(ax); hidespines!(ax); limits!(ax, X..., Y...)
contour!(ax, ED.xs, ED.ys, Z; levels = sort(vcat(-levels, levels)), color = :white, linewidth = lw)
dir = joinpath("output", "banner_still")
out = joinpath(dir, @sprintf("banner_t%05.2fd_%s.png", ζs.times[n] / 86400, tag))
save(out, fig; px_per_unit = 1)
# Erase lines drawn over land and sealed lakes (nearest-neighbour map from canvas pixels to the fine grid)
img = CairoMakie.FileIO.load(out)                        # rows top → bottom
H, Wc = size(img); nfx, nfy = size(mask)
for row in 1:H, col in 1:Wc
    i = clamp(ceil(Int, col / Wc * nfx), 1, nfx); j = clamp(ceil(Int, (H - row + 1) / H * nfy), 1, nfy)
    mask[i, j] && (img[row, col] = typeof(img[row, col])(0, 0, 0, 0))
end
CairoMakie.FileIO.save(out, img)
@info "Saved $out"
