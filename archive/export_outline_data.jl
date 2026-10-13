# Export what an external renderer needs to draw vorticity outlines around the logo:
#   <dir>/zeta.f32       Float32 fields (vorticity / f, or with FIELD=b the buoyancy anomaly in 10⁻³ m s⁻²,
#                        land filled from neighbouring water), frame-major, each frame nx × ny (x fastest),
#                        upsampled 2×, land (and sealed lakes, unless SHOW_LAKES=1) set to 0
#   <dir>/meta.txt       nx ny nframes x0 x1 y0 y1 W H, then one line of times (days)
#   <dir>/letters.png    the styled letters alone on a transparent W × H canvas
#
# Usage: julia --project=. scripts/export_outline_data.jl [name] [surface|bottom] [first day] [last day] [width px]
using Oceananigans, CairoMakie, Printf
include(joinpath(@__DIR__, "..", "src", "plotting.jl"))
include(joinpath(@__DIR__, "..", "src", "logo_config.jl"))
include(joinpath(@__DIR__, "..", "src", "letter_styles.jl"))
include(joinpath(@__DIR__, "..", "src", "emblem_overlay.jl"))

name  = length(ARGS) ≥ 1 ? ARGS[1] : "heated_text_694x231x8"
level = length(ARGS) ≥ 2 ? ARGS[2] : "bottom"
day₁  = length(ARGS) ≥ 3 ? parse(Float64, ARGS[3]) : 0.0
day₂  = length(ARGS) ≥ 4 ? parse(Float64, ARGS[4]) : 16.0
W     = length(ARGS) ≥ 5 ? parse(Int, ARGS[5]) : 1100

# FIELD = zeta (vorticity / f, default) | b (buoyancy anomaly b − b(t = 0), in 10⁻³ m s⁻²)
#       | ow (minus the Okubo–Weiss parameter, −W/f² = (ζ² − s_n² − s_s²)/f²: positive in eddy cores)
fieldname = get(ENV, "FIELD", "zeta")
ζs = FieldTimeSeries(joinpath("output", "$(name)_$(level).jld2"), "ζ"; backend = OnDisk())
bs = fieldname == "b" ? FieldTimeSeries(joinpath("output", "$(name)_$(level).jld2"), "b"; backend = OnDisk()) : nothing
us = fieldname == "ow" ? FieldTimeSeries(joinpath("output", "$(name)_$(level).jld2"), "u"; backend = OnDisk()) : nothing
vs = fieldname == "ow" ? FieldTimeSeries(joinpath("output", "$(name)_$(level).jld2"), "v"; backend = OnDisk()) : nothing
frames = findall(t -> day₁ ≤ t / 86400 ≤ day₂ + 1e-6, ζs.times)
grid = ζs.grid; ug = grid.underlying_grid
Nx, Ny = size(grid, 1), size(grid, 2)
xf = xnodes(ug, Face())[1:Nx] ./ 1e3; yf = ynodes(ug, Face())[1:Ny] ./ 1e3
Δx_km = xf[2] - xf[1]
land = Array(grid.immersed_boundary.bottom_height[1:Nx, 1:Ny, 1]) .== 0
lakes = enclosed_water(land; periodic = (true, true))
show_lakes = get(ENV, "SHOW_LAKES", "0") == "1"           # keep the water sealed inside letters
r = 2
ED = emblem_overlay_data(; cells_per_Ld = 10 / Δx_km, Δx_km, x₀_km = xf[1], y₀_km = yf[1], layout = layout_from_name(name), r)

function upsample(A, r)
    nx, ny = size(A); B = zeros(r * nx, r * ny)
    for J in 1:r*ny, I in 1:r*nx
        x = (I - 0.5) / r + 0.5; y = (J - 0.5) / r + 0.5
        i₀ = floor(Int, x); j₀ = floor(Int, y); a = x - i₀; b = y - j₀
        i₁, i₂ = mod1(i₀, nx), mod1(i₀ + 1, nx); j₁, j₂ = mod1(j₀, ny), mod1(j₀ + 1, ny)
        B[I, J] = (1-a)*(1-b)*A[i₁, j₁] + a*(1-b)*A[i₂, j₁] + (1-a)*b*A[i₁, j₂] + a*b*A[i₂, j₂]
    end
    B
end

dir = joinpath("output", "outline_data_$(level)" * (fieldname == "zeta" ? "" : "_" * fieldname)); mkpath(dir)

# −W/f² at cell centers from C-grid velocities (u at x-faces, v at y-faces), periodic differences
function minus_okubo_weiss(u, v, Δ)
    ∂x(A) = (circshift(A, (-1, 0)) .- A) ./ Δ;  ∂y(A) = (circshift(A, (0, -1)) .- A) ./ Δ
    avgx(A) = (A .+ circshift(A, (-1, 0))) ./ 2;  avgy(A) = (A .+ circshift(A, (0, -1))) ./ 2
    sn = ∂x(u) .- ∂y(v)                                     # centers (exact on the C-grid)
    ζ  = ∂x(v) .- ∂y(u)                                     # corners
    ss = ∂x(v) .+ ∂y(u)                                     # corners
    ζc, ssc = avgx(avgy(ζ)), avgx(avgy(ss))                 # corners → centers
    return (ζc .^ 2 .- sn .^ 2 .- ssc .^ 2) ./ 1e-8         # / f², f = 1e-4
end

# Fill land with values of the nearest water (repeated neighbour averaging), so contours of a field
# that is large next to the coast (buoyancy) don't trace the coastline itself. Letters are drawn on top.
function fill_land(A, land)
    A = copy(A); known = .!land
    nx, ny = size(A)
    while !all(known)
        B = copy(A); K = copy(known)
        for j in 1:ny, i in 1:nx
            known[i, j] && continue
            s = 0.0; c = 0
            for (di, dj) in ((1, 0), (-1, 0), (0, 1), (0, -1))
                i2, j2 = mod1(i + di, nx), mod1(j + dj, ny)
                known[i2, j2] && (s += A[i2, j2]; c += 1)
            end
            c > 0 && (B[i, j] = s / c; K[i, j] = true)
        end
        A, known = B, K
    end
    return A
end
X = (xf[1], xf[end] + Δx_km); Y = (yf[1], yf[end] + Δx_km)
Hpx = round(Int, W * Ny / Nx)
open(joinpath(dir, "zeta.f32"), "w") do io
    for n in frames
        if fieldname == "ow"
            u = Float64.(interior(us[n], 1:Nx, 1:Ny, 1)); v = Float64.(interior(vs[n], 1:Nx, 1:Ny, 1))
            ζ = minus_okubo_weiss(u, v, Δx_km * 1e3)
            ζ = fill_land(ζ, show_lakes ? land : (land .| lakes))
        elseif fieldname == "b"
            ζ = (interior(bs[n], 1:Nx, 1:Ny, 1) .- interior(bs[1], 1:Nx, 1:Ny, 1)) .* 1e3
            ζ = fill_land(ζ, show_lakes ? land : (land .| lakes))
        else
            ζ = interior(ζs[n], 1:Nx, 1:Ny, 1) ./ 1e-4
            ζ[show_lakes ? land : (land .| lakes)] .= 0
        end
        write(io, Float32.(upsample(ζ, r)))
    end
end
open(joinpath(dir, "meta.txt"), "w") do io
    println(io, join((r * Nx, r * Ny, length(frames), ED.xs[1], ED.xs[end], ED.ys[1], ED.ys[end], X..., Y..., W, Hpx), " "))
    println(io, join((@sprintf("%.4f", ζs.times[n] / 86400) for n in frames), " "))
end

fig = Figure(size = (W, Hpx), figure_padding = 0, backgroundcolor = :transparent)
ax = Axis(fig[1, 1]; backgroundcolor = :transparent)
hidedecorations!(ax); hidespines!(ax); limits!(ax, X..., Y...)
letter_fill, outc, _ = TEXT_SCHEMES[:white]
draw_styled_letters!(ax, ED.xs, ED.ys, ED.letters; style = :flat, stroke_px = ED.stroke_px,
                     fills = (letter_fill, letter_fill), outline_color = outc, outline_width = 1.2)
save(joinpath(dir, "letters.png"), fig; px_per_unit = 1)
@info "Exported $(length(frames)) frames ($(r*Nx) × $(r*Ny)) and letters ($W × $Hpx) to $dir"
