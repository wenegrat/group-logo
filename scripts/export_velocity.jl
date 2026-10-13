# Export surface (or bottom) velocities at cell centers for particle-trail rendering
# (scripts/particle_trails.py), with the land mask and the styled letters overlay.
#   <dir>/u.f32, v.f32   Float32 [m/s], frame-major, each frame Nx × Ny (x fastest)
#   <dir>/land.u8        UInt8 land mask, Nx × Ny
#   <dir>/zeta.f32, b.f32  vorticity / f and buoyancy anomaly (10⁻³ m s⁻²) at centers, same layout as u
#   <dir>/meta.txt       Nx Ny nframes Δx_km X0 X1 Y0 Y1 W H, then a line of times (days)
#   <dir>/letters.png    letters on a transparent W × H canvas
#
# Usage: julia --project=. scripts/export_velocity.jl [name] [surface|bottom] [first day] [last day] [width px]
using Oceananigans, CairoMakie, Printf
include(joinpath(@__DIR__, "..", "src", "plotting.jl"))
include(joinpath(@__DIR__, "..", "src", "logo_config.jl"))
include(joinpath(@__DIR__, "..", "src", "letter_styles.jl"))
include(joinpath(@__DIR__, "..", "src", "emblem_overlay.jl"))

name  = length(ARGS) ≥ 1 ? ARGS[1] : "heated_text_694x231x8"
level = length(ARGS) ≥ 2 ? ARGS[2] : "surface"
day₁  = length(ARGS) ≥ 3 ? parse(Float64, ARGS[3]) : 0.0
day₂  = length(ARGS) ≥ 4 ? parse(Float64, ARGS[4]) : 17.5
W     = length(ARGS) ≥ 5 ? parse(Int, ARGS[5]) : 1100

file = joinpath("output", "$(name)_$(level).jld2")
us = FieldTimeSeries(file, "u"; backend = OnDisk()); vs = FieldTimeSeries(file, "v"; backend = OnDisk())
frames = findall(t -> day₁ ≤ t / 86400 ≤ day₂ + 1e-6, us.times)
grid = us.grid; ug = grid.underlying_grid
Nx, Ny = size(grid, 1), size(grid, 2)
xf = xnodes(ug, Face())[1:Nx] ./ 1e3; yf = ynodes(ug, Face())[1:Ny] ./ 1e3
Δx_km = xf[2] - xf[1]
land = Array(grid.immersed_boundary.bottom_height[1:Nx, 1:Ny, 1]) .== 0
ED = emblem_overlay_data(; cells_per_Ld = 10 / Δx_km, Δx_km, x₀_km = xf[1], y₀_km = yf[1], layout = layout_from_name(name), r = 2)

dir = joinpath("output", "velocity_data_$(level)"); mkpath(dir)
open(joinpath(dir, "u.f32"), "w") do iu
    open(joinpath(dir, "v.f32"), "w") do iv
        for n in frames
            u = Float64.(interior(us[n], 1:Nx, 1:Ny, 1)); v = Float64.(interior(vs[n], 1:Nx, 1:Ny, 1))
            uc = (u .+ circshift(u, (-1, 0))) ./ 2; vc = (v .+ circshift(v, (0, -1))) ./ 2   # faces → centers
            uc[land] .= 0; vc[land] .= 0
            write(iu, Float32.(uc)); write(iv, Float32.(vc))
        end
    end
end
write(joinpath(dir, "land.u8"), UInt8.(land))
# Coloring fields at cell centers: vorticity / f and buoyancy anomaly (10⁻³ m s⁻²)
ζs = FieldTimeSeries(file, "ζ"; backend = OnDisk()); bs = FieldTimeSeries(file, "b"; backend = OnDisk())
b₀ = Float64.(interior(bs[1], 1:Nx, 1:Ny, 1))
open(joinpath(dir, "zeta.f32"), "w") do iz
    open(joinpath(dir, "b.f32"), "w") do ib
        for n in frames
            ζ = Float64.(interior(ζs[n], 1:Nx, 1:Ny, 1)) ./ 1e-4             # at cell corners
            ζc = (ζ .+ circshift(ζ, (-1, 0)) .+ circshift(ζ, (0, -1)) .+ circshift(ζ, (-1, -1))) ./ 4
            b′ = (Float64.(interior(bs[n], 1:Nx, 1:Ny, 1)) .- b₀) .* 1e3
            ζc[land] .= 0; b′[land] .= 0
            write(iz, Float32.(ζc)); write(ib, Float32.(b′))
        end
    end
end
X = (xf[1], xf[end] + Δx_km); Y = (yf[1], yf[end] + Δx_km)
Hpx = round(Int, W * Ny / Nx)
open(joinpath(dir, "meta.txt"), "w") do io
    println(io, join((Nx, Ny, length(frames), Δx_km, X..., Y..., W, Hpx), " "))
    println(io, join((@sprintf("%.4f", us.times[n] / 86400) for n in frames), " "))
end
fig = Figure(size = (W, Hpx), figure_padding = 0, backgroundcolor = :transparent)
ax = Axis(fig[1, 1]; backgroundcolor = :transparent)
hidedecorations!(ax); hidespines!(ax); limits!(ax, X..., Y...)
letter_fill, outc, _ = TEXT_SCHEMES[:white]
draw_styled_letters!(ax, ED.xs, ED.ys, ED.letters; style = :flat, stroke_px = ED.stroke_px,
                     fills = (letter_fill, letter_fill), outline_color = outc, outline_width = 1.2)
save(joinpath(dir, "letters.png"), fig; px_per_unit = 1)
@info "Exported $(length(frames)) velocity frames ($Nx × $Ny) and letters ($W × $Hpx) to $dir"
