# Export the logo lettering alone as a transparent PNG, aligned with a run's domain, for layering
# over a trails-only video (same framing as scripts/export_velocity.jl).
#
# Usage: julia --project=. scripts/export_letters.jl [name] [width px] [out.png]
using Oceananigans, CairoMakie
include(joinpath(@__DIR__, "..", "src", "plotting.jl"))
include(joinpath(@__DIR__, "..", "src", "logo_config.jl"))
include(joinpath(@__DIR__, "..", "src", "letter_styles.jl"))
include(joinpath(@__DIR__, "..", "src", "emblem_overlay.jl"))

name = length(ARGS) ≥ 1 ? ARGS[1] : "heated_text_694x231x8"
W    = length(ARGS) ≥ 2 ? parse(Int, ARGS[2]) : 2200
out  = length(ARGS) ≥ 3 ? ARGS[3] : joinpath("output", "letters_$(W).png")

us = FieldTimeSeries(joinpath("output", "$(name)_surface.jld2"), "u"; backend = OnDisk())
ug = us.grid.underlying_grid
Nx, Ny = size(ug, 1), size(ug, 2)
xf = xnodes(ug, Face())[1:Nx] ./ 1e3; yf = ynodes(ug, Face())[1:Ny] ./ 1e3
Δx_km = xf[2] - xf[1]
r = max(2, ceil(Int, 2W / Nx))                              # fine mask at ≥ 2 mask pixels per output pixel
ED = emblem_overlay_data(; cells_per_Ld = 10 / Δx_km, Δx_km, x₀_km = xf[1], y₀_km = yf[1], layout = layout_from_name(name), r)
X = (xf[1], xf[end] + Δx_km); Y = (yf[1], yf[end] + Δx_km)
H = round(Int, W * Ny / Nx)
fig = Figure(size = (W, H), figure_padding = 0, backgroundcolor = :transparent)
ax = Axis(fig[1, 1]; backgroundcolor = :transparent)
hidedecorations!(ax); hidespines!(ax); limits!(ax, X..., Y...)
letter_fill, outc, _ = TEXT_SCHEMES[:white]
draw_styled_letters!(ax, ED.xs, ED.ys, ED.letters; style = :flat, stroke_px = ED.stroke_px,
                     fills = (letter_fill, letter_fill), outline_color = outc, outline_width = 1.2 * W / 1100)
save(out, fig; px_per_unit = 1)
@info "Saved $out ($W × $H)"
