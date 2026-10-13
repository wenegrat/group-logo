# Preview the logo layout without running anything: the simulation land mask and the styled
# rendering, which draws the drop shadow exactly on its land.
#
# Usage: julia --project=. scripts/preview_logo.jl [cells per Ld (default 7.5)] [layout: text (default) | emblem]
#   LOGO_SCHEME: see src/emblem_overlay.jl
using CairoMakie, Printf
include(joinpath(@__DIR__, "..", "src", "logo_config.jl"))
include(joinpath(@__DIR__, "..", "src", "letter_styles.jl"))
include(joinpath(@__DIR__, "..", "src", "emblem_overlay.jl"))

cells_per_Ld = isempty(ARGS) ? 7.5 : parse(Float64, ARGS[1])
layout = length(ARGS) ≥ 2 ? Symbol(ARGS[2]) : :text
Δx_km = 10 / cells_per_Ld                                  # Ld = 10 km
m = logo_layout(cells_per_Ld; layout)
@info @sprintf("grid %d × %d (Δx = %.1f km), land %.1f%% of the domain", m.grid..., Δx_km, 100sum(m.land) / length(m.land))

ED = emblem_overlay_data(; cells_per_Ld, Δx_km, x₀_km = 0.0, y₀_km = 0.0, layout)
Lx, Ly = m.grid .* Δx_km
navy = RGBf(0.05, 0.12, 0.25)
fig = Figure(size = (1800, round(Int, 2 * 900 * Ly / Lx) + 120))
ax1 = Axis(fig[1, 1]; aspect = DataAspect(), backgroundcolor = navy, title = "simulation land mask")
heatmap!(ax1, ((1:m.grid[1]) .- 0.5) .* Δx_km, ((1:m.grid[2]) .- 0.5) .* Δx_km, ifelse.(m.land, 1.0, NaN);
         colormap = [:white], colorrange = (0, 1))
ax2 = Axis(fig[1, 2]; aspect = DataAspect(), backgroundcolor = navy, title = "rendered logo")
draw_emblem_logo!(ax2, ED; scheme = Symbol(get(ENV, "LOGO_SCHEME", string(default_scheme(layout)))))
for ax in (ax1, ax2); hidedecorations!(ax); limits!(ax, 0, Lx, 0, Ly); end
mkpath("output")
save("output/logo_preview.png", fig)
@info "Saved output/logo_preview.png"
