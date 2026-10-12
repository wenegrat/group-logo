# Athletic-script logo options (Terps-like treatment of free script fonts), plotting only.
# Each option: slant + rising baseline + heavier strokes; red fill, black outline and offset shadow.
# Usage: julia --project=. scripts/script_options.jl
using CairoMakie, Printf, Statistics
using CairoMakie.Colors: RGB
include(joinpath(@__DIR__, "..", "src", "Logo.jl"))
include(joinpath(@__DIR__, "..", "src", "letter_styles.jl"))

yellowtail = FTFont(fontfile("yellowtail/Yellowtail-Regular.ttf"))
kaushan    = FTFont(fontfile("kaushanscript/KaushanScript-Regular.ttf"))
damion     = FTFont(fontfile("damion/Damion-Regular.ttf"))

function script_lines(font; umd_scale = 0.8, swash = false, one_line = false, tracking = 0.0)
    one_line && return Any[Line(text = "UMD Ocean Dynamics", font = font, size = 200.0, tracking = tracking)]
    lines = Any[two_line(font; line2 = "Ocean Dynamics", umd_scale, tracking, gap = 0.35)...]
    swash && push!(lines, Swash(lines[2], 0.25, 0.35))
    return lines
end
terps_post(B, s) = thicken(slant_rotate(B; shear = 0.12, angle = deg2rad(7)), 0.15s)

options = [("Yellowtail",               script_lines(yellowtail)),
           ("Kaushan Script",           script_lines(kaushan)),
           ("Damion",                   script_lines(damion)),
           ("Yellowtail + swash",       script_lines(yellowtail; swash = true)),
           ("Kaushan Script + swash",   script_lines(kaushan; swash = true)),
           ("Yellowtail, one line",     script_lines(yellowtail; one_line = true))]

md_red, ink = RGB(0.886, 0.094, 0.20), RGB(0.05, 0.05, 0.05)
r, cps, margins = 4, 10, (40, 25)
fig = Figure(size = (1800, 1300))
for (k, (label, lines)) in enumerate(options)
    land = logo_mask(lines; cells_per_stroke = cps, margins, post = terps_post)       # simulation grid
    fine = logo_mask(lines; cells_per_stroke = r * cps, grid = r .* size(land), post = terps_post)
    # Thin-stroke check: 10th percentile of horizontal land run lengths, in grid cells
    runs = Int[]; for j in axes(land, 2); n = 0
        for i in axes(land, 1); if land[i, j]; n += 1 elseif n > 0; push!(runs, n); n = 0 end; end; end
    thin = sort(runs)[max(1, round(Int, 0.1length(runs)))]
    stats = @sprintf("grid %d × %d at 10 cells/stroke; thinnest 10%% of strokes ≤ %d cells", size(land)..., thin)
    @info "$label: $stats"
    ax = Axis(fig[(k + 1) ÷ 2, 2 - k % 2]; aspect = DataAspect(), backgroundcolor = RGBf(0.05, 0.12, 0.25),
              title = "$label\n$stats")
    hidedecorations!(ax); hidespines!(ax)
    Nf = size(fine)
    draw_styled_letters!(ax, collect(1:Nf[1]) ./ r, collect(1:Nf[2]) ./ r, fine; style = :retro_flat,
                         stroke_px = r * cps, fills = (md_red, md_red), outline_color = ink)
    limits!(ax, 0, Nf[1] / r, 0, Nf[2] / r)
end
resize_to_layout!(fig)
save("output/script_options.png", fig)
@info "Saved output/script_options.png"
