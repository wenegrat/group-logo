# Terrapin emblem + "Ocean Dynamics" on one line (plotting only). The emblem is scaled so its
# artwork (not the PNG canvas) is as tall as the text's ink.
#   1. Bubble Vibes, retro flat, pastel gold / black
#   2–3. Kaushan Script and Damion with a Terps-like treatment (heavy, tight, italic, horizontal);
#        the emblem matches the capital height and stands on the baseline
# Usage: julia --project=. scripts/emblem_options.jl [path/to/emblem.png]
using CairoMakie, Printf
using CairoMakie.Colors: RGB, RGBA
include(joinpath(@__DIR__, "..", "src", "Logo.jl"))
include(joinpath(@__DIR__, "..", "src", "letter_styles.jl"))

path = isempty(ARGS) ? expanduser("~/Downloads/maryland_terrapins_logo_primary_19737774.png") : ARGS[1]
emblem, img = load_emblem(path)

vibes   = FTFont(fontfile("bubblevibes/BubbleVibes-Regular.otf"))
kaushan = FTFont(fontfile("kaushanscript/KaushanScript-Regular.ttf"))

# Bubble text
Tb = layout_bitmap([Line(text = "OCEAN DYNAMICS", font = vibes, size = 200.0, tracking = 0.02)])
damion = FTFont(fontfile("damion/Damion-Regular.ttf"))

# Terps-like script on a horizontal baseline: tight spacing, italic shear, heavy strokes
function terps_script(font; tracking = -0.03, shear = 0.18, weight = 0.32)
    line(t) = Line(text = t, font = font, size = 200.0, tracking = tracking)
    s₀ = median_stroke(layout_bitmap([line("Ocean Dynamics")]))
    post(B) = thicken(slant_rotate(B; shear, angle = 0.0), weight * s₀)
    T = post(layout_bitmap([line("Ocean Dynamics")]))
    # Cap height and descender depth of the processed text, from single glyphs
    cap = size(post(layout_bitmap([line("D")])), 2)
    desc = size(post(layout_bitmap([line("Dy")])), 2) - cap
    return T, cap, desc
end
Tk, capk, desck = terps_script(kaushan)
Td, capd, descd = terps_script(damion)

cps, r = 10, 4
margins_Ld = (9.4, 7.7)
gold, ink = RGB(0.98, 0.84, 0.45), RGB(0.08, 0.08, 0.08)
md_red = RGB(0.886, 0.094, 0.20)

logo_red, white = RGB(0.927, 0.106, 0.18), RGB(1.0, 1.0, 1.0)     # red sampled from the emblem artwork
# (label, fills (left word, right word), outline, shadow, outline width, split between words?)
options = [("red / black (like the M)",                         (logo_red, logo_red), ink, ink, 2.0, false),
           ("white / black (like the turtle)",                  (white, white), ink, ink, 2.0, false),
           ("red with white rim, black shadow",                 (logo_red, logo_red), white, ink, 3.0, false),
           ("white OCEAN, red DYNAMICS, black",                 (white, logo_red), ink, ink, 2.0, true)]

fig = Figure(size = (1800, 900))
T, cap, desc = Tb, nothing, 0
for (k, (label, fills, outc, shc, ow, split)) in enumerate(options)
    stroke = median_stroke(T)
    C, Ct, (i₀, j₀, w, h) = emblem_row(emblem, T; height_ratio = 1.0, gap = 0.12)
    scale = stroke / cps
    nx, ny = floor.(Int, size(C) ./ scale)
    grid = (nx, ny) .+ 2 .* round.(Int, margins_Ld .* 5)
    land, _ = place_bitmap(C, scale, grid)
    stats = @sprintf("grid %d × %d (5 cells per Ld), land %.1f%%", grid..., 100sum(land) / length(land))
    @info "$label: $stats"

    fine = r .* grid
    Lf, (ox, oy) = place_bitmap(C, scale / r, fine)
    Tf, _ = place_bitmap(Ct, scale / r, fine)
    ax = Axis(fig[(k + 1) ÷ 2, 2 - k % 2]; aspect = DataAspect(), backgroundcolor = RGBf(0.05, 0.12, 0.25), title = "$label\n$stats")
    hidedecorations!(ax); hidespines!(ax)
    xs = collect(1:fine[1]) ./ r; ys = collect(1:fine[2]) ./ r
    # Split between words at the middle of the widest gap in the text columns
    cols = [any(Tf[i, :]) for i in 1:fine[1]]
    on = findall(cols); gaps = [(on[m+1] - on[m], m) for m in 1:length(on)-1]
    gw, gm = maximum(gaps); sx = (on[gm] + on[gm+1]) ÷ 2
    draw_styled_letters!(ax, xs, ys, Tf; style = :retro_flat, stroke_px = r * cps, fills,
                         outline_color = outc, shadow_color = shc, outline_width = ow,
                         split_x = split ? sx : nothing)
    x0 = (ox + i₀ * r / scale) / r; y0 = (oy + j₀ * r / scale) / r
    image!(ax, (x0, x0 + w / scale), (y0, y0 + h / scale), img)
    limits!(ax, 0, grid[1], 0, grid[2])
end
resize_to_layout!(fig)
save("output/emblem_colors.png", fig)
@info "Saved output/emblem_colors.png"
