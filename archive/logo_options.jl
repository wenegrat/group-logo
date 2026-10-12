# Preview lettering/layout options for the logo, rasterized onto the simulation grid
# (no simulation is run). Each option is scaled so the median stroke of its thinnest line is
# `cells_per_stroke` cells wide, and placed in a domain with the same margins as the 512×128 Eady run.
#
# Usage: julia --project=. scripts/logo_options.jl [first | block | bubble | vibes | vibes2 | inflated] [block height]
using CairoMakie, FreeTypeAbstraction, Statistics, Printf
include(joinpath(@__DIR__, "..", "src", "Logo.jl"))
include(joinpath(@__DIR__, "..", "src", "plotting.jl"))
include(joinpath(@__DIR__, "..", "src", "TextMask.jl"))
using .TextMask: gaussian_smooth

cells_per_stroke = parse(Int, get(ENV, "CELLS_PER_STROKE", "5"))
margins = (80, 45)                      # cells of open water around the text (as in the 512×128 run)

heros = "TeX Gyre Heros Makie Bold"
first_options = [
    ("A. Current: Heros Bold, one line", [Line(text = "UMD Ocean Dynamics", font = heros)]),
    ("B. Arial Black caps, tracked", [Line(text = "UMD OCEAN DYNAMICS", font = "Arial Black", tracking = 0.08)]),
    ("C. Stacked: UMD over OCEAN DYNAMICS (Arial Black / Avenir Next Bold)",
        [Line(text = "UMD", font = "Arial Black", size = 520, y = 120, tracking = 0.06, width = 2600),
         Line(text = "OCEAN DYNAMICS", font = "Avenir Next Bold", size = 200, y = -260, tracking = 0.12, width = 2600)]),
    ("D. Stacked mixed case (Avenir Next Bold)",
        [Line(text = "UMD", font = "Avenir Next Bold", size = 520, y = 120, tracking = 0.04, width = 2000),
         Line(text = "Ocean Dynamics", font = "Avenir Next Bold", size = 200, y = -260, tracking = 0.04, width = 2600)]),
    ("E. Wave baseline (Arial Rounded Bold)",
        [Line(text = "UMD Ocean Dynamics", font = "Arial Rounded MT Bold", tracking = 0.04,
              amplitude = 0.03, wavelength = 0.55)]),
    ("F. Three-line block (DIN Alternate Bold)",
        [Line(text = "UMD", font = "DIN Alternate Bold", y = 330, tracking = 0.25, width = 1500),
         Line(text = "OCEAN", font = "DIN Alternate Bold", y = -10, tracking = 0.12, width = 1500),
         Line(text = "DYNAMICS", font = "DIN Alternate Bold", y = -350, tracking = 0.03, width = 1500)]),
]


block_options = [
    ("Futura Bold",            block(face("$SUP/Futura.ttc", 2))),
    ("Avenir Next Demi Bold",  block(face("$SYS/Avenir Next.ttc", 2))),
    ("Gill Sans SemiBold",     block(face("$SUP/GillSans.ttc", 4))),
    ("Optima Bold",            block(face("$SYS/Optima.ttc", 1))),
    ("Helvetica Neue Bold",    block(face("$SYS/HelveticaNeue.ttc", 1))),
    ("Seravek Bold",           block(face("$SUP/Seravek.ttc", 9))),
]

# "Bubble" letters: blur the glyph mask and threshold below 0.5. This fattens every stroke and
# rounds every corner; σ is in units of the font's own stroke width, lower `threshold` = fatter.
function downsample(B, n)
    nx, ny = size(B) .÷ n
    return [sum(@view B[(i-1)*n+1:i*n, (j-1)*n+1:j*n]) / n^2 for i in 1:nx, j in 1:ny]
end

function inflate(B; σ, threshold, stroke)
    pad = ceil(Int, 4σ * stroke)
    P = zeros(size(B) .+ 2pad); P[pad+1:end-pad, pad+1:end-pad] .= B
    n = max(1, round(Int, stroke / 8))                 # keep ~8 px per stroke for speed
    A = gaussian_smooth(downsample(P, n), σ * stroke / n)
    M = Float64.(A .> threshold)
    xs = findall(vec(maximum(M, dims = 2)) .> 0.5); ys = findall(vec(maximum(M, dims = 1)) .> 0.5)
    return M[first(xs):last(xs), first(ys):last(ys)]
end

rounded = "Arial Rounded MT Bold"
synthetic_bubble_options = [
    ("Arial Rounded Bold, light inflation",  block(rounded; umd_width = 0.45, tracking = 0.08), (σ = 0.35, threshold = 0.3)),
    ("Arial Rounded Bold, strong inflation", block(rounded; umd_width = 0.45, tracking = 0.08), (σ = 0.6, threshold = 0.2)),
    ("Futura Bold, inflated",                block(face("$SUP/Futura.ttc", 2); umd_width = 0.45, tracking = 0.08), (σ = 0.5, threshold = 0.25)),
    ("Avenir Next Heavy, inflated",          block(face("$SYS/Avenir Next.ttc", 8); umd_width = 0.45, tracking = 0.08), (σ = 0.5, threshold = 0.25)),
    ("Chalkboard SE Bold, inflated",         block(face("$SUP/ChalkboardSE.ttc", 2); umd_width = 0.45, tracking = 0.08), (σ = 0.5, threshold = 0.25)),
    ("Arial Rounded Bold, one line, strong", [Line(text = "UMD Ocean Dynamics", font = rounded, tracking = 0.1)], (σ = 0.6, threshold = 0.2)),
]

# Rounded / bubble display fonts from Google Fonts (see fonts/*/ for licenses)
bubble_block(file) = block(FTFont(fontfile(file)); umd_width = 0.45, tracking = 0.06)
bubble_options = [
    ("Fredoka Bold",        bubble_block("fredoka/Fredoka-Bold.ttf")),
    ("Baloo 2 ExtraBold",   bubble_block("baloo2/Baloo2-ExtraBold.ttf")),
    ("Titan One",           bubble_block("titanone/TitanOne-Regular.ttf")),
    ("Lilita One",          bubble_block("lilitaone/LilitaOne-Regular.ttf")),
    ("Chewy",               bubble_block("chewy/Chewy-Regular.ttf")),
    ("Luckiest Guy",        bubble_block("luckiestguy/LuckiestGuy-Regular.ttf")),
    ("Sniglet ExtraBold",   bubble_block("sniglet/Sniglet-ExtraBold.ttf")),
    ("Arial Rounded Bold",  block(rounded; umd_width = 0.45, tracking = 0.06)),
]

# Bubble Vibes (1001fonts, personal use only; kept out of git)
vibes = "bubblevibes/BubbleVibes-Regular.otf"
vibes_font() = FTFont(fontfile(vibes))
vibes_options = [
    ("Bubble Vibes: natural OCEAN, UMD 70%",          block(vibes_font(); ocean_natural = true, umd_width = 0.35, tracking = 0.02, gap = 0.2)),
    ("Bubble Vibes: natural OCEAN, UMD 85%",          block(vibes_font(); ocean_natural = true, umd_scale = 0.85, umd_width = 0.42, tracking = 0.02, gap = 0.2)),
    ("Bubble Vibes: natural OCEAN, tighter, UMD 70%", block(vibes_font(); ocean_natural = true, umd_width = 0.32, tracking = -0.02, gap = 0.12)),
    ("Bubble Vibes: natural OCEAN, UMD 100%",         block(vibes_font(); ocean_natural = true, umd_scale = 1.0, umd_width = 0.5, tracking = 0.02, gap = 0.2)),
]

vibes_two_line_options = [
    ("Bubble Vibes: UMD / OCEAN DYNAMICS, UMD 100%",  two_line(vibes_font())),
    ("Bubble Vibes: UMD / OCEAN DYNAMICS, UMD 150%",  two_line(vibes_font(); umd_scale = 1.5, gap = 0.15)),
    ("Bubble Vibes: UMD / Ocean Dynamics, UMD 150%",  two_line(vibes_font(); line2 = "Ocean Dynamics", umd_scale = 1.5, gap = 0.15)),
    ("Bubble Vibes: UMD / OCEAN DYNAMICS, UMD 200%, tight", two_line(vibes_font(); umd_scale = 2.0, tracking = -0.02, gap = 0.1)),
]

set = isempty(ARGS) ? "first" : ARGS[1]
# Optional 2nd argument: scale every option to the same total height (cells) instead of by stroke width
block_height = length(ARGS) ≥ 2 ? parse(Int, ARGS[2]) : nothing
options, outfile = set == "block"  ? (block_options, "output/logo_options_block.png") :
                    set == "bubble" ? (bubble_options, "output/logo_options_bubble.png") :
                    set == "vibes"  ? (vibes_options, "output/logo_options_vibes.png") :
                    set == "vibes2" ? (vibes_two_line_options, "output/logo_options_vibes2.png") :
                    set == "inflated" ? (synthetic_bubble_options, "output/logo_options_inflated.png") :
                                      (first_options, "output/logo_options.png")

ocean = RGBf(0.05, 0.12, 0.25)
fig = Figure(size = (1800, 500 * cld(length(options), 2)))
for (k, opt) in enumerate(options)
    label, lines = opt[1], opt[2]
    B = layout_bitmap(lines)
    # Resolution is set by the thinnest line so every line gets ≥ cells_per_stroke cells per stroke
    stroke = minimum(median_stroke(layout_bitmap([l])) for l in lines)
    if length(opt) == 3
        B = inflate(B; opt[3]..., stroke)
        # inflate() returns a downsampled bitmap, so re-measure the stroke on inflated single lines
        stroke = minimum(median_stroke(inflate(layout_bitmap([l]); opt[3]..., stroke)) for l in lines)
    end
    scale = isnothing(block_height) ? stroke / cells_per_stroke : size(B, 2) / block_height
    stroke_cells = stroke / scale
    M = rasterize(B, scale; margins)
    Nx, Ny = size(M)
    cost = Nx * Ny / (512 * 128)
    ax = Axis(fig[(k + 1) ÷ 2, 2 - k % 2]; aspect = DataAspect(), backgroundcolor = ocean,
              title = @sprintf("%s\ngrid %d × %d (%.1f× current cost), stroke ≈ %.1f cells", label, Nx, Ny, cost, stroke_cells))
    hidedecorations!(ax); hidespines!(ax)
    draw_letters!(ax, 1:Nx, 1:Ny, M; outline = 1.5)
    limits!(ax, 0.5, Nx + 0.5, 0.5, Ny + 0.5)
    @info @sprintf("%-60s grid %4d × %4d  cost %.2f×", label, Nx, Ny, cost)
end
resize_to_layout!(fig)
mkpath("output")
isnothing(block_height) || (outfile = replace(outfile, ".png" => "_h$(block_height).png"))
save(outfile, fig)
@info "Saved $outfile"
