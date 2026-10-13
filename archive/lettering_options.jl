# Compare lettering for "OCEAN DYNAMICS" next to the emblem (plotting only; nothing is simulated).
# Each option is shown over a strip of a background image (e.g. the website's header photo) as it
# would be drawn in the figures: flat fill with a thin outline. For each, the grid needed at the
# run's resolution (7.5 cells per Ld, median stroke = 2 Ld) and the width of its thinnest strokes
# are reported, since thin strokes and serifs still have to be resolved as land.
#
# Usage: julia --project=. scripts/lettering_options.jl [background image]
using CairoMakie, Printf, Statistics
using CairoMakie.Colors: RGB, RGBA
include(joinpath(@__DIR__, "..", "src", "logo_config.jl"))
include(joinpath(@__DIR__, "..", "src", "letter_styles.jl"))

bg_path = isempty(ARGS) ? joinpath(REPO, "output", "website_mock", "site", "img", "Landsat.jpg") : ARGS[1]
cells_per_Ld, r = 7.5, 3

# (label, font, tracking)
options = [("Montserrat ExtraBold (the website's font)", FTFont(fontfile("montserrat/Montserrat-ExtraBold.ttf")), 0.06),
           ("Futura Bold",                         face("$SUP/Futura.ttc", 2), 0.08),
           ("Avenir Next Demi Bold",               face("$SYS/Avenir Next.ttc", 2), 0.08),
           ("Optima Bold",                         face("$SYS/Optima.ttc", 1), 0.08),
           ("Cinzel ExtraBold (Roman capitals)",   FTFont(fontfile("cinzel/Cinzel-ExtraBold.ttf")), 0.06),
           ("Playfair Display ExtraBold",          FTFont(fontfile("playfairdisplay/PlayfairDisplay-ExtraBold.ttf")), 0.05),
           ("DM Serif Display",                    FTFont(fontfile("dmserifdisplay/DMSerifDisplay-Regular.ttf")), 0.05),
           ("Baskerville Bold",                    face("$SUP/Baskerville.ttc", 1), 0.06)]

bg = CairoMakie.FileIO.load(bg_path)
bg = permutedims(reverse(bg, dims = 1))                     # [x, y], y up

function thin_cells(land)                                   # 10th percentile of land run lengths (both directions)
    runs = Int[]
    for M in (land, permutedims(land)), j in axes(M, 2)
        n = 0
        for i in axes(M, 1); if M[i, j]; n += 1 elseif n > 0; push!(runs, n); n = 0 end; end
    end
    return sort(runs)[max(1, round(Int, 0.1length(runs)))]
end

fig = Figure(size = (1800, 1700))
for (k, (label, font, tracking)) in enumerate(options)
    args = (; emblem_path = emblem_path(), font, cells_per_stroke = STROKE_PER_LD * cells_per_Ld,
            margins = round.(Int, MARGINS_IN_LD .* cells_per_Ld), shadow = nothing, text = "OCEAN DYNAMICS", tracking)
    m = emblem_text_masks(; args...)
    f = emblem_text_masks(; args..., r)
    stats = @sprintf("grid %d × %d at 7.5 cells/Ld; thinnest strokes ≈ %d cells", m.grid..., thin_cells(m.land))
    @info "$label: $stats"

    ax = Axis(fig[(k + 1) ÷ 2, 2 - k % 2]; aspect = DataAspect(), title = "$label\n$stats")
    hidedecorations!(ax); hidespines!(ax)
    # Background: a strip of the photo with the domain's aspect ratio
    bw = size(bg, 1); bh = min(size(bg, 2), round(Int, bw * m.grid[2] / m.grid[1])); j₀ = (size(bg, 2) - bh) ÷ 2
    image!(ax, (0, m.grid[1]), (0, m.grid[2]), bg[:, j₀+1:j₀+bh])
    xs = collect(1:f.grid[1]) ./ r; ys = collect(1:f.grid[2]) ./ r
    draw_styled_letters!(ax, xs, ys, f.letters; style = :flat, stroke_px = r * STROKE_PER_LD * cells_per_Ld,
                         fills = (RGB(1, 1, 1), RGB(1, 1, 1)), outline_color = RGB(0.1, 0.1, 0.1), outline_width = 1.2)
    x₀, y₀, w, h = f.emblem_box ./ r
    image!(ax, (x₀, x₀ + w), (y₀, y₀ + h), f.img)
    limits!(ax, 0, m.grid[1], 0, m.grid[2])
end
resize_to_layout!(fig)
mkpath("output")
save("output/lettering_options.png", fig)
@info "Saved output/lettering_options.png"
