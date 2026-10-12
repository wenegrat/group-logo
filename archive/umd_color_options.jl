# Retro-flat letters with alternative UMD colors chosen to match the visual weight
# (CIELAB lightness and chroma) of the peach OCEAN DYNAMICS line. Plotting only.
# Usage: julia --project=. scripts/umd_color_options.jl [name] [frame]
using Oceananigans, CairoMakie, Printf
using CairoMakie.Colors: Lab, LCHab, RGB, convert
include(joinpath(@__DIR__, "..", "src", "Logo.jl"))
include(joinpath(@__DIR__, "..", "src", "letter_styles.jl"))

name = length(ARGS) ≥ 1 ? ARGS[1] : "eady_vibes_448x128x8"
ζs = FieldTimeSeries("output/$(name)_surface.jld2", "ζ"; backend = OnDisk())
n = length(ARGS) ≥ 2 ? parse(Int, ARGS[2]) : length(ζs.times)
grid = ζs.grid; ug = grid.underlying_grid
Nx, Ny = size(grid, 1), size(grid, 2)
Δx = xspacings(ug, Center())[1] / 1e3
xf = xnodes(ug, Face())[1:Nx] ./ 1e3; yf = ynodes(ug, Face())[1:Ny] ./ 1e3
ζ = interior(ζs[n], 1:Nx, 1:Ny, 1) ./ 1e-4

r = 4; cells_per_stroke = 10
vibes = FTFont(fontfile("bubblevibes/BubbleVibes-Regular.otf"))
Mf = logo_mask(two_line(vibes); cells_per_stroke = r * cells_per_stroke, grid = (r * Nx, r * Ny))
xfine = collect(range(xf[1] + Δx / (2r), step = Δx / r, length = r * Nx))
yfine = collect(range(yf[1] + Δx / (2r), step = Δx / r, length = r * Ny))
navy = cgrad([RGBf(0.02, 0.06, 0.18), RGBf(0.10, 0.30, 0.55), RGBf(0.75, 0.88, 0.95)])

peach = RGB(0.97, 0.78, 0.55)
lch(c) = convert(LCHab, c)
P = lch(peach)
# Same lightness and chroma as the peach, different hue
samelc(h) = convert(RGB, LCHab(P.l, P.c, h))
lighter(c, dL) = (L = lch(c); convert(RGB, LCHab(L.l + dL, L.c, L.h)))
maroon, ink = RGB(0.42, 0.13, 0.20), RGB(0.08, 0.08, 0.08)
md_red, md_gold = RGB(0.886, 0.094, 0.20), RGB(1.0, 0.824, 0.0)        # Maryland flag red #E21833, gold #FFD200
# Pastel Maryland colors: the flag hues at retro-like lightness and chroma
pastel(c, l, ch) = convert(RGB, LCHab(l, ch, lch(c).h))
# (label, UMD fill, OCEAN DYNAMICS fill, outline/shadow color, UMD pattern)
options = [("soft pink",                        samelc(5),          peach, maroon, :none),
           ("rose",                             samelc(20),         peach, maroon, :none),
           ("coral",                            samelc(40),         peach, maroon, :none),
           ("salmon (between coral and apricot)", samelc(50),       peach, maroon, :none),
           ("apricot",                          samelc(60),         peach, maroon, :none),
           ("coral, a little lighter",          lighter(samelc(40), 4), peach, maroon, :none),
           ("Maryland pastel: red / gold",      pastel(md_red, 68, 50), pastel(md_gold, 86, 60), ink, :none),
           ("Maryland bold: red / gold",        md_red,             md_gold, ink, :none),
           ("Maryland flag: Calvert UMD / gold", md_gold,           pastel(md_gold, 86, 60), ink, :calvert)]
for (label, c, c2, _, _) in options
    L = lch(c); @info @sprintf("%-38s UMD L* = %5.1f C* = %5.1f h = %5.1f | line 2 L* = %5.1f", label, L.l, L.c, L.h, lch(c2).l)
end

clip(c) = RGB(clamp(c.r, 0, 1), clamp(c.g, 0, 1), clamp(c.b, 0, 1))
fig = Figure(size = (1800, 1300))
for (k, (label, c, c2, oc, pat)) in enumerate(options)
    L = lch(c)
    ax = Axis(fig[(k - 1) ÷ 3 + 1, (k - 1) % 3 + 1]; aspect = DataAspect(),
              title = pat == :none ? @sprintf("%s  (UMD L* %.0f, C* %.0f)", label, L.l, L.c) : label)
    hidedecorations!(ax); hidespines!(ax)
    heatmap!(ax, xf, yf, abs.(ζ); colormap = navy, colorrange = (0, 2))
    draw_styled_letters!(ax, xfine, yfine, Mf; style = :retro_flat, stroke_px = r * cells_per_stroke,
                         fills = (clip(c), clip(c2)), outline_color = oc, umd_pattern = pat)
end
resize_to_layout!(fig)
save("output/umd_color_options2.png", fig)
@info "Saved output/umd_color_options2.png  (OCEAN DYNAMICS peach: L* $(round(P.l, digits=1)), C* $(round(P.c, digits=1)))"
