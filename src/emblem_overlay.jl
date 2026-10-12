# Overlay the styled emblem + OCEAN DYNAMICS logo (land shadow included) on an axis whose
# coordinates are km, for a run built with `emblem_text_masks` (plotting only).
using CairoMakie.Colors: RGB

# Color schemes from output/emblem_colors.png: (fills (OCEAN, DYNAMICS), outline, shadow, outline width, split words?)
const LOGO_RED = RGB(0.927, 0.106, 0.18)
const EMBLEM_SCHEMES = Dict(
    :red       => ((LOGO_RED, LOGO_RED), RGB(0.08, 0.08, 0.08), RGB(0.08, 0.08, 0.08), 2.0, false),
    :white     => ((RGB(1, 1, 1), RGB(1, 1, 1)), RGB(0.08, 0.08, 0.08), RGB(0.08, 0.08, 0.08), 2.0, false),
    :red_rim   => ((LOGO_RED, LOGO_RED), RGB(1, 1, 1), RGB(0.08, 0.08, 0.08), 3.0, false),
    :two_tone  => ((RGB(1, 1, 1), LOGO_RED), RGB(0.08, 0.08, 0.08), RGB(0.08, 0.08, 0.08), 2.0, true))

"""
    emblem_overlay_data(; cells_per_Ld, Δx_km, x₀_km, y₀_km, r = 4)

Fine-resolution masks and emblem box (in km) for the logo layout (src/logo_config.jl), consistent
with the simulation mask. (x₀_km, y₀_km) is the domain's lower-left corner.
"""
function emblem_overlay_data(; cells_per_Ld, Δx_km, x₀_km, y₀_km, r = 4)
    f = logo_layout(cells_per_Ld; r)
    xs = x₀_km .+ ((1:f.grid[1]) .- 0.5) ./ r .* Δx_km
    ys = y₀_km .+ ((1:f.grid[2]) .- 0.5) ./ r .* Δx_km
    bx, by, bw, bh = f.emblem_box
    box = (x₀_km + bx / r * Δx_km, y₀_km + by / r * Δx_km, bw / r * Δx_km, bh / r * Δx_km)
    shadow = f.land .& .!f.letters
    i₁, i₂ = clamp(ceil(Int, bx), 1, f.grid[1]), clamp(floor(Int, bx + bw), 1, f.grid[1])
    shadow[i₁:i₂, :] .= false                                  # the emblem region shows its artwork
    return (; xs = collect(xs), ys = collect(ys), letters = f.letters, shadow, img = f.img, box, r,
            stroke_px = STROKE_PER_LD * cells_per_Ld * r)
end

function draw_emblem_logo!(ax, d; scheme = :red)
    fills, outc, shc, ow, split = EMBLEM_SCHEMES[scheme]
    cols = [any(d.letters[i, :]) for i in axes(d.letters, 1)]
    on = findall(cols); gaps = [(on[m+1] - on[m], m) for m in 1:length(on)-1]
    _, gm = maximum(gaps); sx = (on[gm] + on[gm+1]) ÷ 2
    draw_styled_letters!(ax, d.xs, d.ys, d.letters; style = :retro_flat, stroke_px = d.stroke_px, fills,
                         outline_color = outc, shadow_color = shc, outline_width = ow,
                         shadow = d.shadow, split_x = split ? sx : nothing)
    x, y, w, h = d.box
    image!(ax, (x, x + w), (y, y + h), d.img)
    return nothing
end
