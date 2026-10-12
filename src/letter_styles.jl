# Decorative rendering of the logo islands for figures and animations (plotting only).
# Letters are drawn from a finer-resolution mask of the same layout, so outlines are smooth;
# the simulation still sees the grid-resolution mask.
using CairoMakie, Statistics
using CairoMakie.Colors: RGBA, RGB, @colorant_str, weighted_color_mean

# Euclidean-ish distance (in pixels) from each land pixel to the nearest water pixel (two-pass chamfer)
function inside_distance(M)
    nx, ny = size(M)
    d = [M[i, j] ? Inf : 0.0 for i in 1:nx, j in 1:ny]
    a, b = 1.0, sqrt(2.0)
    for j in 1:ny, i in 1:nx
        M[i, j] || continue
        i > 1 && (d[i, j] = min(d[i, j], d[i-1, j] + a))
        j > 1 && (d[i, j] = min(d[i, j], d[i, j-1] + a))
        i > 1 && j > 1 && (d[i, j] = min(d[i, j], d[i-1, j-1] + b))
        i < nx && j > 1 && (d[i, j] = min(d[i, j], d[i+1, j-1] + b))
    end
    for j in ny:-1:1, i in nx:-1:1
        M[i, j] || continue
        i < nx && (d[i, j] = min(d[i, j], d[i+1, j] + a))
        j < ny && (d[i, j] = min(d[i, j], d[i, j+1] + a))
        i < nx && j < ny && (d[i, j] = min(d[i, j], d[i+1, j+1] + b))
        i > 1 && j < ny && (d[i, j] = min(d[i, j], d[i-1, j+1] + b))
    end
    return d
end

# Connected land components (4-connectivity), numbered left to right within top and bottom halves
function letter_components(M)
    nx, ny = size(M)
    lab = zeros(Int, nx, ny); n = 0
    for j in 1:ny, i in 1:nx
        (M[i, j] && lab[i, j] == 0) || continue
        n += 1; lab[i, j] = n; stack = [(i, j)]
        while !isempty(stack)
            a, c = pop!(stack)
            for (da, dc) in ((1, 0), (-1, 0), (0, 1), (0, -1))
                a2, c2 = a + da, c + dc
                (1 ≤ a2 ≤ nx && 1 ≤ c2 ≤ ny && M[a2, c2] && lab[a2, c2] == 0) || continue
                lab[a2, c2] = n; push!(stack, (a2, c2))
            end
        end
    end
    cents = [(mean(I[1] for I in findall(==(k), lab)), mean(I[2] for I in findall(==(k), lab))) for k in 1:n]
    ymid = mean(c[2] for c in cents)
    order = sortperm(cents; by = c -> (c[2] < ymid, c[1]))     # top line first, then left to right
    rank = zeros(Int, n); rank[order] .= 1:n
    return map(l -> l == 0 ? 0 : rank[l], lab), n
end

"""
    draw_styled_letters!(ax, x, y, M; style, stroke_px, shadow = nothing)

Draw land mask `M` (pixel centers at `x`, `y`) as decorated letters. Styles:
`:plain`, `:retro`, `:candy`, `:bubblegum`, `:sunset`, `:aqua`.
"""
function draw_styled_letters!(ax, x, y, M; style = :retro, stroke_px = 40, shadow = nothing,
                              fills = (RGB(1.0, 0.60, 0.68), RGB(0.97, 0.78, 0.55)),
                              outline_color = nothing, umd_pattern = :none, outline_width = nothing,
                              shadow_color = nothing, split_x = nothing)
    style == :plain && return begin
        heatmap!(ax, x, y, ifelse.(M, 1.0, NaN); colormap = [:white], colorrange = (0, 1))
        contour!(ax, x, y, Float64.(M); levels = [0.5], color = :black, linewidth = 1.5)
    end
    nx, ny = size(M)
    lab, n = letter_components(M)

    # Puffy height field: rounded profile rising over half a stroke from the edge
    d = inside_distance(M)
    R = stroke_px / 2
    h = [M[i, j] ? R * sqrt(1 - (1 - min(d[i, j] / R, 1))^2) : 0.0 for i in 1:nx, j in 1:ny]
    hx = similar(h); hy = similar(h)
    for j in 1:ny, i in 1:nx
        hx[i, j] = (h[min(i+1, nx), j] - h[max(i-1, 1), j]) / 2
        hy[i, j] = (h[i, min(j+1, ny)] - h[i, max(j-1, 1)]) / 2
    end
    L = [-0.5, 0.6, 0.62]; L ./= sqrt(sum(abs2, L))                       # light from the upper left
    Hv = L .+ [0, 0, 1]; Hv ./= sqrt(sum(abs2, Hv))

    # Vertical position within each text line (0 = bottom, 1 = top), for per-line gradients
    rows = [any(M[:, j]) for j in 1:ny]
    yl = zeros(ny); runs = UnitRange{Int}[]; j = 1
    while j ≤ ny
        if rows[j]; j2 = j; while j2 < ny && rows[j2+1]; j2 += 1; end; push!(runs, j:j2); j = j2 + 1; else j += 1; end
    end
    for rr in runs, jj in rr; yl[jj] = (jj - first(rr)) / max(1, length(rr) - 1); end

    base(i, j) = begin
        k = lab[i, j]; yfrac = (j - 1) / (ny - 1)
        style in (:retro, :retro_flat) && k ≤ 3 && umd_pattern == :calvert ? calvert(i, j) :
        style in (:retro, :retro_flat) && !isnothing(split_x) ? (i ≤ split_x ? fills[1] : fills[2]) :   # by position
        style in (:retro, :retro_flat) ? (k ≤ 3 ? fills[1] : fills[2]) :          # (UMD, OCEAN DYNAMICS)
        style in (:sunset_lines, :sunset_gloss) ? weighted_color_mean(yl[j], colorant"#FFE66D", colorant"#FF4F87") :
        style == :terps     ? (k ≤ 3 ? colorant"#E21833" : colorant"#FFD200") :
        style == :candy     ? [colorant"#FF5E8A", colorant"#FFA62B", colorant"#FFD93D", colorant"#6BCB77",
                               colorant"#4D96FF", colorant"#B983FF"][mod1(k, 6)] :
        style == :bubblegum ? RGB(1.0, 0.56, 0.78) :
        style == :sunset    ? weighted_color_mean(clamp(yfrac, 0, 1), colorant"#FFE66D", colorant"#FF4F87") :
                              weighted_color_mean(clamp(yfrac, 0, 1), colorant"#B5FFF1", colorant"#1FA2D6")   # :aqua
    end
    # Simplified Calvert arms (Maryland flag): paly gold and black, counterchanged across a bend
    umd_px = findall(l -> 1 ≤ l ≤ 3, lab)
    ci = isempty(umd_px) ? (0, 0, 1, 1) : (minimum(I[1] for I in umd_px), minimum(I[2] for I in umd_px),
                                           maximum(I[1] for I in umd_px), maximum(I[2] for I in umd_px))
    calvert(i, j) = begin
        w = (ci[3] - ci[1] + 1) / 12                          # 12 stripes across "UMD"
        stripe = isodd(floor(Int, (i - ci[1]) / w))
        above = (j - ci[2]) / (ci[4] - ci[2] + 1) > 1 - (i - ci[1]) / (ci[3] - ci[1] + 1)   # bend from top-left
        xor(stripe, above) ? RGB(0.08, 0.08, 0.08) : RGB(1.0, 0.82, 0.0)
    end
    glossy = style in (:candy, :bubblegum, :aqua, :sunset_gloss)
    flat = style == :retro_flat

    img = fill(RGBA(0.0, 0.0, 0.0, 0.0), nx, ny)
    for j in 1:ny, i in 1:nx
        M[i, j] || continue
        nrm = [-hx[i, j], -hy[i, j], 1.0]; nrm ./= sqrt(sum(abs2, nrm))
        diffuse = max(0, sum(nrm .* L))
        spec = glossy ? max(0, sum(nrm .* Hv))^40 : max(0, sum(nrm .* Hv))^80 * 0.4
        c = base(i, j)
        s = 0.55 + 0.55 * diffuse
        if flat   # mostly flat fill with a soft rim shade and a small highlight, like printed bubble letters
            s = 0.88 + 0.12 * diffuse
            spec = max(0, sum(nrm .* Hv))^120 * 0.5
        end
        r, g, b = clamp(c.r * s + spec, 0, 1), clamp(c.g * s + spec, 0, 1), clamp(c.b * s + spec, 0, 1)
        img[i, j] = RGBA(r, g, b, 1.0)
    end

    dx, dy = x[2] - x[1], y[2] - y[1]
    outline = !isnothing(outline_color) ? outline_color :
              style in (:retro, :retro_flat, :sunset, :sunset_lines, :sunset_gloss) ? RGB(0.42, 0.13, 0.20) :
              style == :terps ? RGB(0.05, 0.05, 0.05) : style == :aqua ? RGB(0.03, 0.20, 0.35) : RGB(1.0, 1.0, 1.0)

    # Retro offset shadow (down and to the right), drawn under the letters. Pass `shadow` (a mask on
    # the same pixels) to draw a shadow that is part of the land instead of a decorative offset.
    if style in (:retro, :retro_flat, :sunset, :sunset_lines, :sunset_gloss, :terps)
        sh = round(Int, 0.35stroke_px)
        S = falses(nx, ny); S[1+sh:end, 1:end-sh] .= M[1:end-sh, 1+sh:end]
        isnothing(shadow) || (S = shadow)
        shc = isnothing(shadow_color) ? outline : shadow_color
        image!(ax, (x[1] - dx/2, x[end] + dx/2), (y[1] - dy/2, y[end] + dy/2),
               [S[i, j] && !M[i, j] ? RGBA(shc.r, shc.g, shc.b, 1.0) : RGBA(0, 0, 0, 0) for i in 1:nx, j in 1:ny];
               interpolate = false)
        contour!(ax, x, y, Float64.(S .| M); levels = [0.5], color = shc,
                 linewidth = isnothing(outline_width) ? 1.5 : outline_width)
    end
    image!(ax, (x[1] - dx/2, x[end] + dx/2), (y[1] - dy/2, y[end] + dy/2), img; interpolate = false)
    contour!(ax, x, y, Float64.(M); levels = [0.5], color = outline,
             linewidth = !isnothing(outline_width) ? outline_width : style == :candy ? 2.5 : 1.8)

    # Sparkles on a few letters (as in the font's sample art)
    if style in (:aqua, :bubblegum, :candy)
        tops = [findmax(j -> any(lab[:, j] .== k) ? j : 0, 1:ny) for k in 1:n]
        for k in (1, 3, 5, 9, 12)
            k ≤ n || continue
            I = findall(==(k), lab)
            j★ = maximum(I2[2] for I2 in I); i★ = round(Int, mean(I2[1] for I2 in I if I2[2] == j★))
            scatter!(ax, [x[i★] + 0.15stroke_px * dx], [y[j★] - 0.15stroke_px * dy]; marker = '✦',
                     markersize = 18, color = :white)
        end
    end
    return nothing
end
