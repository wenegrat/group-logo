# Logo layouts: render lines of text with a given font, spacing and arrangement to a bitmap,
# and rasterize it onto a simulation grid as a land mask.
using CairoMakie, FreeTypeAbstraction, Statistics
const FTA = FreeTypeAbstraction

fontfile(name) = joinpath(@__DIR__, "..", "fonts", name)

# A line of text: string, font, fontsize (canvas px), center height, tracking (fraction of em),
# and an optional sine-wave baseline (amplitude and wavelength in units of the line width).
Base.@kwdef struct Line
    text :: String
    font :: Any                                  # font name or FTFont
    size :: Float64 = 200
    y    :: Float64 = 0
    tracking  :: Float64 = 0
    amplitude :: Float64 = 0
    wavelength :: Float64 = 1
    width :: Union{Nothing, Float64} = nothing   # if set, rescale `size` so the line spans this width
end

advances(font, line, size) = [FreeTypeAbstraction.hadvance(FreeTypeAbstraction.get_extent(font, c)) * size +
                              line.tracking * size for c in line.text]

# Load face `i` of a font collection (.ttc); Makie's name lookup only sees face 0
face(file, i) = (r = FTA.newface_mmapped(file, i); FTFont(r[1], true, r[2]))
const SYS, SUP = "/System/Library/Fonts", "/System/Library/Fonts/Supplemental"

function draw_line!(ax, line)
    font = Makie.to_font(line.font)
    size = line.size
    if !isnothing(line.width)
        size *= line.width / sum(advances(font, line, size))
    end
    adv = advances(font, line, size)
    W = sum(adv)
    xs = -W / 2 .+ cumsum(adv) .- adv ./ 2
    k = 2π / (line.wavelength * W)
    ys = line.y .+ line.amplitude * W .* sin.(k .* xs)
    θ = atan.(line.amplitude * W * k .* cos.(k .* xs))
    for (c, x, y, r) in zip(line.text, xs, ys, θ)
        c == ' ' && continue
        text!(ax, x, y; text = string(c), font, fontsize = size, rotation = r,
              align = (:center, :baseline), color = :white, markerspace = :data)
    end
end

# Three-line block: UMD / OCEAN / DYNAMICS. OCEAN and DYNAMICS share a font size; OCEAN is
# letter-spaced to DYNAMICS' width unless `ocean_natural`, in which case it keeps `tracking`; UMD has cap height `umd_scale` × theirs, spaced to `umd_width`.
function block(font; umd_scale = 0.7, umd_width = 0.6, W = 1500.0, tracking = 0.04, gap = 0.35,
               ocean_natural = false)
    font = Makie.to_font(font)
    unit(text, tr) = sum(FTA.hadvance(FTA.get_extent(font, c)) + tr for c in text)
    capH = FTA.get_extent(font, 'H').scale[2]
    s = W / unit("DYNAMICS", tracking)
    tr_ocean = ocean_natural ? tracking : tracking + (W / s - unit("OCEAN", tracking)) / length("OCEAN")
    s_umd = umd_scale * s
    tr_umd = tracking + (umd_width * W / s_umd - unit("UMD", tracking)) / length("UMD")
    y_umd = 0.0
    y_ocean = y_umd - capH * s - gap * capH * s_umd
    y_dyn = y_ocean - capH * s * (1 + gap)
    return [Line(text = "UMD", font = font, size = s_umd, y = y_umd, tracking = tr_umd),
            Line(text = "OCEAN", font = font, size = s, y = y_ocean, tracking = tr_ocean),
            Line(text = "DYNAMICS", font = font, size = s, y = y_dyn, tracking = tracking)]
end

# Two lines: UMD over `line2`, both at natural spacing; UMD cap height is `umd_scale` × line 2's.
function two_line(font; line2 = "OCEAN DYNAMICS", umd_scale = 1.0, tracking = 0.02, gap = 0.2, size = 200.0)
    font = Makie.to_font(font)
    capH = FTA.get_extent(font, 'H').scale[2]
    s_umd = umd_scale * size
    return [Line(text = "UMD", font = font, size = s_umd, y = 0.0, tracking = tracking),
            Line(text = line2, font = font, size = size, y = -capH * (size + gap * min(size, s_umd)), tracking = tracking)]
end

function layout_bitmap(lines; canvas = (6000, 2400))
    fig = Figure(size = canvas, backgroundcolor = :black, figure_padding = 0)
    ax = Axis(fig[1, 1], backgroundcolor = :black)
    hidedecorations!(ax); hidespines!(ax)
    limits!(ax, -canvas[1] / 2, canvas[1] / 2, -canvas[2] / 2, canvas[2] / 2)
    foreach(l -> l isa Swash ? draw_swash!(ax, l) : draw_line!(ax, l), lines)
    img = Makie.colorbuffer(fig; px_per_unit = 1)
    A = permutedims(reverse(Float64.(Makie.Colors.Gray.(img)), dims = 1))
    xs = findall(vec(maximum(A, dims = 2)) .> 0.5); ys = findall(vec(maximum(A, dims = 1)) .> 0.5)
    return A[first(xs):last(xs), first(ys):last(ys)]
end

function median_stroke(B)
    runs = Int[]
    for j in axes(B, 2)
        n = 0
        for i in axes(B, 1)
            if B[i, j] > 0.5; n += 1 elseif n > 0; push!(runs, n); n = 0 end
        end
    end
    return median(runs)
end

# Nearest-neighbor sample of the bitmap onto a grid with `scale` bitmap pixels per cell, with
# `margins` cells of open water on each side; grid dimensions are rounded up to multiples of 32.
function rasterize(B, scale; margins = (80, 45))
    margin_x, margin_y = margins
    nx, ny = floor.(Int, size(B) ./ scale)
    Nx = 32 * cld(nx + 2margin_x, 32); Ny = 32 * cld(ny + 2margin_y, 32)
    M = falses(Nx, Ny)
    ox, oy = (Nx - nx) ÷ 2, (Ny - ny) ÷ 2
    for j in 1:ny, i in 1:nx
        M[ox + i, oy + j] = B[clamp(round(Int, (i - 0.5) * scale), 1, size(B, 1)),
                              clamp(round(Int, (j - 0.5) * scale), 1, size(B, 2))] > 0.5
    end
    return M
end

# Union of a bitmap with a copy shifted by (ox, oy) pixels; the canvas grows to hold both.
# With `letters_only = true`, return just the unshifted letters on the same enlarged canvas.
function add_shadow(B, (ox, oy); letters_only = false)
    nx, ny = size(B)
    S = zeros(nx + abs(ox), ny + abs(oy))
    i₀, j₀ = max(0, -ox), max(0, -oy)
    S[i₀+1:i₀+nx, j₀+1:j₀+ny] .= B
    letters_only && return S
    S[i₀+ox+1:i₀+ox+nx, j₀+oy+1:j₀+oy+ny] .= max.(S[i₀+ox+1:i₀+ox+nx, j₀+oy+1:j₀+oy+ny], B)
    return S
end

# Land mask for a layout, scaled so the thinnest line's median stroke is `cells_per_stroke` cells.
# Either pad by `margins` (rounded up to multiples of 32) or center the logo in a fixed `grid = (Nx, Ny)`.
# `shadow = (sx, sy)` adds a drop shadow, offset by (sx, sy) stroke widths (sy < 0 is downward), to the land.
# `letters_only = true` returns just the letters, positioned exactly as in the shadowed mask.
# `post(B, stroke_px)` optionally transforms the rendered bitmap (e.g. slant, thicken) before scaling.
function logo_mask(lines; cells_per_stroke, margins = (80, 45), grid = nothing, shadow = nothing, letters_only = false,
                   post = nothing)
    textlines = filter(l -> l isa Line, lines)
    stroke₀ = minimum(median_stroke(layout_bitmap([l])) for l in textlines)
    stroke = isnothing(post) ? stroke₀ : minimum(median_stroke(post(layout_bitmap([l]), stroke₀)) for l in textlines)
    B = layout_bitmap(lines)
    isnothing(post) || (B = post(B, stroke₀))
    isnothing(shadow) || (B = add_shadow(B, round.(Int, shadow .* stroke); letters_only))
    scale = stroke / cells_per_stroke
    isnothing(grid) && return rasterize(B, scale; margins)
    nx, ny = floor.(Int, size(B) ./ scale)
    M = falses(grid...)
    ox, oy = (grid[1] - nx) ÷ 2, (grid[2] - ny) ÷ 2
    M₀ = rasterize(B, scale; margins = (0, 0))          # logo centered in a 32-rounded grid
    o₀x, o₀y = (size(M₀, 1) - nx) ÷ 2, (size(M₀, 2) - ny) ÷ 2
    M[ox+1:ox+nx, oy+1:oy+ny] .= M₀[o₀x+1:o₀x+nx, o₀y+1:o₀y+ny]
    return M
end

#####
##### Athletic-script treatment: slant, rising baseline, heavier strokes, swash underline
#####

# Distance (pixels, two-pass chamfer) from each pixel to the nearest pixel where `M` is true
function distance_to(M)
    nx, ny = size(M)
    d = [M[i, j] ? 0.0 : Inf for i in 1:nx, j in 1:ny]
    a, b = 1.0, sqrt(2.0)
    for j in 1:ny, i in 1:nx
        i > 1 && (d[i, j] = min(d[i, j], d[i-1, j] + a))
        j > 1 && (d[i, j] = min(d[i, j], d[i, j-1] + a))
        i > 1 && j > 1 && (d[i, j] = min(d[i, j], d[i-1, j-1] + b))
        i < nx && j > 1 && (d[i, j] = min(d[i, j], d[i+1, j-1] + b))
    end
    for j in ny:-1:1, i in nx:-1:1
        i < nx && (d[i, j] = min(d[i, j], d[i+1, j] + a))
        j < ny && (d[i, j] = min(d[i, j], d[i, j+1] + a))
        i < nx && j < ny && (d[i, j] = min(d[i, j], d[i+1, j+1] + b))
        i > 1 && j < ny && (d[i, j] = min(d[i, j], d[i-1, j+1] + b))
    end
    return d
end

# Grow every stroke by `r` pixels on each side
thicken(B, r) = (pad = ceil(Int, r) + 1;
                 P = zeros(size(B) .+ 2pad); P[pad+1:end-pad, pad+1:end-pad] .= B;
                 Float64.(distance_to(P .> 0.5) .≤ r))

# Shear (x += shear * y, more italic) then rotate counterclockwise by `angle` (rising baseline)
function slant_rotate(B; shear = 0.0, angle = 0.0)
    nx, ny = size(B)
    c, s = cos(angle), sin(angle)
    T = [c -s; s c] * [1 shear; 0 1]
    Ti = inv(T)
    corners = [T * [x, y] for x in (0, nx), y in (0, ny)]
    xmin, xmax = extrema(p[1] for p in corners); ymin, ymax = extrema(p[2] for p in corners)
    mx, my = ceil(Int, xmax - xmin), ceil(Int, ymax - ymin)
    out = zeros(mx, my)
    for j in 1:my, i in 1:mx
        p = Ti * [xmin + i - 0.5, ymin + j - 0.5]
        a, b = floor(Int, p[1]) + 1, floor(Int, p[2]) + 1
        (1 ≤ a ≤ nx && 1 ≤ b ≤ ny) && (out[i, j] = B[a, b])
    end
    return out
end

# A swash under a text line: a tapered quadratic curve from below the line's left end to its right
# end. `drop` is how far below the baseline it dips (in cap heights), `thickness` its maximum
# width (in cap heights), tapering to a point at the right.
struct Swash
    line      :: Line
    drop      :: Float64
    thickness :: Float64
end

function draw_swash!(ax, sw::Swash)
    font = Makie.to_font(sw.line.font)
    W = sum(advances(font, sw.line, sw.line.size))
    capH = FTA.get_extent(font, 'H').scale[2] * sw.line.size
    y₀ = sw.line.y
    P0 = (-0.50W, y₀ - 0.15capH); P1 = (0.0, y₀ - (0.15 + 2sw.drop) * capH); P2 = (0.52W, y₀ - 0.05capH)
    ts = range(0, 1, length = 200)
    pt(t) = ((1-t)^2 * P0[1] + 2(1-t)*t * P1[1] + t^2 * P2[1], (1-t)^2 * P0[2] + 2(1-t)*t * P1[2] + t^2 * P2[2])
    w(t) = sw.thickness * capH * (1 - t)^0.8 * min(1, 8t + 0.3)       # rounded start, tapering tail
    up, lo = Point2f[], Point2f[]
    for t in ts
        p, q = pt(max(t - 1e-3, 0)), pt(min(t + 1e-3, 1))
        tx, ty = q[1] - p[1], q[2] - p[2]; nrm = sqrt(tx^2 + ty^2)
        nxv, nyv = -ty / nrm, tx / nrm
        c = pt(t)
        push!(up, Point2f(c[1] + w(t) / 2 * nxv, c[2] + w(t) / 2 * nyv))
        push!(lo, Point2f(c[1] - w(t) / 2 * nxv, c[2] - w(t) / 2 * nyv))
    end
    poly!(ax, vcat(up, reverse(lo)); color = :white, strokewidth = 0)
end

#####
##### Image logos (e.g. a PNG emblem) combined with text
#####

# Load an RGBA image as (silhouette, image), both indexed [x, y] with y upward.
# The silhouette is alpha ≥ 0.5, with enclosed transparent holes filled.
function load_emblem(path)
    img = CairoMakie.FileIO.load(path)                       # rows top → bottom
    img = permutedims(reverse(img, dims = 1))                 # [x, y], y up
    S = [Float64(c.alpha) ≥ 0.5 for c in img]
    # Crop both to the silhouette's bounding box so sizes refer to the artwork, not the PNG canvas
    I = findall(S)
    is, js = minimum(i[1] for i in I):maximum(i[1] for i in I), minimum(i[2] for i in I):maximum(i[2] for i in I)
    return Float64.(fill_holes(S[is, js])), img[is, js]
end

# Fill regions of false that are not connected to the array border
function fill_holes(S)
    nx, ny = size(S)
    outside = falses(nx, ny)
    stack = [(i, j) for i in 1:nx, j in (1, ny) if !S[i, j]]
    append!(stack, [(i, j) for i in (1, nx), j in 1:ny if !S[i, j]])
    for (i, j) in stack; outside[i, j] = true; end
    while !isempty(stack)
        a, c = pop!(stack)
        for (da, dc) in ((1, 0), (-1, 0), (0, 1), (0, -1))
            a2, c2 = a + da, c + dc
            (1 ≤ a2 ≤ nx && 1 ≤ c2 ≤ ny && !S[a2, c2] && !outside[a2, c2]) || continue
            outside[a2, c2] = true; push!(stack, (a2, c2))
        end
    end
    return S .| .!outside
end

# Nearest-neighbor resize of a matrix to (mx, my)
resize_nearest(A, mx, my) = [A[clamp(ceil(Int, i * size(A, 1) / mx), 1, size(A, 1)),
                               clamp(ceil(Int, j * size(A, 2) / my), 1, size(A, 2))] for i in 1:mx, j in 1:my]

"""
    emblem_row(emblem, text_bitmap; height_ratio, gap, raise = 0)

Place an emblem (silhouette matrix) to the left of a text bitmap on one canvas. The emblem is
scaled to `height_ratio` × the text height and bottom-aligned with the text (shifted up by
`raise` × text height); `gap` is in units of the text height. Returns (canvas with both,
canvas with text only, emblem placement (i₀, j₀, w, h) in canvas pixels).
"""
function emblem_row(emblem, T; height_ratio = 1.5, gap = 0.25, raise = 0.0)
    tw, th = size(T)
    h = round(Int, height_ratio * th)
    w = round(Int, h * size(emblem, 1) / size(emblem, 2))
    g = round(Int, gap * th)
    ty = max(0, h - th) ÷ 2 + round(Int, -raise * th)      # text vertically centered on the emblem
    ty = clamp(ty, 0, max(h, th) - th)
    H = max(h, th)
    C = zeros(w + g + tw, H); Ct = zeros(size(C))
    C[1:w, 1:h] .= resize_nearest(emblem, w, h)
    C[w+g+1:end, ty+1:ty+th] .= max.(C[w+g+1:end, ty+1:ty+th], T)
    Ct[w+g+1:end, ty+1:ty+th] .= T
    return C, Ct, (0, 0, w, h)
end

# Center bitmap B, scaled by `scale` bitmap pixels per cell, in a grid; returns the mask and the
# cell offset (ox, oy) such that bitmap pixel p maps to cell ox + p / scale.
function place_bitmap(B, scale, grid)
    nx, ny = floor.(Int, size(B) ./ scale)
    M = falses(grid...)
    ox, oy = (grid[1] - nx) ÷ 2, (grid[2] - ny) ÷ 2
    for j in 1:ny, i in 1:nx
        M[ox + i, oy + j] = B[clamp(round(Int, (i - 0.5) * scale), 1, size(B, 1)),
                              clamp(round(Int, (j - 0.5) * scale), 1, size(B, 2))] > 0.5
    end
    return M, (ox, oy)
end

"""
    emblem_row_baseline(emblem, T; cap_px, desc_px, gap)

Like `emblem_row`, but the emblem is `cap_px` tall and sits on the text baseline, which is
`desc_px` above the bottom of the text bitmap (the descender depth). `gap` is in units of `cap_px`.
"""
function emblem_row_baseline(emblem, T; cap_px, desc_px, gap = 0.15)
    tw, th = size(T)
    h = round(Int, cap_px); w = round(Int, h * size(emblem, 1) / size(emblem, 2))
    g = round(Int, gap * cap_px)
    e₀ = round(Int, desc_px)                                   # emblem bottom = baseline
    H = max(th, e₀ + h)
    C = zeros(w + g + tw, H); Ct = zeros(size(C))
    C[1:w, e₀+1:e₀+h] .= resize_nearest(emblem, w, h)
    C[w+g+1:end, 1:th] .= max.(C[w+g+1:end, 1:th], T)
    Ct[w+g+1:end, 1:th] .= T
    return C, Ct, (0, e₀, w, h)
end

"""
    emblem_text_masks(; emblem_path, font, cells_per_stroke, margins, shadow = nothing, r = 1, text, tracking)

Land mask for an emblem followed by one line of text. The emblem artwork (cropped to its
silhouette) is as tall as the letters and bottom-aligned with them. `shadow = (sx, sy)` adds a
drop shadow under the text only, offset by (sx, sy) stroke widths, as part of the land.
`margins` are cells of open water on each side at the simulation resolution; `r` refines the
output grid for plotting (r = 1 gives the simulation mask). Returns a named tuple with the land
mask, the letters-only mask on the same grid, the emblem image and its box (x₀, y₀, width,
height) in output cells, and the grid size.
"""
function emblem_text_masks(; emblem_path, font, cells_per_stroke, margins, shadow = nothing, r = 1,
                           text = "OCEAN DYNAMICS", tracking = 0.02, gap = 0.12)
    emblem, img = load_emblem(emblem_path)
    T₀ = layout_bitmap([Line(; text, font, size = 200.0, tracking)])
    stroke = median_stroke(T₀)
    off = isnothing(shadow) ? (0, 0) : round.(Int, shadow .* stroke)
    Tl = add_shadow(T₀, off)                                  # letters + shadow (land)
    Tt = add_shadow(T₀, off; letters_only = true)             # letters only, same canvas
    lx, ly = max(0, -off[1]), max(0, -off[2])                 # letters' offset within that canvas

    h = size(T₀, 2); w = round(Int, h * size(emblem, 1) / size(emblem, 2)); g = round(Int, gap * h)
    H = max(size(Tl, 2), ly + h)
    C = zeros(w + g + size(Tl, 1), H); Ct = zeros(size(C))
    C[1:w, ly+1:ly+h] .= resize_nearest(emblem, w, h)        # emblem bottom = letters' bottom
    C[w+g+1:end, 1:size(Tl, 2)] .= max.(C[w+g+1:end, 1:size(Tl, 2)], Tl)
    Ct[w+g+1:end, 1:size(Tt, 2)] .= Tt

    scale = stroke / cells_per_stroke
    coarse = floor.(Int, size(C) ./ scale) .+ 2 .* margins
    grid = r .* coarse
    land, (ox, oy) = place_bitmap(C, scale / r, grid)
    letters, _ = place_bitmap(Ct, scale / r, grid)
    s = scale / r
    return (; land, letters, img, emblem_box = (ox, oy + ly / s, w / s, h / s), grid)
end

"""
    text_masks(; font, text, tracking, cells_per_stroke, margins, shadow = nothing, r = 1)

Land mask for one line of text alone (no emblem), with the same conventions as
`emblem_text_masks`. Returns (land, letters, img = nothing, emblem_box = nothing, grid).
"""
function text_masks(; font, text, tracking = 0.06, cells_per_stroke, margins, shadow = nothing, r = 1)
    T₀ = layout_bitmap([Line(; text, font, size = 200.0, tracking)])
    stroke = median_stroke(T₀)
    off = isnothing(shadow) ? (0, 0) : round.(Int, shadow .* stroke)
    Tl = add_shadow(T₀, off); Tt = add_shadow(T₀, off; letters_only = true)
    scale = stroke / cells_per_stroke
    coarse = floor.(Int, size(Tl) ./ scale) .+ 2 .* margins
    grid = r .* coarse
    land, _ = place_bitmap(Tl, scale / r, grid)
    letters, _ = place_bitmap(Tt, scale / r, grid)
    return (; land, letters, img = nothing, emblem_box = nothing, grid)
end
