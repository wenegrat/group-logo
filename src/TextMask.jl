"""
Rasterize a text string with Makie and sample it onto simulation grid nodes.
"""
module TextMask

using CairoMakie

export text_bitmap, text_mask, gaussian_smooth

"""
    text_bitmap(str; font=:bold, fontsize=200)

Render `str` (white on black) and return a cropped Float64 array `B[i, j]` in [0, 1],
with `i` increasing in x (left → right) and `j` increasing in y (bottom → top).
"""
function text_bitmap(str; font = :bold, fontsize = 200)
    W = round(Int, 1.2 * fontsize * length(str)) + 2fontsize
    H = 3fontsize
    fig = Figure(size = (W, H), backgroundcolor = :black, figure_padding = 0)
    ax = Axis(fig[1, 1], backgroundcolor = :black)
    hidedecorations!(ax); hidespines!(ax)
    limits!(ax, 0, W, 0, H)
    text!(ax, W / 2, H / 2; text = str, font, fontsize, color = :white,
          align = (:center, :center))
    img = Makie.colorbuffer(fig; px_per_unit = 1)       # rows top → bottom
    A = Float64.(Makie.Colors.Gray.(img))
    A = permutedims(reverse(A, dims = 1))              # A[x, y], y upward

    xs = findall(vec(maximum(A, dims = 2)) .> 0.5)
    ys = findall(vec(maximum(A, dims = 1)) .> 0.5)
    return A[first(xs):last(xs), first(ys):last(ys)]
end

"""
    text_mask(str, x, y; width, center=(0, 0), kw...)

Sample the text bitmap at nodes `x`, `y` (vectors), scaled so the text spans `width`
in x (aspect ratio preserved) and is centered at `center`. Returns `M[i, j]` in [0, 1].
"""
function text_mask(str, x, y; width, center = (0, 0), kw...)
    B = text_bitmap(str; kw...)
    nbx, nby = size(B)
    Δ = width / nbx                      # physical size of one bitmap pixel
    height = nby * Δ
    x₀ = center[1] - width / 2
    y₀ = center[2] - height / 2

    M = zeros(length(x), length(y))
    for (j, yⱼ) in enumerate(y), (i, xᵢ) in enumerate(x)
        p = (xᵢ - x₀) / Δ + 0.5
        q = (yⱼ - y₀) / Δ + 0.5
        M[i, j] = bilinear(B, p, q)
    end
    return M
end

function bilinear(B, p, q)
    nx, ny = size(B)
    (1 ≤ p ≤ nx && 1 ≤ q ≤ ny) || return 0.0
    i = min(floor(Int, p), nx - 1); a = p - i
    j = min(floor(Int, q), ny - 1); b = q - j
    i = max(i, 1); j = max(j, 1)
    return (1 - a) * (1 - b) * B[i, j]   + a * (1 - b) * B[i+1, j] +
           (1 - a) * b       * B[i, j+1] + a * b       * B[i+1, j+1]
end

"""
    gaussian_smooth(M, σ)

Separable Gaussian smoothing of `M` with standard deviation `σ` grid cells
(zero padding at the edges; the text is assumed to sit away from the boundaries).
"""
function gaussian_smooth(M, σ)
    σ > 0 || return copy(M)
    r = ceil(Int, 4σ)
    k = exp.(-((-r:r) .^ 2) ./ (2σ^2)); k ./= sum(k)
    conv1(A, dim) = begin
        out = zero(A)
        n = size(A, dim)
        for I in CartesianIndices(A)
            s = 0.0
            for (m, km) in zip(-r:r, k)
                l = I[dim] + m
                1 ≤ l ≤ n || continue
                J = Base.setindex(Tuple(I), l, dim)
                s += km * A[J...]
            end
            out[I] = s
        end
        out
    end
    return conv1(conv1(M, 1), 2)
end

end # module
