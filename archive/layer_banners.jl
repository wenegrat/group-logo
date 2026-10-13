# Composite two transparent banner layers (front over back), then crop symmetrically about the
# center (where the text is) to the visible content. Plotting only.
#
# Usage: julia --project=. scripts/layer_banners.jl front.png back.png out.png
using CairoMakie
using CairoMakie.Colors: RGBA

front, back, out = ARGS
F = CairoMakie.FileIO.load(front); B = CairoMakie.FileIO.load(back)
over(f, b) = begin
    af, ab = Float64(f.alpha), Float64(b.alpha)
    a = af + ab * (1 - af)
    a == 0 && return RGBA(0.0, 0.0, 0.0, 0.0)
    mix(x, y) = (x * af + y * ab * (1 - af)) / a
    RGBA(mix(Float64(f.r), Float64(b.r)), mix(Float64(f.g), Float64(b.g)), mix(Float64(f.b), Float64(b.b)), a)
end
C = over.(F, B)
A = [c.alpha > 0.02 for c in C]
rows = findall(vec(any(A, dims = 2))); cols = findall(vec(any(A, dims = 1)))
pad = 8; cr, cc = size(C, 1) ÷ 2, size(C, 2) ÷ 2
hr = max(cr - first(rows), last(rows) - cr) + pad; hc = max(cc - first(cols), last(cols) - cc) + pad
C = C[max(cr - hr, 1):min(cr + hr, size(C, 1)), max(cc - hc, 1):min(cc + hc, size(C, 2))]
CairoMakie.FileIO.save(out, C)
println(out, " ", size(C, 2), " ", size(C, 1))
