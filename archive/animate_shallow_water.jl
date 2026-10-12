# Animate a field from a JLD2 output file.
# Usage: julia --project=. scripts/animate.jl output/<name>.jld2 <field> [clim]
using Oceananigans, CairoMakie, Printf

file, name = ARGS[1], ARGS[2]
clim = length(ARGS) ≥ 3 ? parse(Float64, ARGS[3]) : nothing

ts = FieldTimeSeries(file, name)
times = ts.times
x, y, _ = nodes(ts)
Lx, Ly = x[end] - x[1], y[end] - y[1]

n = Observable(1)
φ = @lift interior(ts[$n], :, :, 1)

symmetric = name in ("ζ", "ω")
colormap = symmetric ? :balance : :deep
colorrange = isnothing(clim) ? (symmetric ? (-1, 1) .* maximum(abs, ts[end]) : extrema(ts[1])) :
             (symmetric ? (-clim, clim) : (0, clim))

fig = Figure(size = (1600, round(Int, 1600 * Ly / Lx) + 60), backgroundcolor = :black)
ax = Axis(fig[1, 1], aspect = DataAspect(), backgroundcolor = :black)
hidedecorations!(ax); hidespines!(ax)
heatmap!(ax, x, y, φ; colormap, colorrange)
label = @lift @sprintf("t = %.0f / f", times[$n])
Label(fig[0, 1], label, color = :white, tellwidth = false)

out = replace(file, ".jld2" => "_$(name).mp4")
record(fig, out, eachindex(times); framerate = 24) do i
    n[] = i
end
@info "Saved $out"

# Snapshot strip of a few frames for quick inspection
idx = unique(round.(Int, range(1, length(times), length = 4)))
fig2 = Figure(size = (1600, 4 * round(Int, 1600 * Ly / Lx) + 100), backgroundcolor = :black)
for (r, i) in enumerate(idx)
    a = Axis(fig2[r, 1], aspect = DataAspect(), backgroundcolor = :black,
             title = @sprintf("t = %.0f / f", times[i]), titlecolor = :white)
    hidedecorations!(a); hidespines!(a)
    heatmap!(a, x, y, interior(ts[i], :, :, 1); colormap, colorrange)
end
png = replace(out, ".mp4" => "_snapshots.png")
save(png, fig2)
@info "Saved $png"
