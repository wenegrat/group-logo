include(joinpath(@__DIR__, "..", "src", "TextMask.jl"))
using .TextMask, CairoMakie
B = text_bitmap("UMD Ocean Dynamics")
@show size(B) size(B, 1) / size(B, 2)
x = range(-64, 64, length = 768); y = range(-16, 16, length = 192)
M = gaussian_smooth(text_mask("UMD Ocean Dynamics", x, y; width = 110), 1.5)
fig = Figure(size = (1200, 360)); ax = Axis(fig[1, 1], aspect = DataAspect())
heatmap!(ax, x, y, M); mkpath("output"); save("output/check_mask.png", fig)
