using CairoMakie

"""
    draw_letters!(ax, xc, yc, land; outline = 2)

Overlay the land mask as white letters with a black outline.
"""
function draw_letters!(ax, xc, yc, land; outline = 2)
    heatmap!(ax, xc, yc, ifelse.(land, 1.0, NaN); colormap = [:white], colorrange = (0, 1))
    contour!(ax, xc, yc, Float64.(land); levels = [0.5], color = :black, linewidth = outline)
end

"""
    enclosed_water(land; periodic = (true, false))

Ocean cells not connected to the largest body of water (4-connectivity), i.e. water trapped
inside letters. `periodic` sets whether connectivity wraps in x and y.
"""
function enclosed_water(land; periodic = (true, false))
    Nx, Ny = size(land)
    label = zeros(Int, Nx, Ny)
    sizes = Int[]
    for j in 1:Ny, i in 1:Nx
        (land[i, j] || label[i, j] > 0) && continue
        push!(sizes, 0)
        id = length(sizes)
        label[i, j] = id
        stack = [(i, j)]
        while !isempty(stack)
            a, c = pop!(stack)
            sizes[id] += 1
            for (da, dc) in ((1, 0), (-1, 0), (0, 1), (0, -1))
                a2, c2 = a + da, c + dc
                a2 = periodic[1] ? mod1(a2, Nx) : a2
                c2 = periodic[2] ? mod1(c2, Ny) : c2
                (1 ≤ a2 ≤ Nx && 1 ≤ c2 ≤ Ny) || continue
                (land[a2, c2] || label[a2, c2] > 0) && continue
                label[a2, c2] = id
                push!(stack, (a2, c2))
            end
        end
    end
    return (label .> 0) .& (label .!= argmax(sizes))
end
