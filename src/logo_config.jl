# The logo layouts used by the simulation and by the plotting scripts, defined once so both always
# build the same land mask.
#
#   :text    "OCEAN DYNAMICS" in Montserrat ExtraBold (SIL Open Font License, included)   ← current
#   :emblem  the Terrapin emblem + "OCEAN DYNAMICS" in Bubble Vibes with a land drop shadow
#            (earlier runs; neither the emblem nor the font is distributed, see README.md)

include(joinpath(@__DIR__, "Logo.jl"))

const REPO = normpath(joinpath(@__DIR__, ".."))

emblem_path() = get(ENV, "EMBLEM_PATH", joinpath(REPO, "assets", "emblem.png"))
bubble_font_path() = get(ENV, "LOGO_FONT_PATH", joinpath(REPO, "fonts", "bubblevibes", "BubbleVibes-Regular.otf"))
text_font_path() = joinpath(REPO, "fonts", "montserrat", "Montserrat-ExtraBold.ttf")

const STROKE_PER_LD   = 2               # deformation radius Ld = stroke width / 2
const MARGINS_IN_LD   = (9.4, 7.7)      # :emblem layout: open water east/west and north/south of the logo
const TEXT_MARGINS_IN_LD = (17.5, 19.5) # :text layout: ≈ 3:1 desktop banner, text ≈ 75% of the width
const SHADOW_OFFSET   = (0.35, -0.35)   # :emblem layout only: drop shadow (stroke widths), right and down

"""
    logo_layout(cells_per_Ld; r = 1, layout = :text, margins_in_Ld)

Land mask and plotting data for a logo layout on a grid with `cells_per_Ld` cells per deformation
radius. `r > 1` refines the grid by that factor (for smooth plotting); `r = 1` is the simulation mask.
"""
function logo_layout(cells_per_Ld; r = 1, layout = :text,
                     margins_in_Ld = layout == :text ? TEXT_MARGINS_IN_LD : MARGINS_IN_LD)
    cells_per_stroke = STROKE_PER_LD * cells_per_Ld
    margins = round.(Int, margins_in_Ld .* cells_per_Ld)
    if layout == :text
        return text_masks(; font = FTFont(text_font_path()), text = "OCEAN DYNAMICS", tracking = 0.06,
                          cells_per_stroke, margins, r)
    elseif layout == :emblem
        for p in (emblem_path(), bubble_font_path())
            isfile(p) || error("$p not found; set EMBLEM_PATH / LOGO_FONT_PATH (see README.md)")
        end
        return emblem_text_masks(; emblem_path = emblem_path(), font = FTFont(bubble_font_path()),
                                 cells_per_stroke, margins, shadow = SHADOW_OFFSET, r)
    else
        error("unknown layout $layout")
    end
end

# Layout of a run from its name, e.g. "heated_text_985x169x8" → :text
layout_from_name(name) = occursin("emblem", name) ? :emblem : occursin("_text_", name) ? :text : nothing
