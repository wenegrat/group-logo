# The logo layout used by the simulation and by the plotting scripts, defined once so both always
# build the same land mask.
#
# Layout: the Terrapin emblem followed by "OCEAN DYNAMICS" in Bubble Vibes, with a drop shadow
# under the text that is part of the land. Neither the emblem artwork nor the font is distributed
# with this repository; see README.md.

include(joinpath(@__DIR__, "Logo.jl"))

const REPO = normpath(joinpath(@__DIR__, ".."))

emblem_path() = get(ENV, "EMBLEM_PATH", joinpath(REPO, "assets", "emblem.png"))
logo_font_path() = get(ENV, "LOGO_FONT_PATH", joinpath(REPO, "fonts", "bubblevibes", "BubbleVibes-Regular.otf"))

function check_logo_inputs()
    isfile(emblem_path()) || error("Emblem image not found at $(emblem_path()); set EMBLEM_PATH (see README.md)")
    isfile(logo_font_path()) || error("Font not found at $(logo_font_path()); set LOGO_FONT_PATH (see README.md)")
end

const STROKE_PER_LD   = 2               # deformation radius Ld = stroke width / 2
const MARGINS_IN_LD   = (9.4, 7.7)      # open water east/west and north/south of the logo
const SHADOW_OFFSET   = (0.35, -0.35)   # drop shadow (in stroke widths): right and down

"""
    logo_layout(cells_per_Ld; r = 1)

Land mask and plotting data for the logo on a grid with `cells_per_Ld` cells per deformation
radius. `r > 1` refines the grid by that factor (for smooth plotting); `r = 1` is the simulation mask.
"""
function logo_layout(cells_per_Ld; r = 1)
    check_logo_inputs()
    return emblem_text_masks(; emblem_path = emblem_path(), font = FTFont(logo_font_path()),
                             cells_per_stroke = STROKE_PER_LD * cells_per_Ld,
                             margins = round.(Int, MARGINS_IN_LD .* cells_per_Ld),
                             shadow = SHADOW_OFFSET, r)
end
