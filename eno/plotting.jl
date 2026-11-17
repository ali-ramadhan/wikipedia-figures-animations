# Colors, theme and text shared by the figures and animations, and the functions that save them to the output
# directory.

using CairoMakie
using CairoMakie.Makie.LaTeXStrings: latexstring   # a LaTeX label built from parts, as L"…" can't be
using Printf

include("video_encoding.jl")

output_directory = mkpath(get(ENV, "OUTPUT_DIRECTORY", joinpath(@__DIR__, "output")))
output_file(name) = joinpath(output_directory, name)

#####
##### Colors
#####

const INK = colorant"#0b0b0b"
const MUTED = colorant"#52514e"
const FAINT = colorant"#c3c2b7"
const PAPER = colorant"#ffffff"

# The diverging pair of the dataviz reference palette (light mode), and its gray midpoint.
const RED = colorant"#e34948"
const GRAY = colorant"#a9a8a0"
const BLUE = colorant"#2a78d6"

# Slots of the dataviz reference palette (light mode). Each scheme keeps its color in every figure.
const COLORS = (
        exact = INK,
        eno = colorant"#2a78d6",       # blue
        fixed = colorant"#eb6834",     # orange
        tvd = colorant"#1baf7a",       # aqua
        upwind = colorant"#6250d6",    # violet
        # Raspberry, from outside the palette: its free slots are too close to the fixed stencils' orange (magenta,
        # yellow and red, ΔE 7–14) or confusable with it in protan vision (green, ΔE 3.2). This one is at least ΔE 20
        # from every color in the advection figure in normal vision and 14.5 in simulated color-vision deficiency.
        centered = colorant"#b02a7f",
        average = MUTED
)

# The shift r of a third-order stencil, for flow to the right: fully upwind (r = 2), one cell either side (r = 1) and
# downwind-biased (r = 0), as the two poles of the diverging pair and its gray midpoint.
const SHIFT_COLORS = (RED, GRAY, BLUE)   # r = 0, 1, 2

#####
##### Theme and text
#####

const ARTICLE_THEME = Theme(
        fontsize = 22,
        fonts = (; regular = "DejaVu Sans", bold = "DejaVu Sans Bold", italic = "DejaVu Sans Oblique"),
        backgroundcolor = PAPER,
        Axis = (;
                xgridvisible = false,
                ygridvisible = false,
                topspinevisible = false,
                rightspinevisible = false,
                xticklabelcolor = MUTED,
                yticklabelcolor = MUTED,
                titlefont = :bold,
                titlealign = :left
        ),
        Legend = (; framevisible = false, labelcolor = INK)
)

# The panels of the animations have a full frame, with ticks on all four sides.
const FULL_FRAME = (; topspinevisible = true, rightspinevisible = true, xticksmirrored = true, yticksmirrored = true)

italic(s) = rich(s; font = :italic)

# Numbers as tick labels, with true minus signs and no trailing zeros.
tick_label(value) = replace(string(value), r"\.0$" => "", "-" => "−")
labelled_ticks(values) = (values, tick_label.(values))

"""
    use_latin_modern!()

Typeset LaTeX strings in Latin Modern at its 8-point design size, from fonts/ (GUST Font License). LaTeX uses its
smaller design sizes for small text because they are wider and sturdier; Makie bundles only the 10-point New Computer
Modern, which looks narrow in a figure. MathTeXEngine takes letters from the italic font and words, digits and
punctuation from the regular one; symbols such as = and − still come from New Computer Modern Math. The text is drawn
as outlines, since Wikimedia's SVG renderer doesn't have these fonts.
"""
function use_latin_modern!()
    MTE = CairoMakie.Makie.MathTeXEngine
    latin_modern(name) = joinpath(@__DIR__, "fonts", name)
    new_cm(name) = joinpath(pkgdir(MTE), "assets", "fonts", "NewComputerModern", name)
    MTE.set_texfont_family!(regular = latin_modern("lmroman8-regular.otf"),
                            italic = latin_modern("lmroman8-italic.otf"),
                            bold = new_cm("NewCM10-Bold.otf"),
                            bolditalic = new_cm("NewCM10-BoldItalic.otf"),
                            math = new_cm("NewCMMath-Regular.otf"))
    return nothing
end

#####
##### Saving figures and animations
#####

"""Save `figure` to the output directory as `<name>.<format>` for each of `formats`."""
function save_figure(figure, name; formats = ("svg", "png"), px_per_unit = 2)
    for format in formats
        path = output_file("$name.$format")
        save(path, figure; px_per_unit)
        println("Wrote $path")
    end

    return nothing
end

"""
    save_animation(show_frame!, figure, frames, name; framerate, hold_seconds = 0)

Call `show_frame!(frame)` for each of `frames`, then hold the last one for `hold_seconds`, to write the lossless master
`<name>.mkv` and the video `<name>.webm` to the output directory, at twice the figure's size. With STILL_ONLY=1 it does
nothing, so that an animation script writes only its still.
"""
function save_animation(show_frame!, figure, frames, name; framerate, hold_seconds = 0)
    get(ENV, "STILL_ONLY", "0") == "1" && return nothing

    shown = [frames; fill(last(frames), round(Int, hold_seconds * framerate))]
    master = output_file("$name.mkv")
    seconds = @elapsed record_lossless(show_frame!, figure, master, shown; framerate, px_per_unit = 2)
    @printf("Recorded %d frames in %.1f s\n", length(shown), seconds)
    encode_webm(master, output_file("$name.webm"))

    return nothing
end
