# The notation of the article's worked examples of reconstruction on a non-uniform grid (k = 2 and k = 3, in
# "Reconstruction on a fixed stencil"), written to reconstruction_notation.svg and .png. Four cells I_{i-1}, ...,
# I_{i+2} of unequal widths Δ_j, with their averages ū_j; their faces on the x axis above and on the ξ axis below, where
# ξ = x - x_{i+1/2} puts the origin at the face x_{i+1/2}; the two values u^± at that face; brackets for the stencils of
# the two examples, with the cells of the larger one shaded; and marks on the ξ axis at the cell centres that the
# second-order example uses. The cells' widths are arbitrary, chosen only to be visibly unequal. In Latin Modern (see
# use_latin_modern! in plotting.jl), with the text drawn as outlines.
#
#     julia --project reconstruction_notation.jl

include("plotting.jl")
use_latin_modern!()

# The layout, in pixels from the top left.
const WIDTHS = [175.0, 205.0, 245.0, 150.0]          # I_{i-1}, I_i, I_{i+1}, I_{i+2}
const LEFT = 120.0                                   # the face x_{i-3/2}
const FACES = LEFT .+ [0; cumsum(WIDTHS)]            # x_{i-3/2}, ..., x_{i+5/2}
const ORIGIN = FACES[3]                              # x_{i+1/2}, where ξ = 0
const BRACKETS = (k3 = 42.0, k2 = 92.0)              # the heights of the two stencil brackets
const X_AXIS, BOX_TOP, BOX_BOTTOM = 160.0, 160.0, 252.0
const ARROWS, XI_AXIS = 288.0, 330.0
const WIDTH, HEIGHT = 1010.0, 425.0
const STENCIL_FILL, CELL_EDGE = colorant"#e3eefb", GRAY

up(y) = HEIGHT - y
center(k) = (FACES[k] + FACES[k + 1]) / 2

label!(scene, x, y, s; size = 22, align = (:center, :center), color = INK) =
    text!(scene, x, up(y); text = s, fontsize = size, align, color)
segment!(scene, x1, y1, x2, y2; color = INK, linewidth = 1.5, linestyle = :solid) =
    lines!(scene, [x1, x2], [up(y1), up(y2)]; color, linewidth, linestyle)

"""An arrowhead at (x, y) pointing in the direction `dir` = ±1 along the horizontal."""
head!(scene, x, y, dir; color = INK) =
    poly!(scene, Point2f[(x, up(y)), (x - 11dir, up(y) + 4.5), (x - 11dir, up(y) - 4.5)]; color)

"""A bracket over the cells `first`:`last`, at height `y`, labelled `name` at its left."""
function bracket!(scene, first, last, y, name)
    x1, x2 = FACES[first] + 3, FACES[last + 1] - 3
    lines!(scene, [x1, x1, x2, x2], up.([y + 12, y, y, y + 12]); color = COLORS.eno, linewidth = 2)
    label!(scene, x1 - 12, y + 4, name; align = (:right, :center), color = COLORS.eno)
    return nothing
end

set_theme!(ARTICLE_THEME)

fig = Figure(size = (WIDTH, HEIGHT), figure_padding = 0, backgroundcolor = PAPER)
scene = fig.scene

# The stencils of the two examples.
bracket!(scene, 1, 3, BRACKETS.k3, L"k = 3")
bracket!(scene, 2, 3, BRACKETS.k2, L"k = 2")

# The cells, with their names and averages; the stencil cells of the examples are shaded lightly.
for k in 1:4
    poly!(scene, Rect2f(FACES[k], up(BOX_BOTTOM), WIDTHS[k], BOX_BOTTOM - BOX_TOP);
          color = k <= 3 ? STENCIL_FILL : PAPER, strokecolor = CELL_EDGE, strokewidth = 1.5)
end
for (k, (name, average)) in enumerate(((L"I_{i-1}", L"\bar{u}_{i-1}"), (L"I_i", L"\bar{u}_i"),
                                       (L"I_{i+1}", L"\bar{u}_{i+1}"), (L"I_{i+2}", L"\bar{u}_{i+2}")))
    label!(scene, center(k), BOX_TOP + 26, name; size = 23)
    label!(scene, center(k), BOX_BOTTOM - 24, average; size = 23, color = MUTED)
end
# The two values at the face x_{i+1/2}, one from each side.
label!(scene, ORIGIN - 10, (BOX_TOP + BOX_BOTTOM) / 2, L"u^-_{i+1/2}"; align = (:right, :center), size = 21)
label!(scene, ORIGIN + 10, (BOX_TOP + BOX_BOTTOM) / 2, L"u^+_{i+1/2}"; align = (:left, :center), size = 21)

# The x axis along the top of the cells, with the faces.
segment!(scene, FACES[1] - 40, X_AXIS, FACES[end] + 55, X_AXIS)
head!(scene, FACES[end] + 55, X_AXIS, 1)
label!(scene, FACES[end] + 70, X_AXIS, L"x"; align = (:left, :center), size = 24)
for (k, name) in enumerate((L"x_{i-3/2}", L"x_{i-1/2}", L"x_{i+1/2}", L"x_{i+3/2}", L"x_{i+5/2}"))
    label!(scene, FACES[k], X_AXIS - 16, name; size = 21, align = (:center, :bottom))
end

# The widths.
for (k, name) in enumerate((L"\Delta_{i-1}", L"\Delta_i", L"\Delta_{i+1}", L"\Delta_{i+2}"))
    x1, x2 = FACES[k] + 3, FACES[k + 1] - 3
    segment!(scene, x1 + 9, ARROWS, x2 - 9, ARROWS; color = MUTED)
    head!(scene, x1, ARROWS, -1; color = MUTED)
    head!(scene, x2, ARROWS, 1; color = MUTED)
    label!(scene, center(k), ARROWS - 14, name; size = 21, color = MUTED)
end

# The faces carried down to the ξ axis, and the face x_{i+1/2} drawn strongly: ξ = 0 there.
for (k, face) in enumerate(FACES)
    k == 3 && continue
    segment!(scene, face, BOX_BOTTOM, face, XI_AXIS; color = CELL_EDGE, linewidth = 1, linestyle = :dash)
end
segment!(scene, ORIGIN, X_AXIS - 6, ORIGIN, XI_AXIS + 6; color = INK, linewidth = 2.5)

# The ξ axis, with the faces and, below them, the centres of I_i and I_{i+1}.
segment!(scene, FACES[1] - 40, XI_AXIS, FACES[end] + 55, XI_AXIS)
head!(scene, FACES[end] + 55, XI_AXIS, 1)
label!(scene, FACES[end] + 70, XI_AXIS, L"\xi"; align = (:left, :center), size = 24)
xi_names = (L"-\Delta_{i-1} - \Delta_i", L"-\Delta_i", L"0", L"\Delta_{i+1}", L"\Delta_{i+1} + \Delta_{i+2}")
for (k, name) in enumerate(xi_names)
    segment!(scene, FACES[k], XI_AXIS - 6, FACES[k], XI_AXIS + 6)
    label!(scene, FACES[k], XI_AXIS + 14, name; size = 21, align = (:center, :top))
end
for (k, name) in ((2, L"-\frac{1}{2}\Delta_i"), (3, L"\frac{1}{2}\Delta_{i+1}"))
    scatter!(scene, [center(k)], [up(XI_AXIS)]; color = COLORS.eno, markersize = 9)
    label!(scene, center(k), XI_AXIS + 50, name; size = 20, align = (:center, :top), color = COLORS.eno)
end

save_figure(fig, "reconstruction_notation")
