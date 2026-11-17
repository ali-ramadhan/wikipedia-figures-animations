# How ENO chooses a third-order stencil for the cell I_i on a non-uniform grid, drawn as a tree that follows the
# article's procedure. Start from I_i alone, with m = i the index of the leftmost cell; then, for l = 1 and 2, add
# I_{m-1} if |d^l_{m-1}| < |d^l_m|, and otherwise the cell on the right. Here d^l_j = U[x_{j-1/2}, ..., x_{j+l+1/2}] is
# the divided difference of the primitive over the faces of the l + 1 cells I_j, ..., I_{j+l} (Shu 1998, eqs. 2.32 and
# 2.35): d^0_j = ū_j, and d^l_j = (d^{l-1}_{j+1} - d^{l-1}_j) / (x_{j+l+1/2} - x_{j-1/2}) divides by the width the
# stencil spans. Each candidate stencil is labelled with its own d^l_j, so each step keeps the candidate whose label is
# smaller in magnitude. The three possible final stencils have r = i - m = 2, 1 and 0 cells to the left of I_i; r = 1
# can be reached in two ways. The cells' widths are arbitrary, chosen only to be visibly unequal. Written to
# stencil_schematic.svg and .png, in Latin Modern (see use_latin_modern! in plotting.jl), with the text drawn as
# outlines.
#
#     julia --project stencil_schematic.jl

include("plotting.jl")
use_latin_modern!()

# The layout, in pixels from the top left.
const WIDTHS = [62.0, 38.0, 50.0, 70.0, 44.0]        # the cells I_{i-2}, ..., I_{i+2}
const STRIP = sum(WIDTHS)
const CELL_HEIGHT = 38.0
const CENTER = 455.0                                 # the centre of the starting stencil
const SPREAD = (one = 155.0, two = 310.0)            # horizontal offsets of the stencils after steps 1 and 2
const TOPS = (start = 70.0, one = 290.0, two = 530.0)
const NOTES = 905.0                                  # the left edge of the notes on the right
const WIDTH, HEIGHT = 1250.0, 700.0
const STENCIL_FILL, STENCIL_EDGE, CELL_EDGE = colorant"#b7d3f6", COLORS.eno, GRAY

up(y) = HEIGHT - y                                   # Makie's y axis points up
lefts(cx) = cx - STRIP / 2 .+ [0; cumsum(WIDTHS)[1:end-1]]   # the left edges of the five cells

"""A row of five cells centred at `cx` with its top at `top`; `members` are the offsets (-2 to 2) in the stencil."""
function strip!(scene, cx, top, members)
    for (k, offset) in enumerate(-2:2)
        inside = offset in members
        poly!(scene, Rect2f(lefts(cx)[k], up(top + CELL_HEIGHT), WIDTHS[k], CELL_HEIGHT);
              color = inside ? STENCIL_FILL : PAPER, strokecolor = inside ? STENCIL_EDGE : CELL_EDGE, strokewidth = 1.5)
    end
    # The target cell I_i, outlined in every row.
    poly!(scene, Rect2f(lefts(cx)[3], up(top + CELL_HEIGHT), WIDTHS[3], CELL_HEIGHT); color = :transparent,
          strokecolor = INK, strokewidth = 3.5)
    return nothing
end

"""An arrow from `from` to `to`, both (x, y) in pixels from the top left."""
function arrow!(scene, from, to; color = MUTED)
    a, b = Point2f(from[1], up(from[2])), Point2f(to[1], up(to[2]))
    direction = (b - a) / sqrt(sum(abs2, b - a))
    normal = Point2f(-direction[2], direction[1])
    lines!(scene, [a, b - 9direction]; color, linewidth = 1.8)
    poly!(scene, [b, b - 12direction + 5normal, b - 12direction - 5normal]; color)
    return nothing
end

label!(scene, x, y, s; size = 22, align = (:center, :center), color = INK) =
    text!(scene, x, up(y); text = s, fontsize = size, align, color)

# Each stencil in the tree: the centre of its row, its top, the offsets of its cells, and its divided difference.
stencil(x, top, members, difference) = (; x, top, members, difference)

nodes = (
        start = stencil(CENTER, TOPS.start, (0,), L"d^0_i = \bar{u}_i"),
        left = stencil(CENTER - SPREAD.one, TOPS.one, (-1, 0),
                       L"d^1_{i-1} = \frac{d^0_i - d^0_{i-1}}{x_{i+1/2} - x_{i-3/2}}"),
        right = stencil(CENTER + SPREAD.one, TOPS.one, (0, 1),
                        L"d^1_i = \frac{d^0_{i+1} - d^0_i}{x_{i+3/2} - x_{i-1/2}}"),
        r2 = stencil(CENTER - SPREAD.two, TOPS.two, (-2, -1, 0),
                     L"d^2_{i-2} = \frac{d^1_{i-1} - d^1_{i-2}}{x_{i+1/2} - x_{i-5/2}}"),
        r1 = stencil(CENTER, TOPS.two, (-1, 0, 1), L"d^2_{i-1} = \frac{d^1_i - d^1_{i-1}}{x_{i+3/2} - x_{i-3/2}}"),
        r0 = stencil(CENTER + SPREAD.two, TOPS.two, (0, 1, 2),
                     L"d^2_i = \frac{d^1_{i+1} - d^1_i}{x_{i+5/2} - x_{i-1/2}}")
)

below(node) = node.top + CELL_HEIGHT + (node === nodes.start ? 26 : 40)   # where a stencil's difference is written

set_theme!(ARTICLE_THEME)

fig = Figure(size = (WIDTH, HEIGHT), figure_padding = 0, backgroundcolor = PAPER)
scene = fig.scene

# The names of the cells, above the starting stencil.
for (k, name) in enumerate((L"I_{i-2}", L"I_{i-1}", L"I_i", L"I_{i+1}", L"I_{i+2}"))
    label!(scene, lefts(CENTER)[k] + WIDTHS[k] / 2, TOPS.start - 22, name; size = 21)
end

for node in nodes
    strip!(scene, node.x, node.top, node.members)
    label!(scene, node.x, below(node), node.difference; size = 21)
end
for (node, r) in ((nodes.r2, 2), (nodes.r1, 1), (nodes.r0, 0))
    label!(scene, node.x, below(node) + 50, latexstring("r = $r"); size = 23)
end

# Each step leads to two candidates; the one with the smaller divided difference is kept.
for (parent, child, shift) in ((nodes.start, nodes.left, 0), (nodes.start, nodes.right, 0),
                               (nodes.left, nodes.r2, 0), (nodes.left, nodes.r1, -16),
                               (nodes.right, nodes.r1, 16), (nodes.right, nodes.r0, 0))
    arrow!(scene, (parent.x, below(parent) + (parent === nodes.start ? 20 : 36)), (child.x + shift, child.top - 10))
end

# The procedure, step by step, beside the rows. (MathTeXEngine spaces a colon as a relation, so there are none.)
label!(scene, NOTES, TOPS.start - 22, L"$d^l_j = U[x_{j-1/2}, \ldots, x_{j+l+1/2}]$"; align = (:left, :center))
label!(scene, NOTES, TOPS.start + CELL_HEIGHT / 2, L"start with $\{I_i\}$, so $m = i$"; align = (:left, :center))
for (l, top, bottom) in ((1, below(nodes.start), TOPS.one), (2, below(nodes.left) + 20, TOPS.two))
    middle = (top + bottom) / 2 + 6
    label!(scene, NOTES, middle - 45, latexstring("step \$l = $(l)\$"); align = (:left, :center))
    label!(scene, NOTES, middle - 15, latexstring("if \$|d^$(l)_{m-1}| < |d^$(l)_m|\$, add \$I_{m-1}\$");
           align = (:left, :center))
    label!(scene, NOTES, middle + 15, L"and replace $m$ by $m - 1$;"; align = (:left, :center))
    label!(scene, NOTES, middle + 45, latexstring("otherwise add \$I_{m+$(l)}\$"); align = (:left, :center))
end
label!(scene, NOTES, below(nodes.r1) + 50, L"$r = i - m$"; align = (:left, :center))

save_figure(fig, "stencil_schematic")
