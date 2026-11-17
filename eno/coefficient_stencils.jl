# The stencils and coefficients of the article's table "Coefficients c_{r,j} on a uniform grid", as a diagram for
# Wikimedia Commons, written to coefficient_stencils.svg and coefficient_stencils.png.
#
# Each row is a row of the table: a stencil of k = 1, 2 or 3 cells reaching r cells to the left of the cell I_i, drawn
# on the cells I_{i-2}, ..., I_{i+3}. The number in each cell of the stencil is the coefficient c_{r,j} of that cell's
# average in the value at the face x_{i+1/2}, the vertical line. Positive coefficients are shaded blue and negative ones
# red, more deeply for a larger magnitude; these are the two poles of the diverging pair in plotting.jl. The
# coefficients are computed exactly from Shu (1998), eq. (2.21), and checked against the article's table first.
#
# The text is set in New Computer Modern, the LaTeX typeface that Makie bundles, to match the article's formulas.
# Wikimedia's SVG renderer doesn't have it, so the text is drawn as outlines, as in the other CairoMakie figures.
#
#     julia --project coefficient_stencils.jl

include("eno.jl")
include("plotting.jl")

# The article's table, row by row: (k, r) => [c_{r,0}, ..., c_{r,k-1}].
const TABLE = [
        (1, -1) => [1 // 1], (1, 0) => [1 // 1],
        (2, -1) => [3 // 2, -1 // 2], (2, 0) => [1 // 2, 1 // 2], (2, 1) => [-1 // 2, 3 // 2],
        (3, -1) => [11 // 6, -7 // 6, 1 // 3], (3, 0) => [1 // 3, 5 // 6, -1 // 6],
        (3, 1) => [-1 // 6, 5 // 6, 1 // 3], (3, 2) => [1 // 3, -7 // 6, 11 // 6]
]

for ((k, r), row) in TABLE
    computed = [eno_coefficient(k, r, j) for j in 0:k-1]
    computed == row || error("k = $k, r = $r: eq. (2.21) gives $computed, but the table has $row")
end
println("The coefficients from eq. (2.21) match the article's table for k = 1, 2 and 3")

# The layout, in pixels from the top left.
const D = 80.0                          # width of a cell
const H = 44.0                          # height of a cell
const ROW = H + 10                      # from the top of one row to the top of the next
const GROUP_GAP = 24.0                  # extra space between the orders k
const LEFT = 280.0                      # the left edge of the cell I_{i-2}
const BRACKET = 98.0                    # the bracket that groups the rows of each order k
const R_COLUMN, S_COLUMN = 112.0, 208.0 # the left edges of the labels r = … and s = …
const OFFSETS = -2:3                    # the cells I_{i-2}, ..., I_{i+3}
const FACE = LEFT + 3D                  # the face x_{i+1/2}, between I_i and I_{i+1}
const FIRST_ROW = 92.0                  # the top of the first row
const WIDTH = LEFT + length(OFFSETS) * D + 28

"""`color` mixed with white: t = 0 is white and t = 1 is `color` itself."""
mix(color, t) = RGBf((1 + t * (Float32(c) - 1) for c in (color.r, color.g, color.b))...)

"""The fill of a cell whose coefficient is `c`: blue if positive, red if negative, deeper for a larger magnitude."""
shade(c) = mix(c > 0 ? BLUE : RED, 0.15 + 0.45 * min(abs(c) / 2, 1))

"""A coefficient as text, such as 11/6 or −7/6, with a true minus sign."""
function fraction(c)
    magnitude = denominator(c) == 1 ? string(abs(numerator(c))) : "$(abs(numerator(c)))/$(denominator(c))"
    return c < 0 ? "−" * magnitude : magnitude
end

"""
A label such as r = −1: the letter in italic, then the rest as text in the same typeface. Makie's LaTeX typesetting
leaves too little space around = and spaces a leading minus as if it were subtraction, which this avoids.
"""
equation(letter, n) = rich(rich(letter; font = :italic), " = $(n < 0 ? "−" : "")$(abs(n))")

cell_left(offset) = LEFT + (offset - first(OFFSETS)) * D

# The tops of the rows, with a gap between the orders.
tops = Dict{Tuple{Int, Int}, Float64}()
let y = FIRST_ROW
    for k in 1:3
        for r in -1:k-1
            tops[(k, r)] = y
            y += ROW
        end
        y += GROUP_GAP
    end
end
const LAST_BOTTOM = tops[(3, 2)] + H
const LEGEND_Y = LAST_BOTTOM + 40
const TOTAL_HEIGHT = LEGEND_Y + 30

# Makie's y axis points up, so y measured from the top becomes TOTAL_HEIGHT - y.
up(y) = TOTAL_HEIGHT - y
cell!(scene, x, top, w; kwargs...) = poly!(scene, Rect2f(x, up(top + H), w, H); kwargs...)

set_theme!(theme_latexfonts())

fig = Figure(size = (WIDTH, TOTAL_HEIGHT), figure_padding = 0, backgroundcolor = PAPER)
scene = fig.scene

for ((k, r), row) in TABLE
    top = tops[(k, r)]
    members = -r:(-r + k - 1)               # the offsets of I_{i-r}, ..., I_{i-r+k-1}
    for offset in OFFSETS
        offset in members && continue
        cell!(scene, cell_left(offset), top, D; color = PAPER, strokecolor = FAINT, strokewidth = 1.2)
    end
    for (j, offset) in enumerate(members)
        cell!(scene, cell_left(offset), top, D; color = shade(row[j]), strokecolor = MUTED, strokewidth = 1.2)
        text!(scene, cell_left(offset) + D / 2, up(top + H / 2); text = fraction(row[j]), fontsize = 23,
              align = (:center, :center), color = INK)
    end
    cell!(scene, cell_left(first(members)), top, k * D; color = :transparent, strokecolor = INK, strokewidth = 2.4)

    # r and s, in two columns with their equals signs aligned.
    text!(scene, R_COLUMN, up(top + H / 2); text = equation("r", r), fontsize = 23, align = (:left, :center),
          color = INK)
    text!(scene, S_COLUMN, up(top + H / 2); text = equation("s", k - r - 1), fontsize = 23,
          align = (:left, :center), color = INK)
end

# The order k of each group of rows, set off by a bracket.
for k in 1:3
    first_top, last_bottom = tops[(k, -1)], tops[(k, k - 1)] + H
    lines!(scene, [BRACKET + 6, BRACKET, BRACKET, BRACKET + 6], up.([first_top, first_top, last_bottom, last_bottom]);
           color = MUTED, linewidth = 1.5)
    text!(scene, 20, up((first_top + last_bottom) / 2); text = equation("k", k), fontsize = 25,
          align = (:left, :center), color = INK)
end

# The face, through every row, and the names of the cells, with the cell I_i in full ink.
lines!(scene, [FACE, FACE], [up(42), up(LAST_BOTTOM + 8)]; color = INK, linewidth = 2.5)
text!(scene, FACE, up(24); text = L"x_{i+1/2}", fontsize = 25, align = (:center, :center), color = INK)
for (offset, after) in zip(OFFSETS, ("-2", "-1", "", "+1", "+2", "+3"))
    text!(scene, cell_left(offset) + D / 2, up(66); text = latexstring("I_{i$after}"), fontsize = 24,
          align = (:center, :center), color = offset == 0 ? INK : MUTED)
end

# The key to the shading.
key_left = LEFT + length(OFFSETS) * D / 2 - 205
for (n, (c, label)) in enumerate(((1, "positive coefficient"), (-1, "negative coefficient")))
    x = key_left + (n - 1) * 222
    poly!(scene, Rect2f(x, up(LEGEND_Y + 9), 18, 18); color = shade(c), strokecolor = MUTED, strokewidth = 1.2)
    text!(scene, x + 28, up(LEGEND_Y); text = label, fontsize = 21, align = (:left, :center), color = MUTED)
end

save_figure(fig, "coefficient_stencils")
