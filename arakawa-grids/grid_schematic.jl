# The five Arakawa grids as an SVG diagram for Wikimedia Commons, written to arakawa_grids.svg.
#
# Each panel is a 2 × 2 block of square grid cells with the height h (or any scalar: pressure, temperature, ...) at the
# cell centers, as most textbooks and model documentation draw them, so that on grid C the velocity components sit on
# the cell faces normal to them. Grid E is grid B's block rotated by 45°: the same cells, so the same number of heights
# per unit area (Arakawa & Lamb 1977, p. 183).
#
# The diagram is the article's lead image, which Wikipedia shows about 300 px wide (MOS:IMGSIZE), so it is drawn to be
# read at that size and has no words: the caption says what the dots and arrows are. The only text is the panel letters
# and a key to h, u and v, kept as SVG text in fonts Wikimedia's renderer has, so that others can edit it.
#
#     julia grid_schematic.jl

output_directory = mkpath(get(ENV, "OUTPUT_DIRECTORY", joinpath(@__DIR__, "output")))

const D = 72.0                            # side of a grid cell, px
const BLOCK = 2                           # cells per side in each panel
const H = BLOCK * D / 2                   # half the side of a panel
const ARROW = 0.60D                       # arrow length, tail to tip
const SHAFT = 0.065D                      # width of an arrow's line
const HEAD = (length = 0.17D, width = 0.16D)
const DOT = 0.09D                         # radius of a height dot
const RING = 0.025D                       # white ring around each dot, so that on grid A the arrows pass behind it
const LINE = 0.03D                        # width of the cell edges
const TITLE = 0.50D                       # font size of the panel letters
const SYMBOL = 0.55D                      # font size of h, u and v in the key
const COLORS = (ink = "#0b0b0b", u = "#2a78d6", v = "#d4502a", line = "#a3a3a3", paper = "#ffffff")
const SANS = "font-family=\"DejaVu Sans, Arial, Helvetica, sans-serif\""
const SERIF = "font-family=\"DejaVu Serif, Georgia, serif\""

#####
##### The grids
#####

# The points of a panel, relative to its center, with y upward.
const CENTERS = [((i + 0.5) * D - H, (j + 0.5) * D - H) for i in 0:BLOCK-1 for j in 0:BLOCK-1]
const CORNERS = [(i * D - H, j * D - H) for i in 0:BLOCK for j in 0:BLOCK]
const X_FACES = [(i * D - H, (j + 0.5) * D - H) for i in 0:BLOCK for j in 0:BLOCK-1]   # faces normal to x
const Y_FACES = [((i + 0.5) * D - H, j * D - H) for i in 0:BLOCK-1 for j in 0:BLOCK]   # faces normal to y

const GRIDS = (
        A = (u = CENTERS, v = CENTERS, angle = 0.0),     # u, v and h together
        B = (u = CORNERS, v = CORNERS, angle = 0.0),     # u and v together at the corners
        C = (u = X_FACES, v = Y_FACES, angle = 0.0),     # u on the faces normal to x, v on those normal to y
        D = (u = Y_FACES, v = X_FACES, angle = 0.0),     # as C with u and v interchanged
        E = (u = CORNERS, v = CORNERS, angle = π / 4)    # grid B rotated by 45°
)

#####
##### SVG elements
#####

svg = IOBuffer()
emit(parts...) = print(svg, parts..., "\n")
fmt(x) = string(round(x; digits = 2))   # numbers to 0.01 px

"""An arrow of length ARROW centered on (x, y), pointing along (dx, dy) on the page (y downward)."""
function arrow(x, y, dx, dy, color)
    along(a) = (x + dx * a, y + dy * a)
    across(p, b) = (p[1] - dy * b, p[2] + dx * b)
    tip, tail = along(ARROW / 2), along(-ARROW / 2)
    shaft_end = along(ARROW / 2 - 0.6HEAD.length)   # inside the head
    back, notch = along(ARROW / 2 - HEAD.length), along(ARROW / 2 - 0.72HEAD.length)

    emit("<line x1=\"$(fmt(tail[1]))\" y1=\"$(fmt(tail[2]))\" ",
         "x2=\"$(fmt(shaft_end[1]))\" y2=\"$(fmt(shaft_end[2]))\" ",
         "stroke=\"$color\" stroke-width=\"$(fmt(SHAFT))\"/>")

    head = (tip, across(back, HEAD.width / 2), notch, across(back, -HEAD.width / 2))
    emit("<path d=\"M ", join(("$(fmt(p[1])) $(fmt(p[2]))" for p in head), " L "), " Z\" fill=\"$color\"/>")

    return nothing
end

u_arrow(x, y) = arrow(x, y, 1, 0, COLORS.u)
v_arrow(x, y) = arrow(x, y, 0, -1, COLORS.v)

h_point(x, y) = emit("<circle cx=\"$(fmt(x))\" cy=\"$(fmt(y))\" r=\"$(fmt(DOT + RING / 2))\" fill=\"$(COLORS.ink)\" ",
                     "stroke=\"$(COLORS.paper)\" stroke-width=\"$(fmt(RING))\"/>")

edge(p, q) = emit("<line x1=\"$(fmt(p[1]))\" y1=\"$(fmt(p[2]))\" x2=\"$(fmt(q[1]))\" y2=\"$(fmt(q[2]))\" ",
                  "stroke=\"$(COLORS.line)\" stroke-width=\"$(fmt(LINE))\" stroke-linecap=\"square\"/>")

text(x, y, s; size, style) =
    emit("<text x=\"$(fmt(x))\" y=\"$(fmt(y))\" $style font-size=\"$(fmt(size))\" fill=\"$(COLORS.ink)\">$s</text>")

"""A grid's block of cells, with its center at (cx, cy) on the page."""
function panel(cx, cy, grid)
    c, s = cos(grid.angle), sin(grid.angle)
    page((x, y)) = (cx + c * x - s * y, cy - (s * x + c * y))

    for q in 0:BLOCK   # the cell edges
        edge(page((q * D - H, -H)), page((q * D - H, H)))
        edge(page((-H, q * D - H)), page((H, q * D - H)))
    end

    foreach(p -> u_arrow(page(p)...), grid.u)
    foreach(p -> v_arrow(page(p)...), grid.v)
    foreach(p -> h_point(page(p)...), CENTERS)

    return nothing
end

"""A key to the three symbols, centered on (cx, cy)."""
function key(cx, cy)
    for (n, (draw, name)) in enumerate(((h_point, "h"), (u_arrow, "u"), (v_arrow, "v")))
        y = cy + (n - 2) * 0.85D
        draw(cx - 0.30D, y)
        text(cx + 0.25D, y + 0.27SYMBOL, name; size = SYMBOL, style = "$SERIF font-style=\"italic\"")
    end

    return nothing
end

#####
##### Layout
#####

# A, B and C in the first row; D, E and the key in the second. The blocks are evenly spaced, the arrows on their edges
# reach beyond them, and E's diamond hangs from the same line as D's square.
const MARGIN, GAP = 0.12D, 0.30D
const KEY_HALF = 0.60D

# How far each panel reaches from its center, arrows included.
reach = (
        A = H + HEAD.width / 2,
        B = H + ARROW / 2,
        C = H + ARROW / 2,
        D = H + HEAD.width / 2,
        E = sqrt(2) * H + ARROW / 2
)

pitch = max(reach.D + reach.E, reach.B + reach.C, reach.E + KEY_HALF) + GAP   # between the centers of the columns
x1 = MARGIN + max(reach.A, reach.D)
x2, x3 = x1 + pitch, x1 + 2pitch
cap = 0.73TITLE                                    # height of a capital letter
title1 = MARGIN + cap
top1 = title1 + 0.22D + ARROW / 2                  # the top edge of the first row's blocks
title2 = top1 + 2H + ARROW / 2 + 0.38D + cap
top2 = title2 + 0.22D + ARROW / 2
width = x3 + max(reach.C, KEY_HALF) + MARGIN
height = top2 + 2sqrt(2) * H + ARROW / 2 + MARGIN

centers = (A = (x1, top1 + H), B = (x2, top1 + H), C = (x3, top1 + H), D = (x1, top2 + H), E = (x2, top2 + sqrt(2) * H))

#####
##### The diagram
#####

emit("<?xml version=\"1.0\" encoding=\"UTF-8\"?>")
emit("<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"$(fmt(width))\" height=\"$(fmt(height))\" ",
     "viewBox=\"0 0 $(fmt(width)) $(fmt(height))\">")
emit("<title>The five Arakawa grids, A to E</title>")
emit("<rect width=\"100%\" height=\"100%\" fill=\"$(COLORS.paper)\"/>")

for name in keys(GRIDS)
    x, y = centers[name]
    title_baseline = name in (:A, :B, :C) ? title1 : title2
    text(x, title_baseline, name; size = TITLE, style = "$SANS font-weight=\"bold\" text-anchor=\"middle\"")
    panel(x, y, GRIDS[name])
end

key(x3, centers.E[2])
emit("</svg>")

path = joinpath(output_directory, "arakawa_grids.svg")
write(path, take!(svg))
println("wrote $path ($(round(Int, width)) × $(round(Int, height)) px)")
