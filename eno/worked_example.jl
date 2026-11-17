# A worked example of the stencil choice, animated: third-order ENO choosing the stencils of the two cells next to a
# jump, from the cell averages 1.0, 1.2, 1.3, 0.0, 0.1 and 0.3 of cells 1 to 6. Each step shows the two differences
# being compared, and the result is set against the fixed stencil centred on the cell. Writes worked_example.webm, its
# lossless master worked_example.mkv, and a still in three panels, worked_example.svg and worked_example.png.
# STILL_ONLY=1 writes only the still.
#
#     julia --project worked_example.jl

include("eno.jl")
include("plotting.jl")

const AVERAGES = [1.0, 1.2, 1.3, 0.0, 0.1, 0.3]
const k = 3
const DURATION = 44.0
const FRAMERATE = 30

average(j) = AVERAGES[j]

#####
##### The timeline, in seconds
#####

const CAPTIONS = [
        (0, "Six cell averages, with a jump between cells 3 and 4"),
        (3, "Reconstructing cell 3: the stencil starts with the cell itself"),
        (5, "Step 1: compare the first differences on each side"),
        (8, "The left one is smaller in magnitude, so cell 2 joins the stencil"),
        (10, "Step 2: compare the second differences on each side"),
        (13, "The left one is smaller again, so cell 1 joins: the stencil avoids the jump"),
        (15, "The parabola with these three averages gives the values at cell 3's faces"),
        (19, "A fixed stencil, centered on cell 3, crosses the jump and overshoots"),
        (23, "Cell 4: the first differences are −1.3 and 0.1, so cell 5 joins"),
        (27, "The second differences are 1.4 and 0.1, so cell 6 joins"),
        (30, "ENO's parabola for cell 4"),
        (33, "The fixed stencil for cell 4 crosses the jump too"),
        (36, "Each cell takes its information from its own side of the jump, which stays sharp")
]

caption(t) = last(CAPTIONS[findlast(c -> first(c) ≤ t, CAPTIONS)])

# Bands behind the cells, marking the target cell (more strongly) and the stencil as it grows.
const BANDS = [
        (cells = 3:3, appear = 3, disappear = 23, strength = 0.22),
        (cells = 2:3, appear = 8, disappear = 13, strength = 0.14),
        (cells = 1:3, appear = 13, disappear = 23, strength = 0.14),
        (cells = 4:4, appear = 23, disappear = 33, strength = 0.22),
        (cells = 4:5, appear = 25.5, disappear = 28, strength = 0.14),
        (cells = 4:6, appear = 28, disappear = 33, strength = 0.14)
]

# The two candidate stencils at each step, each in its own row of brackets and labelled with the difference that its
# new cell brings in. Each is decided at `decided`, and is then blue if chosen and faint if not.
const CANDIDATES = [
        (cells = 2:3, row = 1, difference = "0.1", appear = 5, decided = 8, disappear = 10, chosen = true),
        (cells = 3:4, row = 2, difference = "−1.3", appear = 5, decided = 8, disappear = 10, chosen = false),
        (cells = 1:3, row = 1, difference = "−0.1", appear = 10, decided = 13, disappear = 15, chosen = true),
        (cells = 2:4, row = 2, difference = "−1.4", appear = 10, decided = 13, disappear = 15, chosen = false),
        (cells = 3:4, row = 1, difference = "−1.3", appear = 23.5, decided = 25.5, disappear = 27, chosen = false),
        (cells = 4:5, row = 2, difference = "0.1", appear = 23.5, decided = 25.5, disappear = 27, chosen = true),
        (cells = 3:5, row = 1, difference = "1.4", appear = 27, decided = 28, disappear = 30, chosen = false),
        (cells = 4:6, row = 2, difference = "0.1", appear = 27, decided = 28, disappear = 30, chosen = true)
]

# The parabolas from ENO's stencils and from the fixed stencils centered on the target cell. Each part across the target
# cell stays until `disappear_cell`, and its two face values are labelled at the sides given by `labels`.
const PARABOLAS = [
        (cell = 3, r = 2, color = COLORS.eno, dashed = false, appear = 15, disappear_rest = 23, disappear_cell = Inf,
         labels = ((:right, :below), (:right, :above))),
        (cell = 3, r = 1, color = COLORS.fixed, dashed = true, appear = 19, disappear_rest = 23, disappear_cell = 23,
         labels = ((:left, :above), (:right, :below))),
        (cell = 4, r = 0, color = COLORS.eno, dashed = false, appear = 30, disappear_rest = 36, disappear_cell = Inf,
         labels = ((:left, :below), (:right, :below))),
        (cell = 4, r = 1, color = COLORS.fixed, dashed = true, appear = 33, disappear_rest = 36, disappear_cell = 36,
         labels = ((:left, :above), (:right, :below)))
]

#####
##### Drawing
#####

"""How visible an element is at time t: it fades in at `appear` and out at `disappear`, over `fade` seconds."""
visibility(t, appear, disappear = Inf; fade = 0.5) = clamp(min((t - appear) / fade, (disappear - t) / fade), 0, 1)

visible(clock, appear, disappear = Inf) = @lift visibility($clock, appear, disappear)

function with_alpha(color, α)
    c = Makie.to_color(color)
    return RGBAf(c.r, c.g, c.b, α)
end

faded(color, α) = @lift with_alpha(color, $α)

function band!(ax, clock; cells, appear, disappear, strength)
    α = @lift strength * visibility($clock, appear, disappear)
    return poly!(ax, Rect(first(cells) - 0.5, -0.1, length(cells), 1.65); color = faded(COLORS.eno, α))
end

"""A candidate stencil: a bracket under its cells, labelled with its difference, gray until it is decided."""
function candidate!(ax, clock; cells, row, difference, appear, decided, disappear, chosen, fontsize)
    y = row == 1 ? -0.3 : -0.58
    x₀, x₁ = first(cells) - 0.45, last(cells) + 0.45

    color = @lift begin
        base = $clock < decided ? MUTED : chosen ? COLORS.eno : FAINT
        with_alpha(base, visibility($clock, appear, disappear))
    end
    linewidth = @lift $clock ≥ decided && chosen ? 4.5 : 2.5

    lines!(ax, [x₀, x₀, x₁, x₁], [y + 0.08, y, y, y + 0.08]; color, linewidth)
    text!(ax, (x₀ + x₁) / 2, y - 0.05; text = difference, align = (:center, :top), fontsize, color, font = :bold)

    return nothing
end

"""The parabola from the stencil of `cell` with shift `r`: thin across its other cells and thick across `cell`."""
function parabola!(ax, clock; cell, r, color, dashed, appear, disappear_rest, disappear_cell, labels, fontsize)
    a = reconstruction_polynomial(average, cell, k, r)

    for j in cell-r:cell-r+k-1
        ξ = range(j - cell, j - cell + 1; length = 40)
        own = j == cell
        α = @lift (own ? 1.0 : 0.45) * visibility($clock, appear, own ? disappear_cell : disappear_rest)
        lines!(ax, (cell - 0.5) .+ ξ, evalpoly.(ξ, Ref(a)); color = faded(color, α), linewidth = own ? 5 : 2.5,
               linestyle = dashed ? :dash : :solid)
    end

    # The face values, each labelled to the left or right of its point and above or below it.
    α = visible(clock, appear + 0.5, disappear_cell)
    for (ξ, (horizontal, vertical)) in zip((0.0, 1.0), labels)
        y = evalpoly(ξ, a)
        scatter!(ax, [cell - 0.5 + ξ], [y]; color = faded(color, α), markersize = 13)
        offset = (horizontal == :left ? -8 : 8, vertical == :above ? 8 : -8)
        text!(ax, cell - 0.5 + ξ, y; text = replace(string(round(y; digits = 3)), "-" => "−"), offset,
              align = (horizontal == :left ? :right : :left, vertical == :above ? :bottom : :top),
              fontsize = fontsize - 2, color = faded(color, α), font = :bold)
    end

    return nothing
end

"""
    worked_example_axis!(position, clock; title, fontsize = 22)

Draw every element of the worked example in a new axis at `position`, and return the axis. The observable `clock`, a
time in seconds on the timeline, sets what is visible.
"""
function worked_example_axis!(position, clock; title, fontsize = 22)
    ax = Axis(position; title, limits = (0.35, 6.65, -0.78, 1.72), xticks = (1:6, ["cell $j" for j in 1:6]),
              yticks = 0:0.5:1.5, ylabel = "u", titlesize = fontsize, titlefont = :regular)

    for band in BANDS
        band!(ax, clock; band...)
    end

    vlines!(ax, 0.5:1:6.5; color = FAINT, linewidth = 1)
    hlines!(ax, 0; color = FAINT, linewidth = 1)

    for (j, ū) in enumerate(AVERAGES)
        lines!(ax, [j - 0.5, j + 0.5], [ū, ū]; color = faded(COLORS.average, @lift(0.6 * $(visible(clock, 0)))),
               linewidth = 9)
        text!(ax, j, ū; text = string(ū), offset = (0, 12), align = (:center, :bottom), fontsize = fontsize - 2,
              color = faded(MUTED, visible(clock, 0.5)))
    end

    for candidate in CANDIDATES
        candidate!(ax, clock; fontsize, candidate...)
    end

    for parabola in PARABOLAS
        parabola!(ax, clock; fontsize, parabola...)
    end

    # The jump between the two reconstructions at the face between cells 3 and 4.
    top, bottom = right_face_value(average, 3, k, 2), left_face_value(average, 4, k, 0)
    α = visible(clock, 36.5)
    arrows2d!(ax, [3.5], [bottom + 0.06], [0.0], [top - bottom - 0.12]; color = faded(INK, α), shaftwidth = 2.5,
              tipwidth = 14, tiplength = 14)
    text!(ax, 3.62, (top + bottom) / 2; text = "difference −1.333\n(averages: −1.3)", align = (:left, :center),
          fontsize, color = faded(INK, α))

    return ax
end

#####
##### The animation and the still
#####

set_theme!(ARTICLE_THEME)

fig = Figure(size = (1280, 720), figure_padding = (14, 24, 10, 14))
clock = Observable(0.0)
worked_example_axis!(fig[1, 1], clock; title = @lift caption($clock))

frames = 0:round(Int, DURATION * FRAMERATE)
save_animation(n -> clock[] = n / FRAMERATE, fig, frames, "worked_example"; framerate = FRAMERATE)

# The still: the two decisions for cell 3, and the result for both cells.
still = Figure(size = (1900, 640), figure_padding = (10, 24, 10, 10))
panels = (("Cell 3, step 1: cell 2 joins", 9.0), ("Cell 3, step 2: cell 1 joins", 14.0),
          ("Result for cells 3 and 4", 42.0))

for (column, (title, t)) in enumerate(panels)
    ax = worked_example_axis!(still[1, column], Observable(t); title, fontsize = 20)
    column > 1 && hideydecorations!(ax; ticks = false, ticklabels = false)
end

colgap!(still.layout, 24)
save_figure(still, "worked_example")
