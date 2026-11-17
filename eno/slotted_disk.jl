# Zalesak's slotted disk turning counterclockwise about the center of the domain −1 ≤ x, y ≤ 1, on 100 × 100 cells, in
# a 2 × 2 grid: first-order upwind (all numerical diffusion) beside the second-order centered scheme (all numerical
# dispersion), above a fixed third-order stencil beside third-order ENO. The equation is
# ∂u/∂t + ∂(au)/∂x + ∂(bu)/∂y = 0 with (a, b) = ω(−y, x), so along each grid row the speed a is constant and along
# each column b is: each scheme is applied one grid line at a time, as in the finite-difference form of ENO in several
# dimensions, with `rotation_tendency!` of eno.jl and third-order SSP Runge–Kutta time stepping. The disk has
# radius 0.3 and center (0, 0.5), with a slot 0.1 wide cut from its edge nearest the center up to y = 0.7, as in
# Zalesak (1979) scaled to this domain. The shading runs from 0 (white) to 1 (gray); values more than TOLERANCE above
# 1 are red and more than TOLERANCE below 0 are blue, so new extrema stand out. The thin black line is the exact
# outline. Prints the errors after each turn, then writes slotted_disk.webm, its lossless master slotted_disk.mkv, and
# the last frame as slotted_disk.png. STILL_ONLY=1 writes only the still.
#
#     julia --project slotted_disk.jl

include("eno.jl")
include("plotting.jl")

const N = 100
const TURNS = 2
const ω = π                  # one turn every 2 time units
const PERIOD = 2π / ω
const CFL = 0.4
const STEPS_PER_FRAME = 6
const FRAMERATE = 25
const HOLD = 2.0
const TOLERANCE = 0.01

const CENTER = (0.0, 0.5)
const RADIUS = 0.3
const SLOT_HALF_WIDTH = 0.05
const SLOT_TOP = 0.7

Δx = 2 / N
xs = [-1 + (i - 0.5) * Δx for i in 1:N]

in_disk(x, y) = x^2 + (y - CENTER[2])^2 ≤ RADIUS^2 && !(abs(x) ≤ SLOT_HALF_WIDTH && y ≤ SLOT_TOP)

"""The cell averages of the slotted disk, from 16 × 16 samples in each cell."""
function initial_condition()
    s = ((1:16) .- 0.5) ./ 16 .- 0.5
    return [count(in_disk(x + Δx * a, y + Δx * b) for a in s, b in s) / 256 for x in xs, y in xs]
end

"""The outline of the slotted disk turned through the angle θ, as a closed polygon."""
function outline(θ)
    yₛ = CENTER[2] - sqrt(RADIUS^2 - SLOT_HALF_WIDTH^2)   # where the slot's sides meet the circle
    φ₀ = atan(yₛ - CENTER[2], SLOT_HALF_WIDTH)            # the angle of the right one, below the center
    φ₁ = φ₀ + 2π - 2(φ₀ + π / 2)                           # the left one, the long way round through the top
    arc = [(CENTER[1] + RADIUS * cos(φ), CENTER[2] + RADIUS * sin(φ)) for φ in range(φ₀, φ₁; length = 400)]
    points = [arc; (-SLOT_HALF_WIDTH, SLOT_TOP); (SLOT_HALF_WIDTH, SLOT_TOP); arc[1]]
    return [Point2f(x * cos(θ) - y * sin(θ), x * sin(θ) + y * cos(θ)) for (x, y) in points]
end

u₀ = initial_condition()
Δt_max = CFL * Δx / (ω * 2maximum(abs, xs))
steps_per_turn = STEPS_PER_FRAME * ceil(Int, PERIOD / (Δt_max * STEPS_PER_FRAME))
Δt = PERIOD / steps_per_turn
steps = TURNS * steps_per_turn

schemes = (("Upwind (1st order)", FirstOrderUpwind()), ("Centered (2nd order)", FixedStencil(2, 0)),
           ("Fixed stencil (3rd order)", FixedStencil(3, 1)), ("ENO (3rd order)", ENO(3)))

println("$steps_per_turn steps of Δt = $(round(Δt; sigdigits = 3)) per turn, $TURNS turns")
println("L1 error (the mean of |u − exact|), highest and lowest value, after each turn:")

fields = map(schemes) do (name, scheme)
    function report_turn(u, n)
        if n > 0 && n % steps_per_turn == 0
            @printf("  %-26s turn %d: L1 %.4f, max %.4f, min %.4f\n", name, n ÷ steps_per_turn,
                    sum(abs.(u .- u₀)) / N^2, maximum(u), minimum(u))
        end
        return nothing
    end

    tendency!(du, u) = rotation_tendency!(du, u, scheme, xs, ω, Δx)
    seconds = @elapsed snapshots = evolve(u₀, tendency!; Δt, steps, save_every = STEPS_PER_FRAME,
                                          callback = report_turn)
    @printf("  %-26s %.1f s\n", name, seconds)

    [Float32.(u) for u in snapshots]
end

times = (0:length(fields[1])-1) .* (STEPS_PER_FRAME * Δt)
colormap = cgrad([PAPER, colorant"#7f7e78"])
const TICKS = labelled_ticks([-1, 0, 1])

set_theme!(ARTICLE_THEME)

fig = Figure(size = (1160, 1120), fontsize = 22, figure_padding = (12, 20, 10, 12))
frame = Observable(1)
t = @lift times[$frame]

Label(fig[0, 1:2], @lift(@sprintf("Rotating slotted disk: %.2f turns", $t / PERIOD)); fontsize = 26, tellwidth = false)

for (n, (name, _)) in enumerate(schemes)
    row, column = fldmod1(n, 2)
    ax = Axis(fig[row, column]; title = name, titlealign = :center, aspect = DataAspect(), limits = (-1, 1, -1, 1),
              xticks = TICKS, yticks = TICKS, xlabel = italic("x"), ylabel = italic("y"), FULL_FRAME...)
    row == 1 && hidexdecorations!(ax; ticks = false)
    column == 2 && hideydecorations!(ax; ticks = false)
    heatmap!(ax, xs, xs, @lift(fields[n][$frame]); colormap, colorrange = (-TOLERANCE, 1 + TOLERANCE), lowclip = BLUE,
             highclip = RED)
    lines!(ax, @lift(outline(ω * $t)); color = COLORS.exact, linewidth = 1.2)
end

Colorbar(fig[1:2, 3]; colormap, limits = (-TOLERANCE, 1 + TOLERANCE), lowclip = BLUE, highclip = RED,
         height = Relative(0.6), ticks = labelled_ticks([0, 0.5, 1]),
         label = rich(italic("u"), @sprintf("; red above %g, blue below −%g", 1 + TOLERANCE, TOLERANCE)))

# The first gap keeps the tick labels −1 and 1 of neighbouring panels apart.
colgap!(fig.layout, 1, 44)
colgap!(fig.layout, 2, 20)
rowgap!(fig.layout, 8)

save_animation(n -> frame[] = n, fig, eachindex(times), "slotted_disk"; framerate = FRAMERATE, hold_seconds = HOLD)

frame[] = lastindex(times)
save_figure(fig, "slotted_disk"; formats = ("png",))
