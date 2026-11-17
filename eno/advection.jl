# Advection of the profile of Jiang & Shu (1996), Example 1, four times round a periodic domain on 200 cells, with
# twelve schemes in three rows: first-order upwind, second-order centered, and the second-order TVD scheme with the
# minmod and superbee limiters; fixed stencils of orders 2 to 5; and ENO of orders 2 to 5. All use the same third-order
# Runge–Kutta time stepping with Δt = 0.4Δx. The y range is kept from −0.38 to 1.38, so the centered scheme's largest
# swings leave its panel. Prints how far each scheme strays outside the range of the exact solution, then writes
# advection.webm, its lossless master advection.mkv, and the last frame as advection.svg and advection.png. STILL_ONLY=1
# writes only the still.
#
#     julia --project advection.jl

include("eno.jl")
include("plotting.jl")

const N = 200
const T = 8.0
const SAVE_EVERY = 3   # steps per frame, 0.012 time units
const FRAMERATE = 25   # so the profile moves 0.3 units, or 30 cells, per second
const HOLD = 2.0       # seconds on the last frame

Δx = 2 / N
x = [-1 + (j - 0.5) * Δx for j in 1:N]
u₀ = jiang_shu_profile.(x)
steps, Δt = time_steps(T, 0.4Δx)

ordinal(k) = k == 1 ? "1st" : k == 2 ? "2nd" : k == 3 ? "3rd" : "$(k)th"

# Each fixed stencil leans upwind, with r = ⌊k/2⌋ of its k cells to the left of the target cell. For odd k it is
# centred on the target cell, which gives the upwind-biased schemes of odd order; for k = 2 it is the classical
# second-order upwind scheme. The centered scheme is the other second-order stencil, r = 0, which is symmetric about the
# face; second-order ENO chooses between the two in each cell.
schemes = [
        ("Upwind (1st order)", FirstOrderUpwind(), COLORS.upwind),
        ("Centered (2nd order)", FixedStencil(2, 0), COLORS.centered),
        ("TVD minmod (2nd order)", Minmod(), COLORS.tvd),
        ("TVD superbee (2nd order)", Superbee(), COLORS.tvd),
        [("Fixed stencil ($(ordinal(k)) order)", FixedStencil(k, k ÷ 2), COLORS.fixed) for k in 2:5]...,
        [("ENO ($(ordinal(k)) order)", ENO(k), COLORS.eno) for k in 2:5]...
]

seconds = @elapsed solutions = [advect(u₀, scheme; Δx, Δt, steps, save_every = SAVE_EVERY)
                                for (_, scheme, _) in schemes]
println("Solved in $(round(seconds; digits = 1)) s")
times = unique([0:SAVE_EVERY:steps; steps]) .* Δt
@assert all(length(s) == length(times) for s in solutions)

# The exact solution stays between 0 and 1, so any value outside that range is a new maximum or minimum. At t = 8 the
# profile is back where it started. Each wave's peak is its largest value at t = 8.
waves = (("Gaussians", -0.8, -0.6), ("square", -0.4, -0.2), ("triangle", 0.0, 0.2), ("ellipse", 0.4, 0.6))
println("\nlargest value and smallest value at any time; L1 error and the peak of each wave at t = 8")
for ((name, _, _), snapshots) in zip(schemes, solutions)
    final = snapshots[end]
    peaks = [maximum(final[j] for j in eachindex(x) if a ≤ x[j] ≤ b) for (_, a, b) in waves]
    @printf("  %-26s %8.4f %9.1e   %.4f   %s\n", name, maximum(maximum, snapshots), minimum(minimum, snapshots),
            sum(abs.(final .- u₀)) / N, join((@sprintf("%.3f", p) for p in peaks), " "))
end

fine = range(-1, 1; length = 2001)
exact(t) = jiang_shu_profile.(mod.(fine .- t .+ 1, 2) .- 1)

set_theme!(ARTICLE_THEME)

fig = Figure(size = (1280, 800), fontsize = 17, figure_padding = (10, 22, 8, 10))
frame = Observable(1)
t = @lift times[$frame]
Label(fig[0, 1:4], @lift(rich("Advection to the right: ", italic("t"), @sprintf(" = %.2f", $t))); fontsize = 22,
      tellwidth = false)

for (n, (name, _, color)) in enumerate(schemes)
    row, column = fldmod1(n, 4)
    ax = Axis(fig[row, column]; title = name, titlealign = :center, limits = (-1, 1, -0.38, 1.38),
              xticks = labelled_ticks(-1:0.5:1), yticks = labelled_ticks(0:0.5:1), xlabel = italic("x"),
              ylabel = italic("u"), FULL_FRAME...)
    row < 3 && hidexdecorations!(ax; ticks = false)
    column > 1 && hideydecorations!(ax; ticks = false)
    lines!(ax, fine, @lift(exact($t)); color = COLORS.exact, linewidth = 1.25)
    lines!(ax, x, @lift(solutions[n][$frame]); color, linewidth = 2.5)
end

rowgap!(fig.layout, 10)
colgap!(fig.layout, 34)

save_animation(n -> frame[] = n, fig, eachindex(times), "advection"; framerate = FRAMERATE, hold_seconds = HOLD)

frame[] = lastindex(times)
save_figure(fig, "advection")
