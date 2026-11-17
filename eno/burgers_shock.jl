# Burgers' equation from the smooth initial data u₀(x) = 1/2 + sin(πx) on the periodic domain −1 ≤ x ≤ 1, through the
# formation of a shock at t = 1/π to t = 2, on 80 cells. A fixed third-order stencil (one cell on either side) is shown
# beside third-order ENO, both finite-volume schemes with the Godunov flux and third-order strong-stability-preserving
# Runge–Kutta time stepping, against the exact solution. The y range reaches 2.25 so that the fixed stencil's largest
# overshoot, 2.11, stays in view. Writes burgers_shock.webm, its lossless master burgers_shock.mkv, and the frame at
# t = 0.67, soon after the shock forms, as burgers_shock.png. STILL_ONLY=1 writes only the still.
#
#     julia --project burgers_shock.jl

include("eno.jl")
include("plotting.jl")

const N = 80
const T = 2.0
const FRAMERATE = 25
const HOLD = 2.0
const T_SHOCK = 1 / π
const STILL_TIME = 0.67   # soon after the shock forms, when the fixed stencil's overshoot is large

Δx = 2 / N
edges = range(-1, 1; length = N + 1)
x = [(edges[j] + edges[j + 1]) / 2 for j in 1:N]
u₀ = [cell_average(x -> 0.5 + sinpi(x), edges[j], edges[j + 1]) for j in 1:N]

"""The exact cell averages at time `t`, with the quadrature split at the shock."""
exact_averages(t) = [cell_average(x -> burgers_exact(x, t), edges[j], edges[j + 1];
                                  jumps = t > T_SHOCK ? (burgers_shock_position(t),) : ()) for j in 1:N]

steps, Δt = time_steps(T, 0.4Δx / 1.5)   # |u| ≤ 1.5
times = (0:steps) .* Δt

schemes = (("Fixed stencil (3rd order)", FixedStencil(3, 1), COLORS.fixed), ("ENO (3rd order)", ENO(3), COLORS.eno))
solutions = [evolve(u₀, (du, u) -> burgers_tendency!(du, u, scheme, Δx); Δt, steps, save_every = 1)
             for (_, scheme, _) in schemes]

# The exact solution stays between −1/2 and 3/2, so anything outside is a new extremum.
exact_final = exact_averages(T)
println("$steps steps of Δt = $(round(Δt; sigdigits = 3)); the shock forms at t = $(round(T_SHOCK; digits = 3))")
for ((name, _, _), snapshots) in zip(schemes, solutions)
    peak = argmax(maximum.(snapshots))
    @printf("  %-26s highest %.4f at t = %.2f\n", name, maximum(snapshots[peak]), times[peak])
    @printf("  %-26s max %.4f, min %.4f at any time; L1 error at t = %.0f: %.2e; total change %.1e\n", name,
            maximum(maximum, snapshots), minimum(minimum, snapshots), T, sum(abs.(snapshots[end] .- exact_final)) * Δx,
            abs(sum(snapshots[end]) - sum(u₀)) * Δx)
end

fine = range(-1, 1; length = 4001)

set_theme!(ARTICLE_THEME)

fig = Figure(size = (1280, 720), fontsize = 24, figure_padding = (12, 28, 10, 12))
frame = Observable(1)
t = @lift times[$frame]

Label(fig[0, 1], @lift(rich("Burgers' equation: ", italic("t"), @sprintf(" = %.2f", $t))); fontsize = 28,
      tellwidth = false)

ax = Axis(fig[2, 1]; xlabel = italic("x"), ylabel = italic("u"), limits = (-1, 1, -0.8, 2.25),
          xticks = labelled_ticks(-1:0.5:1), yticks = labelled_ticks(-0.5:0.5:2), FULL_FRAME...)
lines!(ax, fine, @lift(burgers_exact.(fine, $t)); color = COLORS.exact, linewidth = 1.5)

for ((_, _, color), snapshots) in zip(schemes, solutions)
    lines!(ax, x, @lift(snapshots[$frame]); color, linewidth = 4)
end

entries = [LineElement(color = COLORS.exact, linewidth = 1.5),
           [LineElement(color = c, linewidth = 4) for (_, _, c) in schemes]...]
Legend(fig[1, 1], entries, ["Exact", first.(schemes)...]; orientation = :horizontal, patchsize = (40, 14), colgap = 30)

save_animation(n -> frame[] = n, fig, eachindex(times), "burgers_shock"; framerate = FRAMERATE, hold_seconds = HOLD)

frame[] = argmin(abs.(times .- STILL_TIME))
save_figure(fig, "burgers_shock"; formats = ("png",))
