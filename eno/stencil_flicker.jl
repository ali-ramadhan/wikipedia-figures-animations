# The stencil third-order ENO chooses at each cell, step by step, as a smooth pulse is advected once round a periodic
# domain on 100 cells: on the left the pulse alone, on the right with noise of size 10⁻⁴ added. Above, the solution;
# below, the shift r of each cell's stencil, building up into a map in x and t. In the pulse's nearly flat tails the
# differences being compared are tiny, so there the noise decides the choice. Writes stencil_flicker.webm, its
# lossless master stencil_flicker.mkv, and the last frame as stencil_flicker.png. STILL_ONLY=1 writes only the still.
#
#     julia --project stencil_flicker.jl

include("eno.jl")
include("plotting.jl")

using Random

const N = 100
const T = 2.0
const NOISE = 1e-4
const FRAMERATE = 30
const HOLD = 2.0

Δx = 2 / N
x = [-1 + (j - 0.5) * Δx for j in 1:N]
steps, Δt = time_steps(T, 0.4Δx)
clean = exp.(-60 .* x .^ 2)
noisy = clean .+ NOISE .* rand(Xoshiro(1), N)

"""Advect `u₀` with third-order ENO, and return every cell's stencil shift r and the solution after each step."""
function shifts_and_solutions(u₀)
    shifts = fill(NaN, steps + 1, N)
    solutions = Vector{Vector{Float64}}()

    function record!(u, n)
        shifts[n + 1, :] .= [eno_shift(periodic(u), i, 3) for i in 1:N]
        push!(solutions, copy(u))
        return nothing
    end

    advect(u₀, ENO(3); Δx, Δt, steps, callback = record!)
    return shifts, solutions
end

cases = [("Smooth pulse", shifts_and_solutions(clean)...), ("Plus noise of size 10⁻⁴", shifts_and_solutions(noisy)...)]
times = [n * Δt for n in 0:steps]

for (name, shifts, _) in cases
    changes = count(shifts[2:end, :] .!= shifts[1:end-1, :]) / length(shifts[2:end, :])
    @printf("%s: the stencil changes at %.1f%% of cell steps\n", name, 100changes)
end

set_theme!(ARTICLE_THEME)

fig = Figure(size = (1280, 860), figure_padding = (10, 24, 10, 10))
frame = Observable(1)
Label(fig[0, 1:2], @lift(rich("ENO (3rd order) stencil choice when advecting a smooth pulse: ", italic("t"),
                              @sprintf(" = %.2f", times[$frame]))); fontsize = 24, tellwidth = false)
colormap = cgrad(collect(SHIFT_COLORS), 3; categorical = true)

# Short tick labels, so that the two columns' labels don't run together.
for (column, (name, shifts, solutions)) in enumerate(cases)
    solution_axis = Axis(fig[1, column]; title = name, titlealign = :center, limits = (-1, 1, -0.2, 1.2),
                         xticks = labelled_ticks(-1:0.5:1), yticks = labelled_ticks(0:0.5:1),
                         ylabel = column == 1 ? italic("u") : "", FULL_FRAME...)
    hidexdecorations!(solution_axis; ticks = false)
    lines!(solution_axis, x, @lift(solutions[$frame]); color = INK, linewidth = 3)

    # The stencils chosen so far, as a map in x and t.
    stencil_axis = Axis(fig[2, column]; xlabel = italic("x"), ylabel = column == 1 ? italic("t") : "",
                        limits = (-1, 1, 0, T), xticks = labelled_ticks(-1:0.5:1), yticks = labelled_ticks(0:0.5:2),
                        FULL_FRAME...)
    shown = @lift [n ≤ $frame ? shifts[n, i] : NaN for i in 1:N, n in 1:steps+1]
    heatmap!(stencil_axis, x, times, shown; colormap, colorrange = (-0.5, 2.5), nan_color = :transparent)

    if column == 2
        hideydecorations!(solution_axis; ticks = false)
        hideydecorations!(stencil_axis; ticks = false)
    end
end

entries = [PolyElement(color = SHIFT_COLORS[r + 1]) for r in (2, 1, 0)]
Legend(fig[3, 1:2], entries, ["Cells i−2 to i (r = 2)", "Cells i−1 to i+1 (r = 1)", "Cells i to i+2 (r = 0)"];
       orientation = :horizontal, patchsize = (30, 18), colgap = 30)
rowsize!(fig.layout, 1, Relative(0.22))
rowgap!(fig.layout, 12)
colgap!(fig.layout, 50)

save_animation(n -> frame[] = n, fig, eachindex(times), "stencil_flicker"; framerate = FRAMERATE, hold_seconds = HOLD)

frame[] = lastindex(times)
save_figure(fig, "stencil_flicker"; formats = ("png",))
