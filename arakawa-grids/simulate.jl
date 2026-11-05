# Geostrophic adjustment of a narrow mound of fluid, simulated on each of the five Arakawa grids in two cases, as in
# Arakawa & Lamb (1977, pp. 184–190): the radius of deformation λ twice the grid spacing (well resolved) and half of it
# (poorly resolved). For each case and grid, it saves the height in the plotted window at every animation frame, with
# the height at the center of the mound, to adjustment_<case>.jld2.
#
#     julia --project simulate.jl

include("arakawa_grids.jl")

using JLD2

cases = [
        (name = "fine", λ = 2.0, stop_time = 15.0, frame_interval = 0.05, window = 32.0),
        (name = "coarse", λ = 0.5, stop_time = 40.0, frame_interval = 0.125, window = 20.0)
]

# In units of 1/f, with Runge–Kutta errors of order 10⁻⁹ (verify.jl).
const Δt = 0.005

"""
    window_indices(p, grid, W)

The array indices of the lattice points within [-W, W]² of the center, in each direction, and their coordinates.
"""
function window_indices(p::Problem, grid, W)
    n = lattice_size(p, grid)
    s = spacing(p, grid)
    i₀, j₀ = center(p, grid)
    r = floor(Int, W / s) + 1

    I = [mod1(i₀ + q, n) for q in -r:r]
    J = [mod1(j₀ + q, n) for q in -r:r]
    return I, J, (-r:r) .* s
end

function simulate(case)
    p = Problem(; λ = case.λ)
    steps_per_frame = round(Int, case.frame_interval / Δt)
    @assert steps_per_frame * Δt ≈ case.frame_interval
    frames = round(Int, case.stop_time / case.frame_interval) + 1
    times = collect((0:frames-1) .* case.frame_interval)

    output = Dict{String, Any}(
            "t" => times,
            "λ" => case.λ,
            "σ" => p.σ,
            "N" => p.N,
            "M" => p.M,
            "Δt" => Δt,
            "window" => case.window
    )

    for grid in GRIDS
        I, J, offsets = window_indices(p, grid, case.window)
        i₀, j₀ = center(p, grid)
        state = initial_state(p, grid)
        stepper = Stepper(p, grid)
        E₀ = energy(state, p, grid)

        h = zeros(Float32, length(I), length(J), frames)
        h_center = zeros(frames)

        wall_time = @elapsed for frame in 1:frames
            frame > 1 && advance!(state, stepper, Δt, steps_per_frame)
            h[:, :, frame] .= state.h[I, J]
            h_center[frame] = state.h[i₀, j₀]
        end

        # The run must agree with the exact solution of its semi-discrete equations (verify.jl).
        h_error = maximum(abs, state.h - semidiscrete_h(p, grid, times[end]))
        energy_drift = abs(energy(state, p, grid) - E₀) / E₀

        @printf("%-6s grid %s: %d frames to ft = %g in %.1f s; |h - exact semi-discrete| = %.1e, energy drift %.1e\n",
                case.name,
                grid,
                frames,
                times[end],
                wall_time,
                h_error,
                energy_drift
        )

        output["h_$grid"] = h
        output["offsets_$grid"] = collect(offsets)
        output["center_$grid"] = h_center
    end

    path = output_file("adjustment_$(case.name).jld2")
    jldsave(path; (Symbol(key) => value for (key, value) in output)...)
    println("saved $path")

    return nothing
end

foreach(simulate, cases)
