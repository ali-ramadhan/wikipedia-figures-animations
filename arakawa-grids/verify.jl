# Checks of the five schemes in arakawa_grids.jl:
#
# 1. Every scheme's inertia–gravity eigenmodes rotate at the frequency of its dispersion relation (the relations of
#    Arakawa & Lamb 1977, eqs. 28–31, in two dimensions as in Cushman-Roisin & Beckers 2011, Table 9.1).
# 2. The Runge–Kutta solution for the mound matches the exact solution of the semi-discrete equations, and conserves
#    mass and the semi-discrete energy.
# 3. The exact solution of the continuous equations starts from the mound, and its grid resolves it.
#
#     julia --project verify.jl

include("arakawa_grids.jl")

using LinearAlgebra

results = Bool[]

function check(label, passed)
    push!(results, passed)
    println(passed ? "PASS  " : "FAIL  ", label)
    return nothing
end

#####
##### 1. Eigenmodes rotate at the frequency of the dispersion relation
#####

# The Fourier symbols of the differences in the pressure gradient and the divergence, divided by i, in terms of the
# half-angles a = kΔ/2 and b = lΔ/2. The sum of their squares is the effective_wavenumber² of arakawa_grids.jl.
difference_symbols(::Union{AGrid, EGrid}, a, b, Δ) = (sin(2a), sin(2b)) ./ Δ
difference_symbols(::BGrid, a, b, Δ) = (2sin(a) * cos(b), 2cos(a) * sin(b)) ./ Δ
difference_symbols(::CGrid, a, b, Δ) = (2sin(a), 2sin(b)) ./ Δ
difference_symbols(::DGrid, a, b, Δ) = (sin(2a) * cos(b), cos(a) * sin(2b)) ./ Δ

"""
    symbol_matrix(grid, k, l, p)

The matrix 𝕄 of the semi-discrete equations for a plane wave with wavenumbers (k, l), d/dt (û, v̂, ĥ) = 𝕄 (û, v̂, ĥ),
with each amplitude defined at its own variable's points.
"""
function symbol_matrix(grid, k, l, p::Problem)
    Δ = spacing(p, grid)
    a, b = k * Δ / 2, l * Δ / 2
    α = grid isa Union{CGrid, DGrid} ? cos(a) * cos(b) : 1.0
    px, py = difference_symbols(grid, a, b, Δ)
    f, g, H = p.f, gravity(p), p.H

    return ComplexF64[
            0         f*α       -g*im*px
            -f*α      0         -g*im*py
            -H*im*px  -H*im*py  0
    ]
end

"""
    plane_wave!(state, grid, k, l, amplitudes, ν, t, p)

Fill `state` with the real part of the plane wave (û, v̂, ĥ) exp[i(kx + ly - νt)], at each variable's own points.
"""
function plane_wave!(state, grid, k, l, amplitudes, ν, t, p::Problem)
    u_offset, v_offset = velocity_offsets(grid)
    variables = (
            (state.u, coordinates(p, grid, u_offset), is_velocity_point),
            (state.v, coordinates(p, grid, v_offset), is_velocity_point),
            (state.h, coordinates(p, grid), is_height_point)
    )

    n = lattice_size(p, grid)
    for ((array, (x, y), carries), amplitude) in zip(variables, amplitudes)
        for j in 1:n, i in 1:n
            array[i, j] = carries(grid, i, j) ? real(amplitude * cis(k * x[i] + l * y[j] - ν * t)) : 0.0
        end
    end

    return state
end

"""
    eigenmode_error(grid, k, l, p; Δt = 0.0025, stop_time = 3.0)

Step the inertia–gravity eigenmode of `grid` with wavenumbers (k, l) to `stop_time`, and return the larger of its
error, relative to 1 + max |u|, and the difference between its frequency and `frequency(grid, k, l, p)`.
"""
function eigenmode_error(grid, k, l, p::Problem; Δt = 0.0025, stop_time = 3.0)
    eigenvalues, eigenvectors = eigen(symbol_matrix(grid, k, l, p))
    predicted_ν = frequency(grid, k, l, p)

    # The eigenvalue -iν with ν > 0 is the inertia–gravity wave exp(-iνt).
    mode = argmin(abs.(eigenvalues .+ im * predicted_ν))
    ν = -imag(eigenvalues[mode])
    amplitudes = eigenvectors[:, mode] ./ eigenvectors[3, mode]

    n = lattice_size(p, grid)
    state = plane_wave!(State(n), grid, k, l, amplitudes, ν, 0.0, p)
    advance!(state, Stepper(p, grid), Δt, round(Int, stop_time / Δt))
    exact = plane_wave!(State(n), grid, k, l, amplitudes, ν, stop_time, p)

    misfit = maximum(abs, state.h - exact.h) + maximum(abs, state.u - exact.u) + maximum(abs, state.v - exact.v)
    return max(misfit / (1 + maximum(abs, exact.u)), abs(ν - predicted_ν))
end

println("--- 1. eigenmodes rotate at the predicted frequency (small domain, many wave vectors) ---")

for λ in (2.0, 0.5)
    p = Problem(; λ, N = 32, M = 46)
    L = domain_length(p)

    for grid in GRIDS
        m = lattice_size(p, grid) ÷ 2   # the largest wavenumber is 2πm/L
        modes = ((1, 0), (3, 2), (m ÷ 2, 1), (m - 1, m ÷ 3), (m ÷ 2, m ÷ 2), (m - 2, m - 3))
        worst = maximum(eigenmode_error(grid, 2π * mk / L, 2π * ml / L, p) for (mk, ml) in modes)

        label = @sprintf("λ/d = %.1f, grid %s: eigenmodes follow ν² = α²f² + gH p² (worst error %.1e)",
                λ,
                grid,
                worst
        )

        check(label, worst < 1e-7)
    end
end

#####
##### 2. The mound on each grid
#####

println("--- 2. Runge–Kutta solution from the Gaussian = exact semi-discrete solution (production grids) ---")

for λ in (2.0, 0.5)
    p = Problem(; λ)

    for grid in GRIDS
        state = initial_state(p, grid)
        E₀ = energy(state, p, grid)
        m₀ = total_mass(state, p, grid)

        Δt, stop_time = 0.01, 10.0
        advance!(state, Stepper(p, grid), Δt, round(Int, stop_time / Δt))

        h_error = maximum(abs, state.h - semidiscrete_h(p, grid, stop_time))
        energy_change = abs(energy(state, p, grid) - E₀) / E₀
        mass_change = abs(total_mass(state, p, grid) - m₀) / m₀

        label = @sprintf("λ/d = %.1f, grid %s, t = %g/f: ", λ, grid, stop_time) *
                @sprintf("max |h - h_exact| = %.1e, energy change %.1e, mass change %.1e",
                        h_error,
                        energy_change,
                        mass_change
                )

        check(label, h_error < 1e-7 && energy_change < 1e-8 && mass_change < 1e-12)

        if grid isa EGrid
            no_heights = all(state.h[i, j] == 0 for j in 1:p.M, i in 1:p.M if is_velocity_point(grid, i, j))
            check("grid E: no height ever appears on the velocity points", no_heights)
        end
    end
end

#####
##### 3. The exact solution of the continuous equations
#####

println("--- 3. the exact continuous solution ---")

p = Problem(; λ = 2.0)

x, y, h = exact_h(p, 0.0)
gaussian = [exp(-(xᵢ^2 + yⱼ^2) / 2) for xᵢ in x, yⱼ in y]
check("t = 0 gives back the Gaussian", maximum(abs, h .- gaussian) < 1e-12)

_, _, h₄ = exact_h(p, 10.0; refine = 4)
_, _, h₈ = exact_h(p, 10.0; refine = 8)
change = maximum(abs, h₄ - h₈[1:2:end, 1:2:end])
check(@sprintf("refining the grid from d/4 to d/8 changes h at t = 10/f by %.1e", change), change < 1e-12)

println("\n$(count(results))/$(length(results)) checks passed")
