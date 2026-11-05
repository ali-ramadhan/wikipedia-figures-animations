# Shared code for the linearized shallow-water equations on the five Arakawa grids: the grids, the schemes, the time
# stepping, and the exact solutions that check them. The README states the equations.
#
# Everything is nondimensional: lengths in units of the grid spacing d, time in units of 1/f, and H = 1. Gravity is set
# so that the radius of deformation λ = √(gH)/f takes the chosen value, so λ means λ/d throughout.

using FFTW
using Printf

output_directory = mkpath(get(ENV, "OUTPUT_DIRECTORY", joinpath(@__DIR__, "output")))
output_file(name) = joinpath(output_directory, name)

"""
    Problem(; λ, σ = 1, N = 198, M = 280)

A mound of fluid, of height exp(-r²/2σ²) and at rest, at the center of a doubly periodic square, with the radius of
deformation λ. Grids A–D have `N × N` height points with spacing d = 1, so the square has side L = N. Grid E puts
heights and velocities on alternate points (a checkerboard) of an `M × M` lattice with spacing s = L/M, and its
differences span d* = 2s. With M ≈ √2 N, d* ≈ √2 d and grid E has as many height points as the others (A&L, p. 183):
N = 198 and M = 280 give d* = 1.41429 d, within 0.005% of √2 d. M must be even for the checkerboard to be periodic.
"""
Base.@kwdef struct Problem
    λ::Float64          # radius of deformation √(gH)/f
    σ::Float64 = 1.0    # width of the mound
    N::Int = 198        # height points per side, grids A–D
    M::Int = 280        # lattice points per side, grid E (half of them heights)
    f::Float64 = 1.0
    H::Float64 = 1.0
end

gravity(p::Problem) = p.λ^2 * p.f^2 / p.H
domain_length(p::Problem) = float(p.N)

#####
##### The five grids
#####

"""
    ArakawaGrid{letter}()

One of Arakawa's five grids, A to E, as a type to dispatch on. It prints as its letter.
"""
struct ArakawaGrid{letter} end

const AGrid = ArakawaGrid{:A}
const BGrid = ArakawaGrid{:B}
const CGrid = ArakawaGrid{:C}
const DGrid = ArakawaGrid{:D}
const EGrid = ArakawaGrid{:E}

const GRIDS = (AGrid(), BGrid(), CGrid(), DGrid(), EGrid())

Base.show(io::IO, ::ArakawaGrid{letter}) where letter = print(io, letter)

"""
    Exact()

The continuous equations, whose exact solutions the grids are compared with. It prints as "exact".
"""
struct Exact end

Base.show(io::IO, ::Exact) = print(io, "exact")

# Grid E's lattice has M points per side, with spacing s = L/M. The others have N, with spacing d = 1.
lattice_size(p::Problem, ::ArakawaGrid) = p.N
lattice_size(p::Problem, ::EGrid) = p.M

spacing(p::Problem, ::ArakawaGrid) = 1.0
spacing(p::Problem, ::EGrid) = domain_length(p) / p.M

# The area per height point, and per velocity point, for sums over the domain.
area_per_point(p::Problem, ::ArakawaGrid) = 1.0
area_per_point(p::Problem, grid::EGrid) = 2 * spacing(p, grid)^2

#####
##### Where each variable lives
#####

# The offsets of u and v from the height point with the same array index, in units of the lattice spacing (A&L, Fig. 3).
velocity_offsets(::AGrid) = ((0, 0), (0, 0))
velocity_offsets(::BGrid) = ((1/2, 1/2), (1/2, 1/2))
velocity_offsets(::CGrid) = ((1/2, 0), (0, 1/2))
velocity_offsets(::DGrid) = ((0, 1/2), (1/2, 0))
velocity_offsets(::EGrid) = ((0, 0), (0, 0))

# On grid E, all three arrays span the whole lattice, and the parity of i + j decides which variable a point carries:
# heights where i + j is odd, velocities where it is even.
is_height_point(::ArakawaGrid, i, j) = true
is_height_point(::EGrid, i, j) = isodd(i + j)

is_velocity_point(::ArakawaGrid, i, j) = true
is_velocity_point(::EGrid, i, j) = iseven(i + j)

"""
    center(p, grid)

The array index of the height point at the center of the domain.
"""
center(p::Problem, ::ArakawaGrid) = (p.N ÷ 2 + 1, p.N ÷ 2 + 1)
center(p::Problem, ::EGrid) = (p.M ÷ 2, p.M ÷ 2 + 1)

"""
    coordinates(p, grid, offset = (0, 0))

The x and y of the points of a variable with the given offset, measured from the central height point and wrapped into
[-L/2, L/2).
"""
function coordinates(p::Problem, grid, offset = (0, 0))
    n = lattice_size(p, grid)
    s = spacing(p, grid)
    L = domain_length(p)
    i₀, j₀ = center(p, grid)

    wrap(x) = mod(x + L / 2, L) - L / 2
    x = [wrap((i - i₀ + offset[1]) * s) for i in 1:n]
    y = [wrap((j - j₀ + offset[2]) * s) for j in 1:n]

    return x, y
end

#####
##### The five schemes (A&L, eqs. 5–19)
#####

# The periodic neighbors of i on a lattice of n points: the next one and the previous one.
@inline neighbors(i, n) = ifelse(i == n, 1, i + 1), ifelse(i == 1, n, i - 1)

"""
    tendencies!(du, dv, dh, u, v, h, p, grid)

Evaluate the right-hand sides of the semi-discrete equations on `grid` into `du`, `dv` and `dh`. The height h[i, j] is
at the point (i, j), and u and v are offset from it as in `velocity_offsets`.
"""
function tendencies! end

# Everything at the height points: centered differences over 2d.
function tendencies!(du, dv, dh, u, v, h, p::Problem, grid::AGrid)
    f, g, H = p.f, gravity(p), p.H
    n, Δ = lattice_size(p, grid), spacing(p, grid)

    @inbounds for j in 1:n, i in 1:n
        (i₊, i₋), (j₊, j₋) = neighbors(i, n), neighbors(j, n)
        du[i, j] = f * v[i, j] - g * (h[i₊, j] - h[i₋, j]) / 2Δ
        dv[i, j] = -f * u[i, j] - g * (h[i, j₊] - h[i, j₋]) / 2Δ
        dh[i, j] = -H * ((u[i₊, j] - u[i₋, j]) + (v[i, j₊] - v[i, j₋])) / 2Δ
    end

    return nothing
end

# u and v together at the cell corners (i + ½, j + ½): differences across the cell, averaged along it.
function tendencies!(du, dv, dh, u, v, h, p::Problem, grid::BGrid)
    f, g, H = p.f, gravity(p), p.H
    n, Δ = lattice_size(p, grid), spacing(p, grid)

    @inbounds for j in 1:n, i in 1:n
        (i₊, i₋), (j₊, j₋) = neighbors(i, n), neighbors(j, n)
        du[i, j] = f * v[i, j] - g * ((h[i₊, j] + h[i₊, j₊]) - (h[i, j] + h[i, j₊])) / 2Δ
        dv[i, j] = -f * u[i, j] - g * ((h[i, j₊] + h[i₊, j₊]) - (h[i, j] + h[i₊, j])) / 2Δ
        dh[i, j] = -H * (((u[i, j] + u[i, j₋]) - (u[i₋, j] + u[i₋, j₋])) +
                         ((v[i, j] + v[i₋, j]) - (v[i, j₋] + v[i₋, j₋]))) / 2Δ
    end

    return nothing
end

# u at (i + ½, j) and v at (i, j + ½): the pressure gradient and the divergence are differences over d, but the Coriolis
# force averages each velocity component over the four nearest points.
function tendencies!(du, dv, dh, u, v, h, p::Problem, grid::CGrid)
    f, g, H = p.f, gravity(p), p.H
    n, Δ = lattice_size(p, grid), spacing(p, grid)

    @inbounds for j in 1:n, i in 1:n
        (i₊, i₋), (j₊, j₋) = neighbors(i, n), neighbors(j, n)
        v̄ = (v[i, j] + v[i, j₋] + v[i₊, j] + v[i₊, j₋]) / 4
        ū = (u[i, j] + u[i₋, j] + u[i, j₊] + u[i₋, j₊]) / 4
        du[i, j] = f * v̄ - g * (h[i₊, j] - h[i, j]) / Δ
        dv[i, j] = -f * ū - g * (h[i, j₊] - h[i, j]) / Δ
        dh[i, j] = -H * ((u[i, j] - u[i₋, j]) + (v[i, j] - v[i, j₋])) / Δ
    end

    return nothing
end

# u at (i, j + ½) and v at (i + ½, j): grid C with u and v swapped. The Coriolis force is averaged as on grid C, and the
# pressure gradient and the divergence are differences over d of four-point averages.
function tendencies!(du, dv, dh, u, v, h, p::Problem, grid::DGrid)
    f, g, H = p.f, gravity(p), p.H
    n, Δ = lattice_size(p, grid), spacing(p, grid)

    @inbounds for j in 1:n, i in 1:n
        (i₊, i₋), (j₊, j₋) = neighbors(i, n), neighbors(j, n)
        v̄ = (v[i, j] + v[i₋, j] + v[i, j₊] + v[i₋, j₊]) / 4
        ū = (u[i, j] + u[i, j₋] + u[i₊, j] + u[i₊, j₋]) / 4
        ∂h∂x = ((h[i₊, j] - h[i₋, j]) + (h[i₊, j₊] - h[i₋, j₊])) / 4Δ
        ∂h∂y = ((h[i, j₊] - h[i, j₋]) + (h[i₊, j₊] - h[i₊, j₋])) / 4Δ
        du[i, j] = f * v̄ - g * ∂h∂x
        dv[i, j] = -f * ū - g * ∂h∂y
        dh[i, j] = -H * (((u[i₊, j] - u[i₋, j]) + (u[i₊, j₋] - u[i₋, j₋])) +
                         ((v[i, j₊] - v[i, j₋]) + (v[i₋, j₊] - v[i₋, j₋]))) / 4Δ
    end

    return nothing
end

# A checkerboard: heights where i + j is odd and velocities where it is even, so every neighbor of a height point is a
# velocity point and vice versa. Grid A's scheme on this lattice, whose differences span d* = 2s.
function tendencies!(du, dv, dh, u, v, h, p::Problem, grid::EGrid)
    f, g, H = p.f, gravity(p), p.H
    n, Δ = lattice_size(p, grid), spacing(p, grid)

    @inbounds for j in 1:n, i in 1:n
        (i₊, i₋), (j₊, j₋) = neighbors(i, n), neighbors(j, n)
        if is_height_point(grid, i, j)
            du[i, j] = dv[i, j] = 0.0
            dh[i, j] = -H * ((u[i₊, j] - u[i₋, j]) + (v[i, j₊] - v[i, j₋])) / 2Δ
        else
            du[i, j] = f * v[i, j] - g * (h[i₊, j] - h[i₋, j]) / 2Δ
            dv[i, j] = -f * u[i, j] - g * (h[i, j₊] - h[i, j₋]) / 2Δ
            dh[i, j] = 0.0
        end
    end

    return nothing
end

#####
##### Time stepping
#####

struct State
    u::Matrix{Float64}
    v::Matrix{Float64}
    h::Matrix{Float64}
end

State(n::Int) = State(zeros(n, n), zeros(n, n), zeros(n, n))

tendencies!(tendency::State, state::State, p::Problem, grid) =
    tendencies!(tendency.u, tendency.v, tendency.h, state.u, state.v, state.h, p, grid)

"""
    initial_state(p, grid)

The mound, exp(-r²/2σ²), at the height points of `grid`, with the fluid at rest.
"""
function initial_state(p::Problem, grid)
    n = lattice_size(p, grid)
    x, y = coordinates(p, grid)
    state = State(n)

    for j in 1:n, i in 1:n
        if is_height_point(grid, i, j)
            state.h[i, j] = exp(-(x[i]^2 + y[j]^2) / (2 * p.σ^2))
        end
    end

    return state
end

"""
    Stepper(p, grid)

The work arrays of the classical fourth-order Runge–Kutta scheme on `grid`.
"""
struct Stepper{G <: ArakawaGrid}
    problem::Problem
    grid::G
    stage::State
    k::NTuple{4, State}
end

function Stepper(p::Problem, grid)
    n = lattice_size(p, grid)
    return Stepper(p, grid, State(n), ntuple(_ -> State(n), 4))
end

# The state of a Runge–Kutta stage, state + a k.
function stage!(stage::State, state::State, a, k::State)
    @. stage.u = state.u + a * k.u
    @. stage.v = state.v + a * k.v
    @. stage.h = state.h + a * k.h
    return stage
end

"""
    rk4_step!(state, stepper, Δt)

Advance `state` by one step of the classical Runge–Kutta scheme.
"""
function rk4_step!(state::State, stepper::Stepper, Δt)
    (; problem, grid, stage) = stepper
    k₁, k₂, k₃, k₄ = stepper.k

    tendencies!(k₁, state, problem, grid)
    tendencies!(k₂, stage!(stage, state, Δt / 2, k₁), problem, grid)
    tendencies!(k₃, stage!(stage, state, Δt / 2, k₂), problem, grid)
    tendencies!(k₄, stage!(stage, state, Δt, k₃), problem, grid)

    @. state.u += Δt / 6 * (k₁.u + 2k₂.u + 2k₃.u + k₄.u)
    @. state.v += Δt / 6 * (k₁.v + 2k₂.v + 2k₃.v + k₄.v)
    @. state.h += Δt / 6 * (k₁.h + 2k₂.h + 2k₃.h + k₄.h)

    return state
end

"""
    advance!(state, stepper, Δt, steps)

Advance `state` by `steps` Runge–Kutta steps of `Δt`.
"""
function advance!(state::State, stepper::Stepper, Δt, steps)
    for _ in 1:steps
        rk4_step!(state, stepper, Δt)
    end

    return state
end

#####
##### Diagnostics
#####

"""
    energy(state, p, grid)

The semi-discrete energy, Σ [H(u² + v²) + gh²]/2 times the area per point, which every scheme conserves.
"""
function energy(state::State, p::Problem, grid)
    g, H = gravity(p), p.H
    return area_per_point(p, grid) / 2 * (H * (sum(abs2, state.u) + sum(abs2, state.v)) + g * sum(abs2, state.h))
end

total_mass(state::State, p::Problem, grid) = sum(state.h) * area_per_point(p, grid)

#####
##### Inertia–gravity waves, and the exact solutions for a mound released from rest
#####

# The Coriolis factor α and the squared effective wavenumber p² of each grid, in terms of the half-angles a = kΔ/2 and
# b = lΔ/2, where Δ is the lattice spacing. Grid E's differences span 2Δ, as grid A's do.
coriolis_factor(::ArakawaGrid, a, b) = 1.0
coriolis_factor(::Union{CGrid, DGrid}, a, b) = cos(a) * cos(b)

effective_wavenumber²(::Union{AGrid, EGrid}, a, b, Δ) = (sin(2a)^2 + sin(2b)^2) / Δ^2
effective_wavenumber²(::BGrid, a, b, Δ) = 4 * (sin(a)^2 * cos(b)^2 + cos(a)^2 * sin(b)^2) / Δ^2
effective_wavenumber²(::CGrid, a, b, Δ) = 4 * (sin(a)^2 + sin(b)^2) / Δ^2
effective_wavenumber²(::DGrid, a, b, Δ) = (sin(2a)^2 * cos(b)^2 + sin(2b)^2 * cos(a)^2) / Δ^2

"""
    fourier_symbols(grid, k, l, p)

The Coriolis factor α and the squared effective wavenumber p² of `grid` for the wavenumbers (k, l), so that its
inertia–gravity waves have the frequency ν² = α²f² + gHp² (Cushman-Roisin & Beckers 2011, eq. 9.66 and Table 9.1, with
grid E added). The continuous equations, `Exact()`, have α = 1 and p² = k² + l².
"""
function fourier_symbols(grid::ArakawaGrid, k, l, p::Problem)
    Δ = spacing(p, grid)
    a, b = k * Δ / 2, l * Δ / 2
    return coriolis_factor(grid, a, b), effective_wavenumber²(grid, a, b, Δ)
end

fourier_symbols(::Exact, k, l, p::Problem) = 1.0, k^2 + l^2

"""
    frequency(grid, k, l, p)

The frequency ν of the inertia–gravity waves with wavenumbers (k, l) on `grid`, or in the continuous equations if
`grid = Exact()`.
"""
function frequency(grid, k, l, p::Problem)
    α, p² = fourier_symbols(grid, k, l, p)
    return sqrt(α^2 * p.f^2 + gravity(p) * p.H * p²)
end

"""
    adjusted(ĥ₀, grid, k, l, t, p)

The Fourier amplitude at time t of a height disturbance ĥ₀ with wavenumbers (k, l), released from rest:
ĥ(t) = ĥ₀ [α²f²/ν² + (1 - α²f²/ν²) cos νt]. The first term is the geostrophically balanced part that stays behind; the
second radiates away as inertia–gravity waves.
"""
function adjusted(ĥ₀, grid, k, l, t, p::Problem)
    α, p² = fourier_symbols(grid, k, l, p)
    ν² = α^2 * p.f^2 + gravity(p) * p.H * p²
    ν² == 0 && return ĥ₀   # grid D at kd = ld = π: a stationary wave

    balanced = α^2 * p.f^2 / ν²
    return ĥ₀ * (balanced + (1 - balanced) * cos(sqrt(ν²) * t))
end

# The wavenumbers of the discrete Fourier transform on n points over a length L, in FFTW's order.
wavenumbers(n, L) = 2π / L .* fftfreq(n, n)

"""
    adjust(h₀, grid, t, p)

The height at time t that starts from the height `h₀` at rest, on a periodic square lattice that covers the domain,
from the `adjusted` amplitudes of its Fourier modes.
"""
function adjust(h₀, grid, t, p::Problem)
    n = size(h₀, 1)
    k = wavenumbers(n, domain_length(p))
    ĥ = fft(h₀)

    for j in 1:n, i in 1:n
        ĥ[i, j] = adjusted(ĥ[i, j], grid, k[i], k[j], t, p)
    end

    return real(ifft(ĥ))
end

"""
    semidiscrete_h(p, grid, t)

The exact solution of the semi-discrete equations of `grid` at time t. On grid E, the scheme is grid A's on the lattice
of spacing s, applied to heights that start on alternate points only. The other half of the lattice, with heights on
the even points, then stays at rest (A&L, p. 183, note the equivalence in one dimension).
"""
semidiscrete_h(p::Problem, grid, t) = adjust(initial_state(p, grid).h, grid, t, p)

"""
    exact_h(p, t; refine = 4)

The exact solution of the continuous equations at time t, on a periodic grid `refine` times finer than d, centered on
the mound. Returns (x, y, h).
"""
function exact_h(p::Problem, t; refine = 4)
    L = domain_length(p)
    n = refine * p.N
    x = [mod(i * L / n + L / 2, L) - L / 2 for i in 0:n-1]   # x = 0 at index 1
    h₀ = [exp(-(xᵢ^2 + yⱼ^2) / (2 * p.σ^2)) for xᵢ in x, yⱼ in x]
    h = adjust(h₀, Exact(), t, p)

    order = sortperm(x)
    return x[order], x[order], h[order, order]
end
