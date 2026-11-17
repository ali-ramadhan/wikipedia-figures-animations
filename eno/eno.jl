# The ENO reconstruction, the schemes it is compared with, and the test problems, shared by the figures and animations.
# The notation follows the Wikipedia article and Shu (1998), section 2: k is the order of accuracy, a stencil of k cells
# reaches r cells to the left of its target cell, and the coefficients c_{r,j} turn its k cell averages into the value
# at the target cell's right face. Cell averages are passed as a function `v` of the cell index, such as `periodic(u)`
# for an array `u` on a periodic grid.

using LinearAlgebra

#####
##### Reconstruction from cell averages
#####

"""
    eno_coefficient(k, r, j)

The coefficient ``c_{r,j}`` of a reconstruction of order `k` on a uniform grid, in exact arithmetic: the value at the
right face of cell ``i`` is ``\\sum_{j=0}^{k-1} c_{r,j} \\bar{u}_{i-r+j}``. This is eq. (2.21) of Shu (1998).
"""
function eno_coefficient(k, r, j)
    c = 0 // 1

    for m in j+1:k
        numerator = sum(prod((r - q + 1 for q in 0:k if q != m && q != l); init = 1) for l in 0:k if l != m)
        denominator = prod(m - l for l in 0:k if l != m)
        c += numerator // denominator
    end

    return c
end

# COEFFICIENTS[k][r + 2] holds c_{r,0}, ..., c_{r,k-1}, for r = -1, ..., k - 1.
const COEFFICIENTS = [[Float64[eno_coefficient(k, r, j) for j in 0:k-1] for r in -1:k-1] for k in 1:7]

coefficients(k, r) = COEFFICIENTS[k][r + 2]

"""
    undivided_difference(v, j, l)

The undivided difference ``D^l_j`` of the cell averages `v`, starting at cell `j`: ``D^0_j = v_j`` and
``D^l_j = D^{l-1}_{j+1} - D^{l-1}_j``.
"""
undivided_difference(v, j, l) = sum((-1)^(l - m) * binomial(l, m) * v(j + m) for m in 0:l)

"""
    eno_shift(v, i, k)

The shift `r` of the ENO stencil of cell `i`: how many cells it reaches to the left. Starting from cell `i` alone, the
stencil gains one neighbour at a time, on whichever side gives the smaller undivided difference, and on the right in a
tie.
"""
function eno_shift(v, i, k)
    leftmost = i

    for l in 1:k-1
        left = abs(undivided_difference(v, leftmost - 1, l))
        right = abs(undivided_difference(v, leftmost, l))
        leftmost -= left < right
    end

    return i - leftmost
end

"""The value ``u^-_{i+1/2}`` at the right face of cell `i`, from the stencil of `k` cells with shift `r`."""
right_face_value(v, i, k, r) = sum(coefficients(k, r)[j + 1] * v(i - r + j) for j in 0:k-1)

"""The value ``u^+_{i-1/2}`` at the left face of cell `i`, which uses the row ``r - 1`` of the coefficients."""
left_face_value(v, i, k, r) = sum(coefficients(k, r - 1)[j + 1] * v(i - r + j) for j in 0:k-1)

"""
    reconstruction_polynomial(v, i, k, r)

The coefficients ``a_0, …, a_{k-1}`` of the polynomial ``p(ξ) = \\sum_n a_n ξ^n`` whose averages over the `k` cells of
the stencil with shift `r` equal the cell averages `v`. Here ``ξ = (x - x_{i-1/2})/Δx``, so cell `i` is ``0 ≤ ξ ≤ 1``.
Evaluate it with `evalpoly(ξ, a)`.
"""
function reconstruction_polynomial(v, i, k, r)
    offsets = -r:k-1-r
    A = [((m + 1)^(n + 1) - m^(n + 1)) / (n + 1) for m in offsets, n in 0:k-1]
    return A \ [v(i + m) for m in offsets]
end

"""Index a vector periodically: `periodic(u)(i)` is `u[mod1(i, length(u))]`."""
periodic(u) = i -> @inbounds u[mod1(i, length(u))]

#####
##### Schemes
#####

abstract type Scheme end

"""First-order upwind: the value at a face is the cell average on its upwind side."""
struct FirstOrderUpwind <: Scheme end

"""A linear scheme of order `k` with the same stencil everywhere, reaching `r` cells to the left."""
struct FixedStencil <: Scheme
    k :: Int
    r :: Int
end

"""The second-order TVD scheme with the minmod limiter (MUSCL)."""
struct Minmod <: Scheme end

"""The second-order TVD scheme with the superbee limiter (MUSCL), the most compressive limiter that keeps it TVD."""
struct Superbee <: Scheme end

"""ENO of order `k`."""
struct ENO <: Scheme
    k :: Int
end

"""
The fifth-order WENO scheme of Jiang & Shu (1996): a convex combination of the face values from the three candidate
stencils of third-order ENO, with weights from their smoothness indicators. Where the data are smooth the weights are
close to 1/10, 6/10 and 3/10, which give the fixed five-cell stencil with r = 2. No figure of the ENO article uses it;
it is kept, with its checks in verify.jl, for the WENO article.
"""
struct WENO5 <: Scheme
    ε :: Float64
end

WENO5() = WENO5(1e-6)

minmod(a, b) = a * b > 0 ? sign(a) * min(abs(a), abs(b)) : zero(a)
superbee(a, b) = a * b > 0 ? sign(a) * max(min(2abs(a), abs(b)), min(abs(a), 2abs(b))) : zero(a)

"""The smoothness indicators of Jiang & Shu (1996), eqs. (3.2)–(3.4), for the stencils with r = 2, 1 and 0."""
function smoothness_indicators(v, i)
    return (13 / 12 * (v(i - 2) - 2v(i - 1) + v(i))^2 + 1 / 4 * (v(i - 2) - 4v(i - 1) + 3v(i))^2,
            13 / 12 * (v(i - 1) - 2v(i) + v(i + 1))^2 + 1 / 4 * (v(i - 1) - v(i + 1))^2,
            13 / 12 * (v(i) - 2v(i + 1) + v(i + 2))^2 + 1 / 4 * (3v(i) - 4v(i + 1) + v(i + 2))^2)
end

const WENO5_LINEAR_WEIGHTS = (1 / 10, 6 / 10, 3 / 10)   # for r = 2, 1 and 0

"""The shift `r` of the stencil that a fixed stencil, or ENO, uses for cell `i`."""
stencil_shift(scheme::FixedStencil, v, i) = scheme.r
stencil_shift(scheme::ENO, v, i) = eno_shift(v, i, scheme.k)

# The value at the right face of cell i, for flow to the right.
face_value(::FirstOrderUpwind, v, i) = v(i)
face_value(::Minmod, v, i) = v(i) + minmod(v(i) - v(i - 1), v(i + 1) - v(i)) / 2
face_value(::Superbee, v, i) = v(i) + superbee(v(i) - v(i - 1), v(i + 1) - v(i)) / 2
face_value(scheme::Union{FixedStencil, ENO}, v, i) = right_face_value(v, i, scheme.k, stencil_shift(scheme, v, i))

function face_value(scheme::WENO5, v, i)
    α = WENO5_LINEAR_WEIGHTS ./ (scheme.ε .+ smoothness_indicators(v, i)) .^ 2
    candidates = (right_face_value(v, i, 3, 2), right_face_value(v, i, 3, 1), right_face_value(v, i, 3, 0))
    return sum(α .* candidates) / sum(α)
end

"""
    conservative_tendency!(du, face_flux, Δx)

The semi-discrete conservative scheme ``du_i/dt = -(F_{i+1/2} - F_{i-1/2})/Δx`` on a periodic grid, where
`face_flux(i)` is the flux ``F_{i+1/2}`` through the right face of cell `i`.
"""
function conservative_tendency!(du, face_flux, Δx)
    previous = face_flux(0)

    for i in eachindex(du)
        next = face_flux(i)
        du[i] = -(next - previous) / Δx
        previous = next
    end

    return du
end

"""
    advection_tendency!(du, u, scheme, Δx)

The semi-discrete scheme for ``u_t + u_x = 0`` on a periodic grid, whose flux through each face is the scheme's value
there. Read as cell averages, `u` gives the finite-volume scheme; read as point values, the finite-difference scheme,
which for this equation is the same.
"""
function advection_tendency!(du, u, scheme, Δx)
    v = periodic(u)
    return conservative_tendency!(du, i -> face_value(scheme, v, i), Δx)
end

#####
##### Time stepping
#####

"""
    ssp_rk3!(u, tendency!, Δt, scratch)

One step of the third-order strong-stability-preserving (TVD) Runge–Kutta method of Shu & Osher (1988), with three
arrays like `u` as `scratch`.
"""
function ssp_rk3!(u, tendency!, Δt, (u₁, u₂, du))
    tendency!(du, u)
    @. u₁ = u + Δt * du

    tendency!(du, u₁)
    @. u₂ = 3 / 4 * u + 1 / 4 * (u₁ + Δt * du)

    tendency!(du, u₂)
    @. u = 1 / 3 * u + 2 / 3 * (u₂ + Δt * du)

    return u
end

"""The number of equal steps, and their length, that reach time `T` with steps no longer than `Δt_max`."""
function time_steps(T, Δt_max)
    n = ceil(Int, T / Δt_max - 1e-9)
    return n, T / n
end

"""
    evolve(u₀, tendency!; Δt, steps, save_every = steps, callback = nothing)

Step ``du/dt``, which `tendency!(du, u)` sets, from `u₀` with `ssp_rk3!` for `steps` steps of length `Δt`, and return
the solution at the start, after every `save_every` steps, and after the last step. `callback(u, n)` is called after
each number of steps `n`, from 0 to `steps`.
"""
function evolve(u₀, tendency!; Δt, steps, save_every = steps, callback = nothing)
    u = copy(u₀)
    scratch = (similar(u), similar(u), similar(u))
    snapshots = [copy(u)]

    for n in 1:steps
        isnothing(callback) || callback(u, n - 1)
        ssp_rk3!(u, tendency!, Δt, scratch)
        (n % save_every == 0 || n == steps) && push!(snapshots, copy(u))
    end

    isnothing(callback) || callback(u, steps)

    return snapshots
end

"""
    advect(u₀, scheme; Δx, Δt, steps, save_every = steps, callback = nothing)

Advect the periodic data `u₀` with speed 1, using `scheme`. The other keywords are those of `evolve`.
"""
advect(u₀, scheme; Δx, kwargs...) = evolve(u₀, (du, u) -> advection_tendency!(du, u, scheme, Δx); kwargs...)

#####
##### Initial data
#####

const GAUSS_NODES = (-0.9061798459386640, -0.5384693101056831, 0.0, 0.5384693101056831, 0.9061798459386640)
const GAUSS_WEIGHTS = (0.2369268850561891, 0.4786286704993665, 0.5688888888888889, 0.4786286704993665,
                       0.2369268850561891)

"""
    cell_average(f, a, b; jumps = ())

The average of `f` over ``[a, b]``, by five-point Gauss–Legendre quadrature on each piece between the points `jumps`
where `f` is discontinuous: exact to round-off for the smooth pieces used here.
"""
function cell_average(f, a, b; jumps = ())
    cuts = [a; sort([x for x in jumps if a < x < b]); b]
    total = 0.0

    for (l, r) in zip(cuts[1:end-1], cuts[2:end]), (ξ, w) in zip(GAUSS_NODES, GAUSS_WEIGHTS)
        total += w * (r - l) / 2 * f((l + r) / 2 + ξ * (r - l) / 2)
    end

    return total / (b - a)
end

"""
    jiang_shu_profile(x)

The initial condition of Jiang & Shu (1996), Example 1, on ``-1 ≤ x ≤ 1``: a smooth narrow pulse made of three
Gaussians, a square wave, a triangle and half an ellipse.
"""
function jiang_shu_profile(x)
    a, z, δ, α = 0.5, -0.7, 0.005, 10.0
    β = log(2) / (36δ^2)
    G(x, z) = exp(-β * (x - z)^2)
    F(x, a) = sqrt(max(1 - α^2 * (x - a)^2, 0))

    -0.8 ≤ x ≤ -0.6 && return (G(x, z - δ) + G(x, z + δ) + 4G(x, z)) / 6
    -0.4 ≤ x ≤ -0.2 && return 1.0
    0.0 ≤ x ≤ 0.2 && return 1 - abs(10 * (x - 0.1))
    0.4 ≤ x ≤ 0.6 && return (F(x, a - δ) + F(x, a + δ) + 4F(x, a)) / 6
    return 0.0
end

#####
##### Burgers' equation
#####

"""
    burgers_exact(x, t)

The entropy solution of Burgers' equation ``u_t + (u^2/2)_x = 0`` with ``u_0(x) = 1/2 + \\sin(πx)`` on the periodic
domain ``-1 ≤ x ≤ 1``. In the frame ``ξ = x - t/2`` that moves with the mean speed, ``w = u - 1/2`` starts as
``\\sin(πξ)`` and stays odd in ξ, so the shock that forms at ``t = 1/π`` stays at ``ξ = ±1``. For ``0 < ξ < 1``, w is
carried along the characteristic ``ξ = η + t\\sin(πη)`` from the smallest η that reaches ξ: on the branch where
``η + t\\sin(πη)`` still increases, before the characteristics that have met the shock.
"""
function burgers_exact(x, t)
    ξ = mod(x - t / 2 + 1, 2) - 1
    ξ == 0 && return 0.5
    η_max = t * π > 1 ? acos(-1 / (t * π)) / π : 1.0
    lo, hi = 0.0, η_max

    for _ in 1:60
        η = (lo + hi) / 2
        η + t * sinpi(η) < abs(ξ) ? (lo = η) : (hi = η)
    end

    return 0.5 + sign(ξ) * sinpi((lo + hi) / 2)
end

"""Where the shock of `burgers_exact` is at a time ``t ≥ 1/π``, in ``-1 ≤ x < 1``."""
burgers_shock_position(t) = mod(t / 2 + 2, 2) - 1

"""The Godunov flux for ``f(u) = u^2/2``: the flux of the exact solution of the Riemann problem with states `a`, `b`."""
godunov_burgers(a, b) = max(max(a, 0)^2, min(b, 0)^2) / 2

"""
    burgers_tendency!(du, u, scheme, Δx)

The finite-volume scheme for Burgers' equation on a periodic grid, with a fixed stencil or ENO as `scheme` and the
Godunov flux of the two values at each face: ``u^-_{i+1/2}`` from the stencil of cell ``i`` and ``u^+_{i+1/2}`` from
that of cell ``i+1``.
"""
function burgers_tendency!(du, u, scheme, Δx)
    v, k = periodic(u), scheme.k
    u⁻(i) = right_face_value(v, i, k, stencil_shift(scheme, v, i))
    u⁺(i) = left_face_value(v, i + 1, k, stencil_shift(scheme, v, i + 1))
    return conservative_tendency!(du, i -> godunov_burgers(u⁻(i), u⁺(i)), Δx)
end

#####
##### Solid-body rotation in two dimensions
#####

"""
    line_tendency!(dline, line, speed, scheme, Δx, out)

Add the tendency of ``u_t + c u_s = 0`` along one periodic grid line, with constant speed `speed` = c, to `dline`, with
`out` as scratch. For c < 0 the line is read backwards, so that `advection_tendency!` always sees flow to the right.
"""
function line_tendency!(dline, line, speed, scheme, Δx, out)
    if speed > 0
        advection_tendency!(out, line, scheme, Δx)
        dline .+= speed .* out
    elseif speed < 0
        backwards = reverse(eachindex(line))
        advection_tendency!(out, view(line, backwards), scheme, Δx)
        dline .-= speed .* view(out, backwards)
    end

    return dline
end

"""
    rotation_tendency!(du, u, scheme, xs, ω, Δx)

The tendency of ``u_t + (au)_x + (bu)_y = 0`` with ``(a, b) = ω(-y, x)`` on a square grid with cell centres `xs` in
both directions, where `u[i, j]` belongs to ``(x_i, y_j)``. Along each grid row the speed a is constant, and along each
column b is, so the scheme is applied one grid line at a time, as in the finite-difference form of ENO in several
dimensions.
"""
function rotation_tendency!(du, u, scheme, xs, ω, Δx)
    fill!(du, 0)
    out = similar(xs)   # the tendency along one grid line

    for j in axes(u, 2)
        line_tendency!(view(du, :, j), view(u, :, j), -ω * xs[j], scheme, Δx, out)
    end

    for i in axes(u, 1)
        line_tendency!(view(du, i, :), view(u, i, :), ω * xs[i], scheme, Δx, out)
    end

    return du
end
