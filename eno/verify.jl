# Checks of eno.jl against published values, the numbers in the Wikipedia article, and an earlier independent
# implementation in Python (the prototypes behind the figure plan). Prints each check and stops at the first failure.
#
#     julia --project verify.jl

include("eno.jl")

using Printf
using Random

function check(description, passed)
    println(passed ? "  pass  " : "  FAIL  ", description)
    passed || error("check failed: $description")
    return nothing
end

relative_error(a, b) = abs(a - b) / abs(b)
orders(errors) = [log2(errors[n] / errors[n + 1]) for n in 1:length(errors)-1]

# Numbers as a list, each printed with the Printf format `format`.
listed(numbers, format) = join((Printf.format(Printf.Format(format), x) for x in numbers), ", ")

#####
##### Reconstruction
#####

println("Coefficients c_{r,j}")

# Shu (1998), Table 2.1, the rows for k = 3 and 4.
shu_table = Dict(
        (3, -1) => [11//6, -7//6, 1//3], (3, 0) => [1//3, 5//6, -1//6],
        (3, 1) => [-1//6, 5//6, 1//3], (3, 2) => [1//3, -7//6, 11//6],
        (4, -1) => [25//12, -23//12, 13//12, -1//4], (4, 0) => [1//4, 13//12, -5//12, 1//12],
        (4, 1) => [-1//12, 7//12, 7//12, -1//12], (4, 2) => [1//12, -5//12, 13//12, 1//4],
        (4, 3) => [-1//4, 13//12, -23//12, 25//12]
)

check("k = 3 and 4 match Shu (1998), Table 2.1",
      all([eno_coefficient(k, r, j) for j in 0:k-1] == row for ((k, r), row) in shu_table))
check("every row sums to 1, for k ≤ 7",
      all(sum(eno_coefficient(k, r, j) for j in 0:k-1) == 1 for k in 1:7 for r in -1:k-1))
check("rows are mirror images: c_{r,j} = c_{k-2-r,k-1-j}, for k ≤ 7",
      all(eno_coefficient(k, r, j) == eno_coefficient(k, k - 2 - r, k - 1 - j)
          for k in 1:7 for r in -1:k-1 for j in 0:k-1))

"""Whether the polynomial from the stencil of cell `i` with shift `r` has the averages `v` and the right face values."""
function reconstruction_agrees(v, i, k, r)
    a = reconstruction_polynomial(v, i, k, r)
    offsets = -r:k-1-r
    averages = [sum(a[n + 1] * ((m + 1)^(n + 1) - m^(n + 1)) / (n + 1) for n in 0:k-1) for m in offsets]

    return isapprox(averages, [v(i + m) for m in offsets]; rtol = 1e-12) &&
           isapprox(evalpoly(1.0, a), right_face_value(v, i, k, r); rtol = 1e-12) &&
           isapprox(evalpoly(0.0, a), left_face_value(v, i, k, r); rtol = 1e-12)
end

data = [0.7, -0.4, 1.3, 2.2, 0.1, -1.0, 0.9]
check("the reconstruction polynomial has the given averages and the tabulated face values",
      all(reconstruction_agrees(i -> data[i], 4, k, r) for k in 1:4 for r in 0:k-1))

println("\nThe worked example (worked_example.jl): k = 3, averages 1.0, 1.2, 1.3, 0.0, 0.1, 0.3 in cells 1 to 6")

worked = [1.0, 1.2, 1.3, 0.0, 0.1, 0.3]
w = i -> worked[i]
check("cell 3 chooses r = 2 and cell 4 r = 0", eno_shift(w, 3, 3) == 2 && eno_shift(w, 4, 3) == 0)
check("ENO face values 1.267 and 1.317 (cell 3), -0.017 and 0.033 (cell 4)",
      isapprox(left_face_value(w, 3, 3, 2), 19 / 15) && isapprox(right_face_value(w, 3, 3, 2), 79 / 60) &&
      isapprox(left_face_value(w, 4, 3, 0), -1 / 60) && isapprox(right_face_value(w, 4, 3, 0), 1 / 30))
check("fixed stencil r = 1: 1.483 and 0.883 (cell 3), 0.417 and -0.183 (cell 4)",
      isapprox(left_face_value(w, 3, 3, 1), 89 / 60) && isapprox(right_face_value(w, 3, 3, 1), 53 / 60) &&
      isapprox(left_face_value(w, 4, 3, 1), 5 / 12) && isapprox(right_face_value(w, 4, 3, 1), -11 / 60))

step_data = i -> i ≤ 0 ? 1.0 : 0.0
check("step data: the central stencil overshoots to 7/6 and -1/6",
      isapprox(left_face_value(step_data, 0, 3, 1), 7 / 6) && isapprox(right_face_value(step_data, 1, 3, 1), -1 / 6))

#####
##### Advection
#####

println("\nConvergence on smooth data: sin(2πx) advected once round [0, 1]")

# Where the stencil changes orientation, the coefficient of the error jumps, so the scheme is of order k only in the
# sense of cell averages and its maximum error may converge an order more slowly (Harten et al. 1997, Remarks 1.2 and
# 3.1). For k = 2 this happens at every extremum, where ENO, being TVD, is only first-order accurate.
for k in 2:4
    results = map((40, 80, 160, 320)) do N
        Δx = 1 / N
        exact = [cell_average(x -> sin(2π * x), (j - 1) * Δx, j * Δx) for j in 1:N]
        steps, Δt = time_steps(1.0, 0.5Δx * Δx^((k - 1) / 3))   # keeps the time error below the space error
        u = advect(exact, ENO(k); Δx, Δt, steps)[end]
        sum(abs.(u .- exact)) / N, maximum(abs.(u .- exact))
    end
    l1, l∞ = first.(results), last.(results)
    @printf("    ENO-%d: L1 orders %s; maximum-error orders %s\n", k, listed(orders(l1), "%.2f"),
            listed(orders(l∞), "%.2f"))
    check("ENO-$k converges at order $k in the L1 norm and at least order $(k - 1) in the maximum norm",
          abs(orders(l1)[end] - k) < 0.2 && orders(l∞)[end] > k - 1.1)
end

println("\nLoss of accuracy on sin⁴(πx): third-order ENO, period 2, t = 4, Δt = 0.6Δx, as in Shu (1990)")

sin4_errors = map((160, 320, 640, 1280, 2560)) do N
    Δx = 2 / N
    x = [-1 + (j - 0.5) * Δx for j in 1:N]
    u₀ = sinpi.(x) .^ 4
    steps, Δt = time_steps(4.0, 0.6Δx)
    u = advect(u₀, ENO(3); Δx, Δt, steps)[end]
    sum(abs.(u .- u₀)) / N
end
observed = orders(sin4_errors)
@printf("    L1 errors %s; orders %s\n", listed(sin4_errors, "%.2e"), listed(observed, "%.2f"))
check("the order at the finest grids is near 2, not 3 (Shu 1990 reports about 1.97)", abs(observed[end] - 2) < 0.4)

println("\nJiang & Shu's profile, 200 point values, t = 8, Δt = 0.4Δx, as in the Python prototypes (E3)")

N = 200
Δx = 2 / N
x = [-1 + (j - 0.5) * Δx for j in 1:N]
u₀ = jiang_shu_profile.(x)
steps, Δt = time_steps(8.0, 0.4Δx)
total_variation(u) = sum(abs(u[mod1(j + 1, length(u))] - u[j]) for j in eachindex(u))

python = [
        ("First-order upwind", FirstOrderUpwind(), (0.3707, 0.1150, 0.511, 0.3145)),
        ("MUSCL, minmod", Minmod(), (0.6829, 0.0012, 2.544, 0.2100)),
        ("ENO-3", ENO(3), (0.9991, -0.0002, 6.331, 0.0882))
]

for (name, scheme, expected) in python
    u = advect(u₀, scheme; Δx, Δt, steps)[end]
    found = (maximum(u), minimum(u), total_variation(u), sum(abs.(u .- u₀)) / N)
    @printf("    %-20s max %.4f, min %.4f, total variation %.3f, L1 error %.4f\n", name, found...)
    check("$name matches the Python values to the digits given",
          all(abs(f - e) ≤ 0.6 * 10.0^(-d) for (f, e, d) in zip(found, expected, (4, 4, 3, 4))))
end

println("\nFifth-order WENO of Jiang & Shu (1996), kept for the WENO article")

weights_reproduce_fixed = let combined = zeros(5)
    for (w, r) in zip(WENO5_LINEAR_WEIGHTS, (2, 1, 0))
        combined[(3 - r):(5 - r)] .+= w .* coefficients(3, r)   # cells i-r, ..., i-r+2 of the five cells i-2, ..., i+2
    end
    isapprox(combined, coefficients(5, 2); rtol = 1e-14)
end
check("the linear weights 1/10, 6/10 and 3/10 combine the three stencils into the fixed five-cell stencil, r = 2",
      weights_reproduce_fixed)

# On smooth data the weights differ from the linear ones by O(h²), and the face values by O(h⁵) (Jiang & Shu 1996, 3.9).
departure(a) = (f = i -> sinpi(a * i); abs(face_value(WENO5(), f, 3) - right_face_value(f, 3, 5, 2)))
@printf("    departure from the fixed stencil on sin(πai): %.2e for a = 0.02, %.2e for a = 0.01\n", departure(0.02),
        departure(0.01))
check("on smooth data WENO approaches the fixed five-cell stencil at fifth order",
      abs(log2(departure(0.02) / departure(0.01)) - 5) < 0.3)

# Borges et al. (2008), Table 4: WENO-JS on this profile at t = 8, with Δt ≤ Δx^(5/3) and the L1 error Δx Σ|e|, on N
# points that include x = -1. Their errors are 2–3% above ours at every N; a wrong weight or indicator would differ by
# more.
function borges_error(N)
    Δx = 2 / N
    u₀ = [jiang_shu_profile(-1 + (j - 1) * Δx) for j in 1:N]
    steps, Δt = time_steps(8.0, min(0.4Δx, Δx^(5 / 3)))
    return Δx * sum(abs.(advect(u₀, WENO5(); Δx, Δt, steps)[end] .- u₀))
end

for (N, published) in ((50, 4.64e-1), (100, 2.11e-1), (200, 9.32e-2), (400, 4.01e-2))
    l1_error = borges_error(N)
    @printf("    N = %3d: L1 error %.3e (Borges et al.: %.2e)\n", N, l1_error, published)
    check("WENO-5 matches Borges et al. (2008), Table 4, to 4% at N = $N", relative_error(l1_error, published) < 0.04)
end

#####
##### Burgers' equation
#####

println("\nBurgers' equation from u₀ = 1/2 + sin(πx) (burgers_shock.jl)")

# The exact solution is carried along characteristics, so u = u₀(x − ut) wherever it is smooth, before the shock and
# after; across the shock u drops (the entropy condition), and the mean of the two sides is the shock speed, 1/2.
characteristic_residual(t) = maximum(abs(u - (0.5 + sinpi(x - u * t))) for x in range(-0.999, 0.999; length = 999)
                                     for u in (burgers_exact(x, t),) if abs(x - burgers_shock_position(t)) > 1e-6)
@printf("    largest residual of u = u₀(x − ut): %.1e at t = 0.2, %.1e at t = 1.5\n", characteristic_residual(0.2),
        characteristic_residual(1.5))
check("the exact solution satisfies u = u₀(x − ut) before and after the shock forms",
      characteristic_residual(0.2) < 1e-12 && characteristic_residual(1.5) < 1e-12)

uₗ, uᵣ = burgers_exact(burgers_shock_position(1.5) - 1e-9, 1.5), burgers_exact(burgers_shock_position(1.5) + 1e-9, 1.5)
@printf("    at t = 1.5 the shock drops from %.4f to %.4f\n", uₗ, uᵣ)
check("the shock satisfies the entropy condition and moves at the mean of its two sides, 1/2",
      uₗ > uᵣ && abs((uₗ + uᵣ) / 2 - 0.5) < 1e-9)

function burgers_averages(N, t)
    edges = range(-1, 1; length = N + 1)
    jumps = t > 1 / π ? (burgers_shock_position(t),) : ()
    return [cell_average(x -> burgers_exact(x, t), edges[j], edges[j + 1]; jumps) for j in 1:N]
end

check("the exact solution conserves the integral of u after the shock",
      abs(sum(burgers_averages(80, 1.5)) / 80 - 0.5) < 1e-9)

function burgers_error(scheme, N, T)
    Δx = 2 / N
    steps, Δt = time_steps(T, 0.4Δx / 1.5)
    u = evolve(burgers_averages(N, 0.0), (du, u) -> burgers_tendency!(du, u, scheme, Δx); Δt, steps)[end]
    return Δx * sum(abs.(u .- burgers_averages(N, T)))
end

# Before the shock the fixed stencil converges at third order, which checks the exact solution and the flux. ENO-3,
# with the same flux, reaches only about 2.7 here and nearer 2 on finer grids: the loss of accuracy from switching
# stencils on smooth data (Rogerson & Meiburg 1990), so it is held only to better than second order.
for (name, scheme, required) in (("fixed stencil", FixedStencil(3, 1), 2.9), ("ENO-3", ENO(3), 2.0))
    errors = [burgers_error(scheme, N, 0.15) for N in (80, 160, 320)]
    @printf("    %-13s L1 errors at t = 0.15 on 80, 160 and 320 cells: %s, orders %s\n", name, listed(errors, "%.2e"),
            listed(orders(errors), "%.2f"))
    check("$name converges at an order above $required before the shock", all(orders(errors) .> required))
end

#####
##### Solid-body rotation
#####

println("\nThe rotating slotted disk (slotted_disk.jl)")

# Mirroring the data in x reverses the rotation, so the tendency for -ω on the mirrored data must be the mirrored
# tendency for ω: this exercises the grid lines with negative speed, which are read backwards.
let N = 24, xs = [-1 + (i - 0.5) * (2 / N) for i in 1:N], u = rand(Xoshiro(3), N, N)
    for scheme in (FixedStencil(3, 1), ENO(3))
        forward = rotation_tendency!(zeros(N, N), u, scheme, xs, π, 2 / N)
        mirrored = rotation_tendency!(zeros(N, N), u[end:-1:1, :], scheme, xs, -π, 2 / N)
        check("mirroring the data and the rotation mirrors the tendency, for $(nameof(typeof(scheme)))",
              maximum(abs.(mirrored[end:-1:1, :] .- forward)) < 1e-12)
    end
end

function rotation_error(scheme, N)
    Δx = 2 / N
    xs = [-1 + (i - 0.5) * Δx for i in 1:N]
    bump(x, y) = exp(-(x^2 + (y - 0.5)^2) / 0.03)
    u₀ = [bump(x, y) for x in xs, y in xs]
    steps, Δt = time_steps(2.0, 0.4Δx / (π * 2maximum(abs, xs)))
    u = evolve(u₀, (du, u) -> rotation_tendency!(du, u, scheme, xs, π, Δx); Δt, steps)[end]
    return sum(abs.(u .- u₀)) / N^2
end

# The fixed stencil approaches third order (2.53, then 2.89 on 160² cells), which checks the line-by-line scheme;
# ENO-3 stays near 2.5 (2.67, then 2.54), the same loss of accuracy from switching stencils as on Burgers' equation.
for (scheme, sizes, required) in ((FixedStencil(3, 1), (80, 160), 2.8), (ENO(3), (40, 80), 2.0))
    errors = [rotation_error(scheme, N) for N in sizes]
    @printf("    %-12s mean error after one turn of a smooth bump on %d² and %d² cells: %.2e, %.2e, order %.2f\n",
            nameof(typeof(scheme)), sizes..., errors..., orders(errors)[1])
    check("$(nameof(typeof(scheme))) turns a smooth bump at an order above $required", orders(errors)[1] > required)
end

println("\nAll checks passed.")
