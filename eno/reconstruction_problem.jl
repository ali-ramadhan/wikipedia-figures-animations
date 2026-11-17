# The reconstruction problem: from the cell averages of an unknown function, build a polynomial in each cell, which
# then gives two values at every face, one from each side. Written to reconstruction_problem.svg and .png.
#
# Five cells I_{i-2}, ..., I_{i+2} of width 1 have the averages of u(x, t) = 0.9 + 0.5 sin(0.8x + 7π/4)
# + 0.25 sin(1.6x + 2π/3), computed exactly, with the cell I_i centred on x = 0. Each cell's polynomial is the
# third-order reconstruction from the centred stencil (k = 3, r = 1), the quadratic whose averages over the cell and its
# two neighbours equal theirs. This is a fixed stencil, not ENO's choice: for these averages ENO would choose r = 0 in
# I_{i-2} and r = 2 in I_{i+1}, so that I_i and I_{i+1} would share a stencil and the two values at x_{i+1/2} would be
# equal. The grid is coarse so that the two values at each face visibly differ, but each polynomial
# stays within 7% of the range of u. The labelled face x_{i+1/2} is a level shoulder of u, where both polynomials run
# nearly level, so each label fits on its own side of the face.
#
# The text is set in Latin Modern at its 8-point design size (see use_latin_modern! in plotting.jl), drawn as outlines.
#
#     julia --project reconstruction_problem.jl

include("eno.jl")
include("plotting.jl")
use_latin_modern!()

const k, r = 3, 1
const CELLS = -3:3                                         # I_{i-3}, ..., I_{i+3}; the outer two only supply averages
const SHOWN = -2:2                                         # the cells drawn
const WAVES = ((0.5, 0.8, 7π / 4), (0.25, 1.6, 2π / 3))   # (amplitude, wavenumber, phase) of each sine

u(x) = 0.9 + sum(a * sin(κ * x + φ) for (a, κ, φ) in WAVES)

"""The exact average of u over the cell [j - ½, j + ½]."""
average(j) = 0.9 + sum(a * (cos(κ * (j - 0.5) + φ) - cos(κ * (j + 0.5) + φ)) / κ for (a, κ, φ) in WAVES)

@assert all(isapprox(average(j), sum(u(j - 0.5 + (n - 0.5) / 4000) for n in 1:4000) / 4000; atol = 1e-6) for j in CELLS)

polynomials = Dict(j => reconstruction_polynomial(average, j, k, r) for j in SHOWN)
right_value(j) = evalpoly(1.0, polynomials[j])   # u^-_{j+1/2}, at the right face of cell j
left_value(j) = evalpoly(0.0, polynomials[j])    # u^+_{j-1/2}, at the left face of cell j

# The face values agree with the table's coefficients: u^- = -ū_{i-1}/6 + 5ū_i/6 + ū_{i+1}/3, and its mirror image.
for j in SHOWN
    @assert right_value(j) ≈ -average(j - 1) / 6 + 5average(j) / 6 + average(j + 1) / 3
    @assert left_value(j) ≈ average(j - 1) / 3 + 5average(j) / 6 - average(j + 1) / 6
end

jump = right_value(0) - left_value(1)
println("At x_{i+1/2}: u^- = $(round(right_value(0); digits = 4)), u^+ = $(round(left_value(1); digits = 4)), ",
        "jump = $(round(jump; digits = 4))")

"""Draw the polynomial `a`, in the local coordinate ξ of the cell of width 1 that starts at `left`, across that cell."""
function cell_polynomial!(ax, a, left; kwargs...)
    ξ = range(0, 1; length = 40)
    return lines!(ax, left .+ ξ, evalpoly.(ξ, Ref(a)); kwargs...)
end

set_theme!(merge(theme_latexfonts(), ARTICLE_THEME))

fig = Figure(size = (1000, 620), fontsize = 26, figure_padding = (10, 30, 10, 20))
ax = Axis(fig[1, 1]; xticks = (collect(SHOWN), [L"I_{i-2}", L"I_{i-1}", L"I_i", L"I_{i+1}", L"I_{i+2}"]),
          xticklabelsize = 28, xticklabelcolor = INK, limits = (-2.6, 2.6, 0.2, 1.5))
hideydecorations!(ax)
hidespines!(ax, :l)

faces = (first(SHOWN) - 0.5):1:(last(SHOWN) + 0.5)
vlines!(ax, faces; color = FAINT, linewidth = 1)
lines!(ax, [0.5, 0.5], [0.2, 1.38]; color = MUTED, linewidth = 1.5)
text!(ax, 0.5, 1.38; text = L"x_{i+1/2}", align = (:center, :bottom), offset = (0, 6), fontsize = 28)

x = range(-2.5, 2.5; length = 500)
lines!(ax, x, u.(x); color = INK, linewidth = 2.5, linestyle = :dash, label = L"solution $u(x, t)$")

for j in SHOWN
    lines!(ax, [j - 0.5, j + 0.5], [average(j), average(j)]; color = (COLORS.average, 0.55), linewidth = 8,
           label = j == first(SHOWN) ? L"cell averages $\bar{u}_j$" : nothing)
end

for j in SHOWN
    cell_polynomial!(ax, polynomials[j], j - 0.5; color = COLORS.eno, linewidth = 4,
                     label = j == first(SHOWN) ? L"polynomials $p_j(x)$" : nothing)
    scatter!(ax, [j - 0.5, j + 0.5], [left_value(j), right_value(j)]; color = COLORS.eno, markersize = 13)
end

# The notation, at the cell I_i and the face x_{i+1/2}: u^- ends p_i on the left, u^+ starts p_{i+1} on the right.
text!(ax, 0.22, average(0); text = L"\bar{u}_i", align = (:center, :top), offset = (0, -12), fontsize = 28,
      color = MUTED)
text!(ax, -0.15, evalpoly(0.35, polynomials[0]); text = L"p_i(x)", align = (:center, :bottom), offset = (0, 16),
      fontsize = 28, color = COLORS.eno)
text!(ax, 1.2, evalpoly(0.7, polynomials[1]); text = L"p_{i+1}(x)", align = (:center, :bottom), offset = (0, 16),
      fontsize = 28, color = COLORS.eno)
text!(ax, 0.5, right_value(0); text = L"u^-_{i+1/2}", align = (:right, :bottom), offset = (-8, 8), fontsize = 28)
text!(ax, 0.5, left_value(1); text = L"u^+_{i+1/2}", align = (:left, :top), offset = (8, -8), fontsize = 28)

# A white background, without a frame, so that the faces' lines stop at the legend instead of crossing its text.
axislegend(ax; position = :lt, patchsize = (36, 12), rowgap = 4, labelsize = 24, framevisible = true,
           framecolor = :transparent, backgroundcolor = PAPER)

save_figure(fig, "reconstruction_problem")
