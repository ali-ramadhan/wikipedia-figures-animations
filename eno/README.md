# ENO and WENO

Figures and animations for the Wikipedia articles
[Essentially non-oscillatory scheme](https://en.wikipedia.org/wiki/Essentially_non-oscillatory_scheme) and
[WENO methods](https://en.wikipedia.org/wiki/WENO_methods). Oceananigans.jl has WENO schemes but not ENO, so the
reconstruction and the schemes it is compared with are coded here directly, in plain Julia on the CPU. Oceananigans is
in the environment for the WENO figures.

- `eno.jl`: the ENO reconstruction (the coefficients $c_{r,j}$ in exact arithmetic, and the stencil choice), the schemes
  it is compared with (first-order upwind, fixed stencils, and the second-order TVD scheme with the minmod or superbee
  limiter), the third-order strong-stability-preserving Runge–Kutta method, the initial data, Burgers' equation (its
  exact solution and the Godunov flux), and solid-body rotation, applied one grid line at a time. It also has the
  fifth-order WENO scheme of Jiang & Shu (1996), which no figure here uses, for the WENO article.
- `verify.jl`: checks the coefficients against Shu (1998), Table 2.1, the worked example, the order of convergence, the
  loss of accuracy reported by Shu (1990), WENO against the errors of Borges et al. (2008), Table 4, the exact solution
  of Burgers' equation and the schemes' order on it, the rotation, and the numbers an earlier Python implementation gave
  for the same tests (about 3 min).
- `reconstruction_problem.jl`: the reconstruction problem, with the cell averages, a third-order polynomial in each
  cell, and the two values at each face, in Latin Modern at its 8-point design size, which is wider than the 10-point
  design that Makie bundles.
- `reconstruction_notation.jl`: the notation of the article's worked examples on a non-uniform grid (cells, widths,
  averages, the faces in $x$ and in $\xi = x - x_{i+1/2}$, and the two examples' stencils), in Latin Modern 8.
- `fonts/`: Latin Modern Roman 8, regular and italic, from CTAN (`fonts/lm`), under the GUST Font License.
- `coefficient_stencils.jl`: the stencils and coefficients of the article's table of $c_{r,j}$ on a uniform grid, for
  $k = 1$ to 3, in the LaTeX typeface (New Computer Modern); it checks the coefficients from eq. (2.21) against the
  table first.
- `stencil_schematic.jl`: how ENO chooses a third-order stencil for the cell $I_i$ on a non-uniform grid, with each
  candidate stencil labelled by its divided difference $d^l_j = U[x_{j-1/2}, \dots, x_{j+l+1/2}]$ and the procedure beside
  the tree, in Latin Modern 8.
- `worked_example.jl`: a worked example of ENO's stencil choice for the two cells beside a jump, animated, and a still
  in three panels.
- `advection.jl`: twelve schemes advecting the profile of Jiang & Shu (1996), Example 1: first-order upwind,
  second-order centered, and TVD with the minmod and superbee limiters; fixed stencils of orders 2 to 5; and ENO of
  orders 2 to 5. It prints how far each scheme goes outside the range of the exact solution.
- `stencil_flicker.jl`: the stencil ENO chooses at each cell over time, for a smooth pulse with and without noise.
- `burgers_shock.jl`: Burgers' equation from a smooth wave through the formation of a shock, with a fixed third-order
  stencil beside third-order ENO and the exact solution.
- `slotted_disk.jl`: Zalesak's slotted disk turning twice, in a 2 × 2 grid: first-order upwind and second-order
  centered (all numerical diffusion and all numerical dispersion), above a fixed third-order stencil and third-order
  ENO, each applied one grid line at a time. Values more than 0.01 outside $[0, 1]$ are drawn in red or blue. It
  prints the errors after each turn.
- `plotting.jl`: the colors, theme and tick labels shared by the figures and the animations, and `save_figure` and
  `save_animation`, which write them to the output directory.
- `video_encoding.jl`: records the animations without loss, then encodes them for Wikimedia Commons as AV1 in WebM,
  from `arakawa-grids`.

```sh
julia --project verify.jl
julia --project reconstruction_problem.jl
julia --project reconstruction_notation.jl
julia --project coefficient_stencils.jl
julia --project stencil_schematic.jl
julia --project worked_example.jl
julia --project advection.jl
julia --project stencil_flicker.jl
julia --project burgers_shock.jl
julia --project slotted_disk.jl
```

Outputs go to `output/`, which git ignores, unless `OUTPUT_DIRECTORY` names another directory. `STILL_ONLY=1` makes the
animation scripts write only their stills, and `SVT_PRESET` sets a faster encoder preset for drafts.

## Reconstruction

On a uniform grid of cells $I_i$ with averages $\bar{u}_i$, a stencil of $k$ cells reaching $r$ cells to the left of
cell $i$ gives the value at its right face

```math
u^-_{i+1/2} = \sum_{j=0}^{k-1} c_{r,j}\, \bar{u}_{i-r+j},
```

with the coefficients of Shu (1998), eq. (2.21). ENO starts from the cell itself and adds one neighbour at a time, on
whichever side gives the smaller undivided difference of the averages.

## Set-ups

- **Worked example:** the averages 1.0, 1.2, 1.3, 0.0, 0.1 and 0.3 of cells 1 to 6, with $k = 3$.
- **Advection:** $u_t + u_x = 0$, periodic, solved with the third-order Runge–Kutta method of Shu & Osher (1988) and
  $\Delta t = 0.4\,\Delta x$. Jiang & Shu's profile on $-1 \le x \le 1$ with 200 point values, to $t = 8$, at 25 frames
  per second with a frame every 3 steps; a Gaussian pulse $e^{-60x^2}$ on 100 cells, to $t = 2$, with and without
  uniform noise in $[0, 10^{-4}]$.
- **Fixed stencils:** each leans upwind, with $r = \lfloor k/2 \rfloor$ of its $k$ cells to the left of the target cell.
  For $k = 2$ and 3 these are the second- and third-order upwind schemes.
- **Burgers' equation:** $u_t + (u^2/2)_x = 0$ with $u_0(x) = 1/2 + \sin \pi x$, periodic on $-1 \le x \le 1$, on 80 cells
  to $t = 2$; the shock forms at $t = 1/\pi$. Both schemes are finite-volume schemes with the Godunov flux and
  $\Delta t = 0.4\,\Delta x / 1.5$. The exact solution follows the characteristics; in the frame moving with speed 1/2
  it is odd, so the shock stays where it forms.
- **Slotted disk:** a disk of radius 0.3 centered at $(0, 0.5)$ with a slot 0.1 wide up to $y = 0.7$ (Zalesak 1979,
  scaled), as cell averages on $100 \times 100$ cells of $-1 \le x, y \le 1$, turned by $(a, b) = \pi(-y, x)$, one turn
  every 2 time units, with CFL 0.4 for the sum of the two directions.
