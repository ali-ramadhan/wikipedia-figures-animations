# Arakawa grids

Figures and animations for the Wikipedia article [Arakawa grids](https://en.wikipedia.org/wiki/Arakawa_grids).
Oceananigans.jl uses only the C grid, so the five grids are coded here directly, in plain Julia on the CPU, with the
same time stepping for all of them.

- `arakawa_grids.jl`: the five grids, as types that the schemes dispatch on, the Runge–Kutta time stepper, and the
  exact solutions used to check them.
- `verify.jl`: checks every scheme against its dispersion relation and against the exact solution of its
  semi-discrete equations (24 checks; about 15 s).
- `simulate.jl`: the two geostrophic-adjustment simulations, on all five grids (about 1 min).
- `animate.jl`: the two animations and a still from each, at 2360 × 2200 (about 6 min).
- `dispersion.jl`: the dispersion relations in one dimension, and in two for λ = 2d and λ = d/2. Contour labels are
  placed clear of the panel edges, of each other and of the other lines.
- `grid_schematic.jl`: the diagram of the five grids, as SVG (no packages needed).
- `plotting.jl`: the colors and line styles shared by the figures and the animations.
- `video_encoding.jl`: records the animations without loss, then encodes them for Wikimedia Commons as AV1 in WebM. It
  is adapted from `taylor-green-vortex-2d` in
  [oceananigans-for-wikipedia](https://github.com/ali-ramadhan/oceananigans-for-wikipedia).

```sh
julia --project verify.jl
julia --project simulate.jl && julia --project animate.jl
julia --project dispersion.jl
julia grid_schematic.jl
```

Outputs go to `output/`, which git ignores, unless `OUTPUT_DIRECTORY` names another directory. `STILL_ONLY=1` makes
`animate.jl` write the stills only.

## Problem

The linearized shallow-water equations on an *f*-plane,

```math
\frac{\partial u}{\partial t} - f v = -g \frac{\partial h}{\partial x}, \qquad
\frac{\partial v}{\partial t} + f u = -g \frac{\partial h}{\partial y}, \qquad
\frac{\partial h}{\partial t} + H \left( \frac{\partial u}{\partial x} + \frac{\partial v}{\partial y} \right) = 0,
```

the system Arakawa & Lamb (1977, pp. 180–183) use to compare the grids. Their inertia–gravity waves have the
frequency $\nu^2 = f^2 + gH(k^2 + l^2)$, and the radius of deformation is $\lambda = \sqrt{gH}/f$.

Each grid uses the simplest second-order scheme for its arrangement of $u$, $v$ and $h$: Arakawa & Lamb's eqs. (5)–(19),
with $\delta_x$ a difference over one grid interval and an overbar an average of neighbors. On grid D the Coriolis
terms are $-f\bar{v}^{xy}$ and $+f\bar{u}^{xy}$; eqs. (14) and (15) print the two averages the wrong way round (see
`checks/verify_dispersion.py` in the Wikipedia project). Grids A–D have $N \times N$ height points with spacing $d$.
Grid E puts heights and velocities on alternate points of an $M \times M$ lattice, so its differences span
$d^* = 2L/M$; with $N = 198$ and $M = 280$, $d^* = 1.41429\,d$, within 0.005% of Arakawa & Lamb's $\sqrt{2}\,d$, and all
five grids have the same number of heights. The domain is doubly periodic, with side $L = 198\,d$.

## Set-up

- **Initial condition:** a Gaussian mound $h = h_0\, e^{-r^2/2\sigma^2}$ with $\sigma = d$, at rest. The fluid
  adjusts to geostrophic balance: inertia–gravity waves carry part of the mound away, and a balanced vortex stays.
- **Two cases,** as in Arakawa & Lamb's conclusion that "the simulation of geostrophic adjustment is best with
  Scheme C, except for abnormal situations in which λ/d is less than or close to 1" (p. 190):
  - *fine:* $\lambda = 2d$, run to $ft = 15$ (2.4 inertial periods);
  - *coarse:* $\lambda = d/2$, run to $ft = 40$ (6.4 inertial periods).
- **Time stepping:** classical fourth-order Runge–Kutta with $f\Delta t = 0.005$, on every grid. Arakawa & Lamb keep
  time continuous in their analysis; with this step the time-stepping error is about $10^{-9}$ of the height, so the
  grids differ only in their space discretization.
- **Exact solution** (the "Exact" panels): each Fourier mode of the height evolves as
  $\hat{h}(t) = \hat{h}_0 \left[ f^2/\nu^2 + (1 - f^2/\nu^2) \cos \nu t \right]$, computed by FFT on a grid four
  times finer than $d$. The same formula with each grid's own $\nu$ and Coriolis factor gives the exact solution of
  that grid's semi-discrete equations, which `verify.jl` compares with the simulations.

## Accuracy

From `verify.jl` and the log of `simulate.jl`:

- Every grid's inertia–gravity eigenmodes rotate at the frequency of its dispersion relation (A&L eqs. 28–31 in one
  dimension; Cushman-Roisin & Beckers 2011, Table 9.1, in two), to within $10^{-8}$; halving $\Delta t$ cuts the
  error 16-fold, as fourth-order time stepping should.
- The simulations agree with the exact solutions of their semi-discrete equations to within $2 \times 10^{-9}$ of the
  initial height, conserve mass to round-off, and change the semi-discrete energy by less than $2 \times 10^{-10}$.
- On grid E, no height ever appears on the velocity points.
- The exact solution changes by less than $10^{-15}$ when its grid is refined from $d/4$ to $d/8$.

## What the animations show

- **Fine ($\lambda = 2d$):** the exact solution radiates a circular front. Grid C comes closest: its front lags a
  little and is slightly square, but little is left behind. On grids A and D much of the energy stays near the mound,
  where it oscillates; grids B and E leave grid-scale noise behind an octagonal (B) or diamond-shaped (E) front. The
  height at the center settles on C as in the exact solution, and keeps oscillating on A and D.
- **Coarse ($\lambda = d/2$):** grid C's waves hardly move, and its height at the center keeps oscillating with almost
  undiminished amplitude. Grids B and E come closest to the exact solution. On grid D, waves travel along the axes
  faster than any wave of the exact equations, $\sqrt{gH}$ (the dashed circle), with their energy moving against their
  crests; this follows from A&L's eq. (31), whose group velocity reaches about $0.7\,fd$ at $\lambda = d/2$, where
  $\sqrt{gH} = 0.5\,fd$. Arakawa met it in the first UCLA model, which used grid D: noise "rapidly propagated in space
  with an obviously nonphysical large group velocity" (Arakawa 2000, *A personal perspective on the early years of
  general circulation modeling at UCLA*, p. 32). Haltiner & Williams (1980, p. 228) note that on grids A and D "the
  short waves propagate energy in the wrong direction".

## Reference

Arakawa, A. & Lamb, V. R. (1977). Computational design of the basic dynamical processes of the UCLA general
circulation model. In *Methods in Computational Physics* 17, *General Circulation Models of the Atmosphere* (ed.
J. Chang), pp. 173–265. Academic Press. [doi:10.1016/B978-0-12-460817-7.50009-4](https://doi.org/10.1016/B978-0-12-460817-7.50009-4)
