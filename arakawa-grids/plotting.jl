# Colors and line styles shared by the figures (dispersion.jl) and the animations (animate.jl).

using CairoMakie

const INK = colorant"#0b0b0b"
const MUTED = colorant"#52514e"

# Categorical slots 1–5 of the dataviz reference palette (light mode), in the fixed order A–E.
const GRID_COLORS = Dict(
        AGrid() => colorant"#2a78d6",
        BGrid() => colorant"#eb6834",
        CGrid() => colorant"#1baf7a",
        DGrid() => colorant"#eda100",
        EGrid() => colorant"#e87ba4"
)

# The line styles of the one-dimensional dispersion figure, which leaves out grid E, and for E the one pattern left.
const GRID_STYLES = Dict(
        AGrid() => :dash,
        BGrid() => :dot,
        CGrid() => :solid,
        DGrid() => :dashdot,
        EGrid() => :dashdotdot
)
