# Dispersion relations of inertia–gravity waves on the five Arakawa grids:
#
# 1. In one dimension (waves along x), after Arakawa & Lamb (1977), Fig. 5 and eqs. 27–31: the frequency ν/f against kd
#    for the radius of deformation λ = 2d and λ = d/2.
# 2. In two dimensions, after their Figs. 8 and 9 and Randall (1994), Fig. 1: contours of ν/f over the wavenumber plane,
#    one figure for λ = 2d (fine) and one for λ = d/2 (coarse), as Randall's Fig. 2 compares λ/d = 2 and 0.1.
#
#     julia --project dispersion.jl

include("arakawa_grids.jl")
include("plotting.jl")

set_theme!(
        fontsize = 17,
        fonts = (; regular = "DejaVu Sans", bold = "DejaVu Sans Bold", italic = "DejaVu Sans Oblique"),
        Axis = (;
                xgridvisible = false,
                ygridvisible = false,
                titlefont = :regular,
                xlabelcolor = INK,
                ylabelcolor = INK,
                xticklabelcolor = MUTED,
                yticklabelcolor = MUTED
        )
)

function save_figure(name, fig)
    save(output_file("$name.svg"), fig)
    save(output_file("$name.png"), fig; px_per_unit = 2)
    return nothing
end

#####
##### One dimension
#####

function dispersion_axis(position, title; ylabel = "")
    return Axis(position;
            title,
            xlabel = "wave number  kd/π",
            ylabel,
            xticks = 0:0.25:1,
            limits = ((0, 1.17), (0, nothing))
    )
end

"""
    dispersion_1d!(ax, λ; nudges)

Plot the frequency ν/f of waves along x against kd/π, for the radius of deformation λ, on grids A–D and in the
continuous equations, with a label where each curve ends at kd = π. `nudges` moves labels up or down, apart from others
where two curves end close together.
"""
function dispersion_1d!(ax, λ; nudges)
    p = Problem(; λ)
    kd = range(0, π; length = 401)
    ν(grid) = [frequency(grid, k, 0.0, p) for k in kd] ./ p.f

    lines!(ax, kd ./ π, ν(Exact()); color = INK, linewidth = 3, label = "exact")

    for grid in (AGrid(), BGrid(), CGrid(), DGrid())
        lines!(ax, kd ./ π, ν(grid);
                color = GRID_COLORS[grid],
                linestyle = GRID_STYLES[grid],
                linewidth = 2.5,
                label = string(grid)
        )
    end

    for grid in (Exact(), AGrid(), BGrid(), CGrid(), DGrid())
        label_height = frequency(grid, π, 0.0, p) / p.f + get(nudges, grid, 0.0)
        text!(ax, 1.015, label_height; text = string(grid), align = (:left, :center), color = INK, fontsize = 16)
    end

    return nothing
end

function dispersion_1d()
    fig = Figure(size = (960, 430), figure_padding = (8, 18, 8, 8))

    resolved = dispersion_axis(fig[1, 1], "λ = 2d (deformation radius resolved)"; ylabel = "frequency  ν/f")
    dispersion_1d!(resolved, 2.0; nudges = Dict(BGrid() => 0.22, CGrid() => -0.22, DGrid() => 0.2))
    axislegend(resolved; position = :lt, framevisible = false, patchsize = (34, 10), rowgap = 0)

    unresolved = dispersion_axis(fig[1, 2], "λ = d/2 (deformation radius not resolved)")
    dispersion_1d!(unresolved, 0.5; nudges = Dict(AGrid() => 0.07, CGrid() => -0.07, DGrid() => 0.07))

    save_figure("dispersion_1d", fig)
    return nothing
end

#####
##### Two dimensions
#####

# The resolvable wavenumbers, in (kd/π, ld/π): the unit square on grids A–D. On grid E, the heights lie on a square
# lattice of spacing d turned through 45°, so they fill the triangle kd/π + ld/π ≤ √2.
side(grid) = 1.0
side(::EGrid) = sqrt(2)

resolvable(grid, k, l) = true
resolvable(::EGrid, k, l) = k + l <= sqrt(2) + 1e-9

# ν/f over the square 0 ≤ kd/π, ld/π ≤ side(grid), beyond grid E's triangle as well.
frequency_map(grid, p, ks) = [frequency(grid, k * π, l * π, p) / p.f for k in ks, l in ks]

#####
##### Contour labels
#####

# Makie places a label at the middle of each contour line and clips it at the edge of the panel. Here, each label goes
# as near the middle of its line as it can without crossing the edge of the panel (on grid E, of its triangle) or
# another label, and the line is broken around it. A label's box is (c, a, b, θ): its center c, half-width a and
# half-height b, turned through θ.

const Contours = Makie.Contours   # Contour.jl, which Makie's own contour plots use

label_text(level) = isinteger(level) ? string(Int(level)) : @sprintf("%.1f", level)

# The width of a label in ems, from the advance widths of DejaVu Sans.
text_width(text) = sum(c == '.' ? 0.318 : 0.636 for c in text)

# The distance of a point q from the edge of the resolvable region, negative outside it.
clearance(grid, q) = min(q[1], q[2], 1 - q[1], 1 - q[2])
clearance(::EGrid, q) = min(q[1], q[2], (sqrt(2) - q[1] - q[2]) / sqrt(2))

"""
    clip_to_region(grid, points)

Split a contour line into its runs inside the resolvable region, ending each where the line leaves it.
"""
clip_to_region(grid, points) = [points]

function clip_to_region(::EGrid, points)
    room(q) = sqrt(2) - q[1] - q[2]   # ≥ 0 inside the triangle
    crossing(a, b) = a + room(a) / (room(a) - room(b)) * (b - a)

    runs = Vector{Point2f}[]
    run = Point2f[]
    for (i, q) in enumerate(points)
        if room(q) >= 0
            i > 1 && room(points[i-1]) < 0 && push!(run, crossing(points[i-1], q))
            push!(run, q)
        elseif !isempty(run)
            push!(run, crossing(points[i-1], q))
            push!(runs, run)
            run = Point2f[]
        end
    end
    isempty(run) || push!(runs, run)

    return runs
end

"""
    split_at_corners(points; max_turn = deg2rad(30))

Split a contour line where it turns sharply. A smooth contour turns by about a degree from one point to the next. It
turns sharply only at a saddle, where two lines of the same level cross and the contour cuts the corner in two turns of
45°. Each arm is then labelled at its own middle.
"""
function split_at_corners(points; max_turn = deg2rad(30))
    arms = Vector{Point2f}[]
    start = 1
    for i in 2:length(points)-1
        u, v = points[i] - points[i-1], points[i+1] - points[i]
        turn = abs(atan(u[1] * v[2] - u[2] * v[1], u[1] * v[1] + u[2] * v[2]))
        if turn > max_turn
            push!(arms, points[start:i])
            start = i
        end
    end
    push!(arms, points[start:end])

    return arms
end

# The arc length along a line at each of its points.
cumulative_length(points) = [0.0; cumsum([hypot((points[i+1] - points[i])...) for i in 1:length(points)-1])]

# The point at arc length t along a line whose points have the arc lengths s.
function point_at(points, s, t)
    t <= 0 && return points[1]
    t >= s[end] && return points[end]
    i = searchsortedlast(s, t)
    return points[i] + Float32((t - s[i]) / (s[i+1] - s[i])) * (points[i+1] - points[i])
end

function corners(c, a, b, θ)
    u, v = Point2f(cos(θ), sin(θ)), Point2f(-sin(θ), cos(θ))
    return (c + a * u + b * v, c + a * u - b * v, c - a * u + b * v, c - a * u - b * v)
end

# Whether two boxes overlap: they do unless an edge direction of one of them separates them.
function overlap(box₁, box₂)
    (c₁, a₁, b₁, θ₁), (c₂, a₂, b₂, θ₂) = box₁, box₂
    half_extent(a, b, θ, n) = a * abs(cos(θ) * n[1] + sin(θ) * n[2]) + b * abs(cos(θ) * n[2] - sin(θ) * n[1])

    for θ in (θ₁, θ₁ + π / 2, θ₂, θ₂ + π / 2)
        n = (cos(θ), sin(θ))
        gap = abs((c₂[1] - c₁[1]) * n[1] + (c₂[2] - c₁[2]) * n[2])
        gap > half_extent(a₁, b₁, θ₁, n) + half_extent(a₂, b₂, θ₂, n) && return false
    end

    return true
end

function in_box(q, (c, a, b, θ))
    x, y = q[1] - c[1], q[2] - c[2]
    return abs(x * cos(θ) + y * sin(θ)) <= a && abs(y * cos(θ) - x * sin(θ)) <= b
end

"""
    contour_runs(grid, ks, ν, levels)

The contour lines of ν at `levels`, as (level, points) pairs, cut into runs inside the resolvable region and split at
sharp corners.
"""
function contour_runs(grid, ks, ν, levels)
    return [(level, arm) for level in levels
                         for curve in Contours.lines(Contours.contour(ks, ks, ν, level))
                         for points in clip_to_region(grid, Point2f.(Contours.coordinates(curve)...))
                         for arm in split_at_corners(points) if length(arm) >= 2]
end

"""
    label_place(grid, points, s, text, fontsize, boxes, other_points)

Where a label fits along a contour line whose points have the arc lengths `s`: as near the middle of the line as it can
be, clear of the edge of the resolvable region, of the other labels' `boxes` and of the other lines' points. Returns
the arc length of its center and its box, or `nothing` if there is no room.
"""
function label_place(grid, points, s, text, fontsize, boxes, other_points)
    L = s[end]
    gap = 0.25fontsize
    a, b = text_width(text) * fontsize / 2 + gap, 0.42fontsize + gap / 2   # half-sizes of the label, with its gap
    L > 2a + fontsize || return nothing

    for t in sort(range(a, L - a; length = 240); by = t -> abs(t - L / 2))
        c = point_at(points, s, t)
        p₁, p₂ = point_at(points, s, t - a), point_at(points, s, t + a)
        θ = atan(p₂[2] - p₁[2], p₂[1] - p₁[1])
        θ > π / 2 && (θ -= π)   # keep the text upright
        θ <= -π / 2 && (θ += π)
        box = (c, a, b, θ)

        minimum(clearance(grid, q) for q in corners(box...)) >= 0.15fontsize || continue
        any(overlap(box, other) for other in boxes) && continue
        any(in_box(q, box) for q in other_points) && continue

        return t, box
    end

    return nothing
end

"""
    labelled_contours!(ax, grid, ks, ν, levels, fontsize)

Contour lines of ν at `levels`, labelled with their values where there is room: clear of the edge, of the other labels
and of the other lines. `fontsize` is the size of the labels in data units.
"""
function labelled_contours!(ax, grid, ks, ν, levels, fontsize)
    runs = contour_runs(grid, ks, ν, levels)
    boxes = []

    for (n, (level, points)) in enumerate(runs)
        s = cumulative_length(points)
        text = label_text(level)
        other_points = (q for (m, (_, other)) in enumerate(runs) if m != n for q in other)
        place = label_place(grid, points, s, text, fontsize, boxes, other_points)

        if isnothing(place)
            lines!(ax, points; color = INK, linewidth = 1.1)
            continue
        end

        t, box = place
        c, a, _, θ = box
        push!(boxes, box)

        # The line, broken around the label.
        lines!(ax, vcat(points[s .< t - a], [point_at(points, s, t - a)]); color = INK, linewidth = 1.1)
        lines!(ax, vcat([point_at(points, s, t + a)], points[s .> t + a]); color = INK, linewidth = 1.1)
        text!(ax, c; text, rotation = θ, align = (:center, :center), fontsize, markerspace = :data, color = INK)
    end

    return nothing
end

"""
    dispersion_2d(name; λ, levels, colorrange, condition)

Contours of ν/f over the wavenumber plane, for the radius of deformation λ, on the five grids and in the continuous
equations, saved as `<name>.svg` and `<name>.png`.
"""
function dispersion_2d(name; λ, levels, colorrange, condition)
    p = Problem(; λ)
    fig = Figure(size = (980, 690))
    colormap = cgrad([colorant"#ffffff", colorant"#9ec5ef", colorant"#2a78d6"])

    for (n, grid) in enumerate((Exact(), GRIDS...))
        row, column = fldmod1(n, 3)
        k_max = side(grid)
        ks = range(0, k_max; length = 501)
        ν = frequency_map(grid, p, ks)

        ax = Axis(fig[row, column];
                title = string(grid),
                aspect = DataAspect(),
                limits = ((0, k_max), (0, k_max)),
                xticks = 0:0.5:1.5,
                yticks = 0:0.5:1.5,
                xlabel = row == 2 ? "kd/π" : "",
                ylabel = column == 1 ? "ld/π" : ""
        )

        # Outside grid E's triangle, the map is left blank.
        shown = [resolvable(grid, k, l) ? ν[i, j] : NaN for (i, k) in enumerate(ks), (j, l) in enumerate(ks)]
        heatmap!(ax, ks, ks, shown; colormap, colorrange, nan_color = :transparent, rasterize = 2)
        labelled_contours!(ax, grid, ks, ν, levels, 0.046k_max)   # the same size on every panel
        grid isa EGrid && lines!(ax, [0, sqrt(2)], [sqrt(2), 0]; color = MUTED, linewidth = 1)
    end

    Label(fig[0, :], "Frequency ν/f of inertia–gravity waves, for $condition"; fontsize = 17, color = INK)
    save_figure(name, fig)

    return nothing
end

dispersion_1d()

# Contours every 0.5 for λ = 2d, and every 0.1 for λ = d/2, where no frequency exceeds 2.5f. The level 1 is left out for
# λ = d/2: there the frequency on grid C is exactly f along both axes, and its contour would trace the edges.
dispersion_2d("dispersion_2d_fine";
        λ = 2.0,
        levels = 1:0.5:9,
        colorrange = (1, 6.5),
        condition = "λ = 2d (deformation radius resolved)"
)

dispersion_2d("dispersion_2d_coarse";
        λ = 0.5,
        levels = [0.1:0.1:0.9; 1.1:0.1:2.5],
        colorrange = (1, 2.5),
        condition = "λ = d/2 (deformation radius not resolved)"
)

println("saved dispersion_1d, dispersion_2d_fine and dispersion_2d_coarse (.svg and .png) to $output_directory")
