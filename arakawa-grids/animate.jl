# Animations of the geostrophic adjustment simulated by simulate.jl, for Wikimedia Commons. Each shows the height on the
# five grids beside the exact solution of the continuous equations and, underneath, the height at the center of the
# mound as it settles towards its geostrophically balanced value (after Arakawa & Lamb 1977, Figs. 6 and 7, which show
# the one-dimensional problem). It writes geostrophic_adjustment_<case>.webm and a still at a chosen time,
# geostrophic_adjustment_<case>_ft<t>.png, both at 2 pixels per unit: 2360 × 2200. STILL_ONLY=1 writes only the stills.
#
#     julia --project animate.jl

include("arakawa_grids.jl")
include("plotting.jl")
include("video_encoding.jl")

using JLD2

const PX_PER_UNIT = 2

# Heights beyond the color range take the colors at its ends.
const BALANCE = Makie.to_colormap(:balance)

set_theme!(
        fontsize = 17,
        fonts = (; regular = "DejaVu Sans", bold = "DejaVu Sans Bold"),
        Axis = (;
                xgridvisible = false,
                ygridvisible = false,
                titlefont = :bold,
                titlesize = 19,
                xticklabelcolor = MUTED,
                yticklabelcolor = MUTED
        )
)

fine = (
        name = "fine",
        subtitle = "Radius of deformation λ = 2d, twice the grid spacing",
        colorrange = (-0.06, 0.06),
        still_time = 10.0,
        pixels_per_d = 12,
        ticks = -20:20:20,
        time_ticks = 0:0.5:2
)

coarse = (
        name = "coarse",
        subtitle = "Radius of deformation λ = d/2, half the grid spacing",
        colorrange = (-0.04, 0.04),
        still_time = 30.0,
        pixels_per_d = 16,
        ticks = -10:10:10,
        time_ticks = 0:1:6
)

#####
##### Maps of the height
#####

# The offsets (di, dj), in lattice points, of the height point nearest to (x, y), on a lattice of spacing s.
nearest_height_point(grid, x, y, s) = round(Int, x / s), round(Int, y / s)

# On grid E, the heights are on the points whose offsets from the central height have an even sum.
function nearest_height_point(::EGrid, x, y, s)
    P, Q = 2round((x + y) / 2s), 2round((x - y) / 2s)
    return round(Int, (P + Q) / 2), round(Int, (P - Q) / 2)
end

"""
    nearest_heights(grid, offsets, W, pixels_per_d)

The pixel coordinates of an image of [-W, W]², and a lookup table from its pixels to the nearest height point in the
saved window of `grid`, whose points have the coordinates `offsets` in each direction. Each height fills a square on
grids A–D and a diamond on grid E.
"""
function nearest_heights(grid, offsets, W, pixels_per_d)
    xs = range(-W, W; length = round(Int, 2W * pixels_per_d) + 1)
    s = offsets[2] - offsets[1]
    r = (length(offsets) - 1) ÷ 2   # the window is offsets[r + 1] = 0 ± r points

    index = Matrix{CartesianIndex{2}}(undef, length(xs), length(xs))
    for (b, y) in enumerate(xs), (a, x) in enumerate(xs)
        di, dj = nearest_height_point(grid, x, y, s)
        index[a, b] = CartesianIndex(clamp(r + 1 + di, 1, 2r + 1), clamp(r + 1 + dj, 1, 2r + 1))
    end

    return xs, index
end

height_heatmap!(ax, x, y, h, colorrange) =
    heatmap!(ax, x, y, h; colormap = :balance, colorrange, highclip = BALANCE[end], lowclip = BALANCE[1])

"""
    height_map!(ax, grid, p, data, frame, case)

Plot the height at output `frame[]` of the saved run `data` on `grid`, or of the exact solution if `grid = Exact()`.
"""
function height_map!(ax, grid::ArakawaGrid, p, data, frame, case)
    xs, index = nearest_heights(grid, data["offsets_$grid"], data["window"], case.pixels_per_d)
    h = data["h_$grid"]
    return height_heatmap!(ax, xs, xs, @lift(Float64.(view(h, :, :, $frame)[index])), case.colorrange)
end

# The exact solution, on a grid eight times finer than d to stay smooth at 2 pixels per unit.
function height_map!(ax, ::Exact, p, data, frame, case)
    times = data["t"]
    W = data["window"]

    window = @lift begin
        x, y, h = exact_h(p, times[$frame]; refine = 8)
        inside = findall(xᵢ -> abs(xᵢ) <= W, x)
        (x[inside], y[inside], h[inside, inside])
    end

    return height_heatmap!(ax, @lift($window[1]), @lift($window[2]), @lift($window[3]), case.colorrange)
end

function map_axis(position, title, W, ticks; row, column)
    return Axis(position;
            title,
            aspect = DataAspect(),
            limits = ((-W, W), (-W, W)),
            xlabel = row == 2 ? "x/d" : "",
            ylabel = column == 1 ? "y/d" : "",
            xticks = ticks,
            yticks = ticks,
            xticklabelsvisible = row == 2,
            yticklabelsvisible = column == 1
    )
end

#####
##### The height at the center of the mound
#####

# The exact height at the center of the mound, x = y = 0.
function exact_center_height(p, t)
    x, y, h = exact_h(p, t; refine = 2)
    return h[findfirst(==(0.0), x), findfirst(==(0.0), y)]
end

"""
    center_height_panel!(layout, p, data, periods, frame, time_ticks)

Plot the height at the center of the mound on each grid and in the exact solution, up to output `frame[]`, against
time in inertial periods.
"""
function center_height_panel!(layout, p, data, periods, frame, time_ticks)
    ax = Axis(layout[1, 1];
            xlabel = "time (inertial periods)",
            ylabel = "h/h₀ at the center",
            limits = ((0, periods[end]), nothing),
            xticks = time_ticks
    )

    exact = [exact_center_height(p, t) for t in data["t"]]
    up_to_frame(h) = @lift Point2f.(periods[1:$frame], h[1:$frame])

    for grid in GRIDS
        lines!(ax, up_to_frame(data["center_$grid"]);
                color = GRID_COLORS[grid],
                linestyle = GRID_STYLES[grid],
                linewidth = 2.5,
                label = string(grid)
        )
    end
    lines!(ax, up_to_frame(exact); color = INK, linewidth = 3, label = "exact")

    heights = vcat(exact, (data["center_$grid"] for grid in GRIDS)...)
    ylims!(ax, minimum(heights) - 0.05, maximum(heights) + 0.05)
    Legend(layout[1, 2], ax; framevisible = false, patchsize = (36, 10), rowgap = 2, labelcolor = INK)

    return nothing
end

#####
##### The animations
#####

"""
    save_animation(figure, frame, times, name; still_time)

Save `<name>_ft<still_time>.png`, a still at `still_time`, and then, unless STILL_ONLY=1, step the observable `frame`
through `times` to write the lossless master `<name>.mkv` and the video `<name>.webm`.
"""
function save_animation(figure, frame, times, name; still_time)
    frame[] = argmin(abs.(times .- still_time))
    save(output_file(@sprintf("%s_ft%g.png", name, still_time)), figure; px_per_unit = PX_PER_UNIT)
    get(ENV, "STILL_ONLY", "") == "1" && return nothing

    master = output_file("$name.mkv")
    seconds = @elapsed record_lossless(figure, master, eachindex(times); framerate = 24, px_per_unit = PX_PER_UNIT) do n
        frame[] = n
    end
    encode_webm(master, output_file("$name.webm"))
    @printf("%s: %d frames in %.0f s\n", name, length(times), seconds)

    return nothing
end

function animate(case)
    data = load(output_file("adjustment_$(case.name).jld2"))
    p = Problem(; λ = data["λ"], σ = data["σ"], N = data["N"], M = data["M"])
    times = data["t"]
    periods = p.f .* times ./ 2π   # time in inertial periods, of 2π/f each
    frame = Observable(1)

    fig = Figure(size = (1180, 1100), figure_padding = (12, 16, 10, 10))
    Label(fig[1, 1:3], "Geostrophic adjustment on the five Arakawa grids"; fontsize = 21, font = :bold, color = INK)
    subtitle = @lift @sprintf("%s.   Time: %.2f inertial periods", case.subtitle, periods[$frame])
    Label(fig[2, 1:3], subtitle; fontsize = 17, color = MUTED)

    # The maps, each with a dashed circle at the front of the fastest waves of the continuous equations, at speed √(gH).
    maps = fig[3, 1:3] = GridLayout()
    c = sqrt(gravity(p) * p.H)
    front = @lift [Point2f(c * times[$frame] * cos(θ), c * times[$frame] * sin(θ)) for θ in range(0, 2π; length = 241)]

    for (n, grid) in enumerate((Exact(), GRIDS...))
        row, column = fldmod1(n, 3)
        ax = map_axis(maps[row, column], uppercasefirst(string(grid)), data["window"], case.ticks; row, column)
        height_map!(ax, grid, p, data, frame, case)
        lines!(ax, front; color = (INK, 0.55), linestyle = :dash, linewidth = 1.2)
    end

    # The triangles at the ends of the color bar mark heights beyond its range, which the maps saturate.
    Colorbar(maps[1:2, 4];
            colormap = :balance,
            limits = case.colorrange,
            label = "height h/h₀",
            highclip = BALANCE[end],
            lowclip = BALANCE[1],
            height = Relative(0.85),
            labelcolor = INK,
            ticklabelcolor = MUTED
    )

    bottom = fig[4, 1:3] = GridLayout()
    center_height_panel!(bottom, p, data, periods, frame, case.time_ticks)
    rowsize!(fig.layout, 4, Relative(0.23))

    save_animation(fig, frame, times, "geostrophic_adjustment_$(case.name)"; still_time = case.still_time)
    return nothing
end

foreach(animate, (fine, coarse))
