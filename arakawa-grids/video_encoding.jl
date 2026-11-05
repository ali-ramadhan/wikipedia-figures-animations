# Record Makie animations without loss, and encode them for Wikimedia Commons.
#
# CairoMakie's `record` encodes WebM in a single fast pass (VP9 at CRF 20) without color tags. Instead,
# `record_lossless` writes every frame without loss to an FFV1 master, and `encode_webm` encodes the master as AV1 with
# SVT-AV1's slowest preset: at equal or better quality, the files are about 30% smaller, and their BT.709 color tags
# make every player show the same colors. Both functions call the system's ffmpeg, which must include the SVT-AV1
# encoder (libsvtav1).

using CairoMakie: colorbuffer

next_even(n) = n + isodd(n)

"""
    record_lossless(update!, figure, path, iterator; framerate = 30, px_per_unit = 1)

Like `CairoMakie.record`, but write the frames without loss, as FFV1 in a Matroska file at `path`. The figure is first
enlarged to even dimensions, which `encode_webm` needs, and `px_per_unit = 2` renders every frame at twice its size.
"""
function record_lossless(update!, figure, path, iterator; framerate = 30, px_per_unit = 1)
    resize!(figure, next_even.(size(figure.scene))...)
    width, height = round.(Int, px_per_unit .* size(figure.scene))

    ffmpeg = `ffmpeg -y -loglevel error
                     -f rawvideo -pixel_format bgr0 -video_size $(width)x$(height) -framerate $framerate -i pipe:0
                     -c:v ffv1 -level 3 -slices 24 $path`

    open(ffmpeg, "w") do io
        for i in iterator
            update!(i)
            # colorbuffer is a height × width matrix of ARGB32, stored as bytes B, G, R, A;
            # its transpose holds the image row by row, which is ffmpeg's bgr0 layout.
            write(io, permutedims(colorbuffer(figure; px_per_unit)))
        end
    end

    return path
end

"""
    encode_webm(master, path; crf = 20, preset = 0)

Encode the lossless `master` as AV1 in WebM, the format Wikimedia Commons prefers, with 10-bit 4:2:0 color converted
and tagged as BT.709. `preset = 0` is SVT-AV1's slowest and most efficient preset, and a lower `crf` gives higher
quality and larger files.
"""
function encode_webm(master, path; crf = 20, preset = 0)
    color_filters = join([
            "scale=out_color_matrix=bt709:out_range=tv:flags=lanczos+accurate_rnd",
            "format=yuv420p10le",
            "setparams=color_primaries=bt709:color_trc=bt709:colorspace=bt709:range=tv"
    ], ",")

    ffmpeg = `ffmpeg -y -loglevel error -i $master -vf $color_filters
                     -c:v libsvtav1 -preset $preset -crf $crf -g 240 -svtav1-params tune=0 $path`

    # SVT-AV1 prints its own log unless told not to.
    run(addenv(ffmpeg, "SVT_LOG" => "1"))

    return path
end
