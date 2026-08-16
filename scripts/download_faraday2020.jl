#!/usr/bin/env julia

using Downloads
using FITSIO
using Moose
using SHA

const SOURCE_URL =
    "https://wwwmpa.mpa-garching.mpg.de/ift/data/faraday2020/faraday2020v2.fits"
const SOURCE_SHA256 =
    "2540555659ddaf2fd9e44143c4ee5778dd04c922350b69f3fb3b759241aa97a0"
const NSIDE = 512

function file_sha256(path::AbstractString)
    return bytes2hex(open(sha256, path))
end

function main()
    output = joinpath(@__DIR__, "..", "data", "faraday2020v2.fits")
    mkpath(dirname(output))

    mktemp() do source, io
        close(io)
        println("Downloading Faraday 2020 data from the MPA…")
        Downloads.download(SOURCE_URL, source)

        actual_hash = file_sha256(source)
        actual_hash == SOURCE_SHA256 || error(
            "Downloaded file checksum mismatch: expected $(SOURCE_SHA256), got $(actual_hash).",
        )

        values = FITS(source) do fits
            Float64.(read(fits[2], "faraday_sky_mean"; case_sensitive=false))
        end
        length(values) == 12 * NSIDE^2 || error(
            "Unexpected map size $(length(values)); expected $(12 * NSIDE^2) pixels for NSIDE=$(NSIDE).",
        )

        Moose.write_healpix_map(
            output,
            values;
            nside=NSIDE,
            order=:ring,
            typechar="D",
            unit="rad/m^2",
            extname="MAP",
            coordsys="G",
            overwrite=true,
        )
    end

    println("Wrote MOOSE-compatible HEALPix map to $(abspath(output))")
end

main()
