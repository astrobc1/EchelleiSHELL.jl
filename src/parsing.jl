
# Read headers
Echelle.read_header(data::iSHELLL0) = read_header(data.filename, 1)
Echelle.read_header(data::iSHELLL1) = read_header(data.filename, 2)
Echelle.read_key(data::iSHELLL0, key::Union{Int, String}) = read_key(data.filename, key, hdu=1)
Echelle.read_key(data::iSHELLL1, key::Union{Int, String}) = read_key(data.filename, key, hdu=2)

# Header parsing
Echelle.get_itime(data::AnyiSHELL) = read_key(data, "ITIME")
Echelle.get_object(data::AnyiSHELL) = read_key(data, "OBJECT")
Echelle.get_utdate(data::AnyiSHELL) = join(split(read_key(data, "DATE_OBS"), '-'))
Echelle.get_sky_coord(data::AnyiSHELL) = ICRSCoords(hms2rad(read_key(data, "TCS_RA")), dms2rad(read_key(data, "TCS_DEC")))
Echelle.get_exposure_start_time(data::AnyiSHELL) = read_key(data, "TCS_UTC") + 2400000.5
Echelle.get_image_num(data::AnyiSHELL) = parse(Int, split(split(read_key(data, "IRAFNAME"), PATHSEP)[end], '.')[5])
get_gascell_pos(data::AnyiSHELL) = read_key(data, "GASCELL")
get_mode(data::AnyiSHELL) = read_key(data, "XDTILT")

# Read in image
function Base.read(data::iSHELLL0; hdu::Int=1, trans=true, mask_edges=true, corr_readmath::Bool=true)
    image = FITS(data.filename) do f
        Float64.(read(f[hdu]))
    end
    if trans
        image .= collect(transpose(image))
    end
    if mask_edges
        image[1:4, :] .= NaN
        image[2018:end, :] .= NaN
        image[:, 1:4] .= NaN
        image[:, end-3:end] .= NaN
    end
    if corr_readmath
        image .= correct_readmath(data, image)
    end
    return image
end

function correct_readmath(data::iSHELLL0, image::Matrix)
    image = copy(image)
    h = read_header(data)
    if "BZERO" in keys(h)
        image .= image .- h["BZERO"]
    end
    if "BSCALE" in keys(h)
        image .= image ./ h["BSCALE"]
    end
    return image
end

# Read in reduced spectrum and error, wavelength separate
function Base.read(
        data::iSHELLL1;
        order::Int, flip=true, hdu::Int=2, pix_range::Union{Vector{Int}, Nothing}=nothing, norm::Union{Real, Nothing}=nothing
    )
    spec, specerr = FITS(data.filename) do file
        d = read(file[hdu], "$order")
        return d[:, 1], d[:, 2]
    end
    if flip
        reverse!(spec)
        reverse!(specerr)
    end
    if !isnothing(pix_range)
        spec[.~(pix_range[1] .< eachindex(spec) .< pix_range[2])] .= NaN
        specerr[.~(pix_range[1] .< eachindex(specerr) .< pix_range[2])] .= NaN
    end
    if !isnothing(norm)
        v = nanquantile(Echelle.quantile_filter(spec, window=5), norm)
        spec ./= v
        specerr ./= v
    end
    return spec, specerr
end