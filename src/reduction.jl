export get_L0_science, get_L0_darks, get_L0_flats

# Detector info: http://irtfweb.ifa.hawaii.edu/~ishell/iSHELL_observing_manual_20210827.pdf
const ECHELLE_ORDERS_KGAS = 213:240
const DETECTOR_GAIN = 1.8
const DETECTOR_READ_NOISE = 10.0
const DETECTOR_DARK_CURRENT = 0.1

# Options
const DEFAULT_REDUCE_OPTIONS = Dict{String, Any}(

    # Start/end column for order tracing
    "trace_xrange" => [20, 2048-19],

    # Start/end row for order tracing in XD direction for center slice (x0)
    # This should be tailored to base_flat_field_file below regardless of the dataset being analyzed
    "trace_yrange" => [110, 1974],

    # The column for the starting slice for order tracing
    "trace_x0" => 1024,

    # Tracing polynomial deg
    "trace_poly_deg" => 4,

    # Which orders to trace
    "trace_orders" => ECHELLE_ORDERS_KGAS,

    # Whether or not to do flat or dark corrections
    "do_flat" => true, "do_dark" => false,

    # Minimum spacing in pixels between echelle orders
    "trace_min_spacing" => 5,

    # Trace position in each col is determined with a tophat convolved with a Gaussian.
    "trace_width_bounds" => [25, 40], # Tophat width bounds
    "trace_σ_bounds" => [1, 5], # Gaussian width bounds
    "trace_flux_smooth_width" => 300, # Smoothing width to remove scattered light
    "trace_fit_window" => 40, # Trace fit window

    # trace_yrange corresponds to this file (can be nothing if same as current dataset)
    "base_flat_field_file" => nothing,

    # Which orders are extracted
    "extract_orders" => ECHELLE_ORDERS_KGAS,

    # Bad pix threshold in extraction
    "extract_badpix_nσ" => 4,

    # Whether or not to remove background in extraction
    "extract_remove_background" => true,

    # Background smoothing width
    "extract_background_smooth_width" => 0,

    # What separates science vs. background flux. Multiples of second moment.
    "extract_yrange_thresh" => 5,

    # Polynomial degree fit to science trace positions
    "extract_trace_pos_deg" => 4,

    # Start/end column for extraction
    "extract_xrange" => [300, 2048-299],

    # Profile is modeled as a 1 or 2D cubic spline
    # Set profile_deg_x to 0 for 1D spline (and profile_knot_spacing_x is then irrelevant)
    "extract_profile_knot_spacing_x" => 2048,
    "extract_profile_knot_spacing_y" => 1,
    "extract_profile_deg_x" => 0,
    "extract_profile_deg_y" => 3
)

function reduce(;
        science::Vector{iSHELLL0},
        flats::Vector{Vector{iSHELLL0}},
        darks::Union{Vector{iSHELLL0}, Nothing}=nothing,
        output_path::String,
        options::Union{Dict, Nothing}=nothing
    )

    # Get options
    options = isnothing(options) ? copy(DEFAULT_REDUCE_OPTIONS) : merge(DEFAULT_REDUCE_OPTIONS, options)

    # Dirs
    if output_path[end] != PATHSEP
        output_path *= PATHSEP
    end
    sci_output_path = output_path * "spectra" * PATHSEP
    cal_output_path = output_path * "calib" * PATHSEP
    mkpath(sci_output_path)
    mkpath(cal_output_path)

    # Gen calib images
    median_dark = !isnothing(darks) ? gen_median_dark(darks, cal_output_path, options) : nothing
    median_flats = !isnothing(flats) ? gen_median_flats(flats, cal_output_path, options, median_dark) : nothing
    
    # Trace orders for appropriate images
    traces = trace(median_flats, cal_output_path, options, median_dark)

    # Remove blaze from flats
    if options["do_flat"]
        median_flats_norm = iSHELLL0[]
        continua = OrderedDict{String, <:Any}[]
        for i in eachindex(median_flats)
            r = remove_slitflat_continuum(median_flats[i], traces[i], cal_output_path, options)
            push!(median_flats_norm, r[1])
            push!(continua, r[2])
        end
        med_flats = median_flats_norm
    else
        med_flats = nothing
        continua = nothing
    end
    
    # Extract
    extract(science, traces, sci_output_path, options; dark=median_dark, flats=med_flats, continua=continua)

end

function get_L0_science(path::String; mode::String="KGAS")
    files = iSHELLL0.(sort(vcat(glob("*icm*.fits", path), glob("*spc*.fits", path))))
    targets = lowercase.(get_object.(files))
    modes = lowercase.(get_mode.(files))
    bad = findall(@. (targets == "dark") || (targets == "darks") || (targets == "flat") || (targets == "qth") || (targets == "arc") || (modes .!= lowercase(mode)))
    deleteat!(files, bad)
    if length(files) == 0
        error("No good science files found in $path")
    end
    jds = get_exposure_start_time.(files)
    ss = sortperm(jds)
    files .= files[ss]
    return files
end


function get_L0_darks(path::String; mode::String="KGAS")
    darks = iSHELLL0.(sort(glob("*dark*.fits", path)))
    targets = lowercase.(get_object.(darks))
    modes = lowercase.(get_mode.(darks))
    good = findall(@. ((targets == "dark") || (targets == "darks")) && (modes == mode))
    if length(good) == 0
        error("No good darks found in $path")
    end
    darks = darks[good]
    itimes = get_itime.(darks)
    ss = sortperm(itimes)
    darks .= darks[ss]
    return darks
end


function get_L0_flats(path::String; mode::String="KGAS")
    flats = iSHELLL0.(sort(glob("*flat*.fits", path)))
    targets = lowercase.(get_object.(flats))
    gaspos = lowercase.(get_gascell_pos.(flats))
    modes = lowercase.(get_mode.(flats))
    good = findall(@. ((targets == "flat") || (targets == "qth")) && (gaspos == "out") && (modes == lowercase(mode)))
    if length(good) == 0
        error("No good flats found in $path")
    end
    flats = flats[good]
    itimes = get_itime.(flats)
    ss = sortperm(itimes)
    flats .= flats[ss]
    flats = group_flats(flats)
    return flats
end


function EchelleReduce.gen_median_dark(darks::Vector{iSHELLL0}, output_path::String, options::Dict{String, <:Any})
    median_dark = iSHELLL0(get_median_dark_filename(darks, output_path))
    println("Generating median dark: $median_dark")
    dark_images = read.(darks)
    median_dark_image = EchelleReduce.gen_median_dark(dark_images)
    FITSIO.fitswrite(median_dark.filename, collect(transpose(med_dark_image)), header=deepcopy(read_header(darks[1])))
    return median_dark
end


function gen_median_flats(flats::Vector{Vector{iSHELLL0}}, output_path::String, options::Dict{String, <:Any}, dark::Union{iSHELLL0, Nothing}=nothing)
    median_flats = iSHELLL0[]
    dark_itime = !isnothing(dark) ? get_itime(dark) : nothing
    for _flats ∈ flats
        median_flat = iSHELLL0(get_median_flat_filename(_flats, output_path))
        println("Generating median flat: $median_flat")
        flat_images = read.(_flats)
        dark_image = !isnothing(dark) ? read(dark) : nothing
        flat_itime = get_itime(_flats[1])
        median_flat_image = gen_median_flat(flat_images; dark_image, q=0.75, flat_itime, dark_itime)
        FITSIO.fitswrite(median_flat.filename, collect(transpose(median_flat_image)), header=deepcopy(read_header(_flats[1])))
        push!(median_flats, median_flat)
    end
    return median_flats
end


function get_median_dark_filename(darks::Vector{iSHELLL0}, output_path::String)
    img_nums = get_image_num.(darks)
    img_start, img_end = minimum(img_nums), maximum(img_nums)
    itime = get_itime(darks[1])
    filename = "$(output_path)median_dark_$(get_utdate(darks[1]))_imgs$(img_start)-$(img_end)_$(itime)s.fits"
    return filename
end


function get_median_flat_filename(flats::Vector{iSHELLL0}, output_path::String)
    img_nums = get_image_num.(flats)
    img_start, img_end = minimum(img_nums), maximum(img_nums)
    filename = "$(output_path)median_flat_$(get_utdate(flats[1]))_imgs$(img_start)-$(img_end).fits"
    return filename
end


function trace(
        flats::Vector{iSHELLL0},
        output_path::String, options::AbstractDict{String, <:Any},
        darks::Union{iSHELLL0, Nothing}=nothing
    )

    # Store trace params for each median flat
    traces = OrderedDict{String, <:NamedTuple}[]

    # Loop over median flats
    for flat ∈ flats

        println("Tracing Orders: $flat")

        # Load image
        flat_image = read(flat, corr_readmath=false)

        # Labels
        labels = string.(options["trace_orders"])

        # Offset
        if "base_flat_field_file" in keys(options)
            offset = get_tracing_offset(flat_image, options["base_flat_field_file"])
            trace_yrange = options["trace_yrange"] .- offset
        else
            trace_yrange = options["trace_yrange"]
        end

        # Trace
        _traces = trace_boxcar(flat_image, labels;
                        xrange=options["trace_xrange"], yrange=trace_yrange, x0=options["trace_x0"],
                        width_bounds=options["trace_width_bounds"], σ_bounds=options["trace_σ_bounds"],
                        min_spacing=options["trace_min_spacing"], deg=options["trace_poly_deg"], fit_window=options["trace_fit_window"], flux_smooth_width=options["trace_flux_smooth_width"]
                   )

        # Plot image
        filename_out = "$(output_path)$(basename(flat)[1:end-5])_tracepos.png"
        plot_tracepos_image(flat_image, getproperty.(values(_traces), :yc), filename_out, qmax=0.9)

        # Save trace params to jld file
        filename_out = "$(output_path)$(splitext(basename(flat))[1])_traces.jld"
        jldsave(filename_out; traces=_traces)

        # Store results
        push!(traces, _traces)

    end

    # Return
    return traces

end


function create_output_dirs(output_path::String)
    mkpath("$(output_path)calib")
    mkpath("$(output_path)spectra")
end


function extract(
        data::iSHELLL0, traces::OrderedDict{String, <:NamedTuple}, options::Dict{String, <:Any};
        dark::Union{iSHELLL0, Nothing}=nothing,
        flat::Union{iSHELLL0, Nothing}=nothing,
        continua::Union{<:AbstractDict{String, <:Any}, Nothing}=nothing
    )

    # Load median dark and flat
    dark_image = !isnothing(dark) ? read(dark) : nothing
    flat_image = !isnothing(flat) ? read(flat, corr_readmath=false) : nothing

    # Load data
    data_image = read(data)

    # Calibrate the image
    itime = get_itime(data)
    dark_itime = !isnothing(dark) ? get_itime(dark) : nothing
    calibrate_image!(data_image; dark_image, itime, dark_itime, flat_image)

    # Fix negatives
    data_image[data_image .< 0] .= NaN

    # Convert to pe
    data_image .*= DETECTOR_GAIN
    
    # Extract traces
    reduced_data = OrderedDict{String, Any}()
    for trace in values(traces)
        if parse(Int, trace.label) in options["extract_orders"]
            extractor = get_extractor(data, traces[trace.label], options)
            reduced_data[trace.label] = extract_trace(data.filename, data_image, trace.label, extractor)
        end
    end

    # Correct continuum
    if options["do_flat"] && !isnothing(continua)
        for trace in values(traces)
            if parse(Int, trace.label) in options["extract_orders"] && !isnothing(reduced_data[trace.label])
                reduced_data[trace.label].spec ./= continua[trace.label]
                reduced_data[trace.label].specerr ./= continua[trace.label]
            end
        end
    end

    # Return
    return reduced_data

end


function extract(
        science::Vector{iSHELLL0},
        traces::Vector{<:OrderedDict{String, <:Any}},
        output_path::String, options::Dict{String, <:Any};
        dark::Union{iSHELLL0, Nothing}=nothing,
        flats::Union{Vector{iSHELLL0}, Nothing}=nothing,
        continua::Union{Vector{<:AbstractDict{String, <:Any}}, Nothing}=nothing,
    )

    # Extract in parallel
    pmap(science) do sci_data

        # Get cals
        dark = options["do_dark"] ? get_dark(sci_data, dark) : nothing
        flat = get_flat(sci_data, flats, options)
        k = findfirst((f) -> f == flat, flats)
        if options["do_flat"]
            _continua = continua[k]
        else
            flat = nothing
            _continua = nothing
        end

        # Get traces for this image
        _traces = traces[k]

        # Extract full image
        reduced_data = extract(sci_data, _traces, options; dark, flat, continua=_continua)

        # Save results
        save_extraction_results(sci_data, reduced_data, output_path)
 
    end
end


function save_extraction_results(data::iSHELLL0, reduced_data::OrderedDict{String, <:Any}, output_path::String)

    # Plot
    target = replace(get_object(data), " " => "_")
    filename_out = "$(output_path)$(splitext(basename(data))[1])_$(target)_reduced.png"
    plot_reduced_spectrum(reduced_data, filename_out)

    # Collect extracted spectra
    nx = 2048
    reduced_data_out = OrderedDict{String, Any}()
    for trace_label ∈ keys(reduced_data)
        if !isnothing(reduced_data[trace_label])
            reduced_data_out[trace_label] = [reduced_data[trace_label].spec reduced_data[trace_label].specerr]
        else
            reduced_data_out[trace_label] = fill(NaN, (nx, 2))
        end
    end

    # Save .fits file
    target = replace(get_object(data), " " => "_")
    fname = "$(output_path)$(splitext(basename(data))[1])_$(target)_reduced.fits"
    header = read_header(data)
    FITSIO.FITS(fname, "w") do f
        write(f, reduced_data_out, header=header)
    end

    # Save auxiliary data
    fname = "$(output_path)$(splitext(basename(data))[1])_$(target)_auxiliary.jld"
    save_auxiliary_data(fname, reduced_data)

end

function group_flats(flats::Vector{iSHELLL0})

    # Flat image nums
    flat_image_nums = get_image_num.(flats)

    # Find gaps in image nums
    s = findall(diff(flat_image_nums) .> 1)

    # One group of flats
    if length(s) == 0
        return [flats]
    end

    # Multiple groups
    flat_groups = Vector{iSHELLL0}[]
    push!(flat_groups, flats[1:s[1]])
    for i=1:length(s)-1
        push!(flat_groups, flats[s[i]+1:s[i+1]])
    end
    push!(flat_groups, flats[s[end]+1:end])
    return flat_groups
end


function get_tracing_offset(flat_image::Matrix{<:Real}, base_flat_field_file::String)
    image1 = FITSIO.FITS(base_flat_field_file) do f
        collect(transpose(Float64.(FITSIO.read(f[1]))))
    end
    ny = size(image1)[1]
    image1 .= quantile_filter(image1, window=(5, 5))
    image2 = quantile_filter(flat_image, window=(5, 5))
    lags = -100:100
    n_lags = length(lags)
    n_slices = 100
    ccfs = fill(NaN, (n_lags, n_slices))
    slices = collect(range(500, ny-500-1, length=n_slices)) .|> floor .|> Int
    for i=1:n_slices
        ii = slices[i]
        s1 = @views image1[:, ii] ./ nanquantile(image1[:, ii], 0.95)
        s2 = @views image2[:, ii] ./ nanquantile(image2[:, ii], 0.95)
        ccfs[:, i] .= cross_correlate_trace_cols(s1, s2, lags)
    end
    ccf = nanmedian(ccfs, dim=2)
    offset = lags[nanargmax(ccf)]
    return offset
end

function cross_correlate_trace_cols(s1, s2, lags)
    ccf = fill(NaN, length(lags))
    for i in eachindex(lags)
        s2s = circshift(s2, lags[i])
        good = findall(isfinite.(s1) .& isfinite.(s2s))
        ccf[i] = @views nansum(s1[good] .* s2s[good]) / nansum(s1[good] .^ 2)
    end
    return ccf
end


function get_dark(data::iSHELLL0, darks::Vector{iSHELLL0})
    itime = get_itime(data)
    for dark in darks
        if get_itime(dark) == itime
            return dark
        end
    end
end

function get_flat(data::iSHELLL0, flats::Vector{iSHELLL0}, options::AbstractDict{String, <:Any})

    # Sci target
    target = get_object(data)

    # Sci alt/az coord
    data_coord_altaz = radec2altaz(get_sky_coord(data), get_exposure_start_time(data), "irtf", spherical=true)
    
    # Median flat targets
    flat_targets = get_object.(flats)

    # Median flat alt/az coords
    flat_coords_altaz = [radec2altaz(get_sky_coord(f), get_exposure_start_time(f), "irtf", spherical=true) for f in flats]

    # Coord diffs
    ψs = [abs(acos(sin(f.θ) * sin(data_coord_altaz.θ) + cos(f.θ) * cos(data_coord_altaz.θ) * cos(data_coord_altaz.ϕ - f.ϕ))) for f ∈ flat_coords_altaz]

    # Matching targets
    good = findall(target .== flat_targets)

    # Use flat closest in angular sep and preferably with matching target name
    if length(good) > 0
        k = argmin(abs.(ψs[good]))
        return flats[good][k]
    else
        k = argmin(abs.(ψs))
        return flats[k]
    end

end


function get_read_noise(data::iSHELLL0)
    rn = DETECTOR_READ_NOISE
    header = read_header(data)
    if "NDRS" in keys(header) && header["NDRS"] > 1
        if header["NDRS"] == 32
            rn = 4.0
        end
    end
    return rn
end


function EchelleReduce.remove_slitflat_continuum(flat::iSHELLL0, traces::OrderedDict{String, <:NamedTuple}, output_path::String, options::Dict{String, <:Any})

    # Load header/image
    flat_image = read(flat, corr_readmath=false)

    # Remove continuum
    flat_image_norm, continua = remove_slitflat_continuum(flat_image, traces; xrange=options["trace_xrange"], deg=6)

    # Save new flat
    flat_norm = iSHELLL0(output_path * basename(flat.filename)[1:end-5] * "_norm.fits")
    FITSIO.fitswrite(
        flat_norm.filename,
        collect(transpose(flat_image_norm)),
        header=read_header(flat)
    )
    jldsave(output_path * basename(flat.filename)[1:end-5] * "_continua.jld"; continua)

    # Return
    return flat_norm, continua

end


function get_extractor(data::iSHELLL0, trace::NamedTuple, options::Dict{String, <:Any})
    read_noise = get_read_noise(data)
    return OptimalExtractor(;
        trace.yc, trace_pos_deg=options["extract_trace_pos_deg"],
        xrange=options["extract_xrange"], yrange=trace.yrange,
        yrange_extract=nothing, yrange_extract_thresh=options["extract_yrange_thresh"],
        profile_knot_spacing_x=options["extract_profile_knot_spacing_x"],
        profile_knot_spacing_y=options["extract_profile_knot_spacing_y"],
        profile_deg_x=options["extract_profile_deg_x"], profile_deg_y=options["extract_profile_deg_y"],
        remove_background=options["extract_remove_background"], background_smooth_width=options["extract_background_smooth_width"],
        max_iterations=10, read_noise,
        badpix_nσ=options["extract_badpix_nσ"]
    )
end

