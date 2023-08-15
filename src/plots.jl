export plotuvcov, plotvis, plotstationgains, plotbandpass, plotpointingerrors, plotelevationangle, plotparallacticangle, plotdterms, plottransmission,
plotmeandelays

"""
    plotuvcov(uvw::Matrix{Float64}, flagrow::Vector{Bool}, chanfreqvec::Vector{Float64}; saveprefix="test_")

Plot uv-coverage of observation.
"""
function plotuvcov(uvw::Matrix{Float64}, flagrow::Vector{Bool}, chanfreqvec::Vector{Float64}; saveprefix="test_")
    @info("Generating uv-coverage plots...")

    maskindices = findall(isequal(true), flagrow)
    muwave = deepcopy(uvw[1,:]) / (299792458.0/mean(chanfreqvec)) / 1e9 # convert to Gλ
    mvwave = deepcopy(uvw[2,:]) / (299792458.0/mean(chanfreqvec)) / 1e9

    muwave[maskindices] .= NaN
    mvwave[maskindices] .= NaN

    p = plot()
    plot!(p, muwave, mvwave, seriestype=:scatter, ls=:dot, ms=1, mc=:blue, msc=:blue, label="", xflip=true)
    plot!(p, xlabel="u (Gλ)", ylabel="v (Gλ)", title="uv-coverage")

    savefig(p, saveprefix*"uvcoverage.png")
    @info("Done 🙆")
end

"""
    plotvis(uvw::Matrix{Float64}, chanfreqvec::Array{Float64,1}, flag::Array{Bool,4}, data::Array{Complex{Float32},4},
    numchan::Int64, times::Vector{Float64}; saveprefix="data_")

Plot complex visibilities against time and projected baseline length.
"""
function plotvis(uvw::Matrix{Float64}, chanfreqvec::Array{Float64,1}, flag::Array{Bool,4}, data::Array{Complex{Float32},4},
    numchan::Int64, times::Vector{Float64}; plotphases::Bool=false, saveprefix="data_")
    @info("Generating visibility plots...")
    uvwave = sqrt.(uvw[1,:].^2 .+ uvw[2,:].^2) / (299792458.0/mean(chanfreqvec)) / 1e9 # in units of Gλ

    maskindices = findall(isequal(true), flag)
    maskeddata = deepcopy(data)
    maskeddata[maskindices] .= NaN # assign NaN values to flagged visibilities

    # plot visibility amplitudes against projected baseline separation
    p = plot()
    # plot first frequency channel
    plot!(p, uvwave, abs.(maskeddata[1,1,1,:]), seriestype=:scatter, ls=:dot, ms=1, mc=:red, msc=:red, label="RR")
    plot!(p, uvwave, abs.(maskeddata[1,2,1,:]), seriestype=:scatter, ls=:dot, ms=1, mc=:cyan, msc=:cyan, label="RL")
    plot!(p, uvwave, abs.(maskeddata[2,1,1,:]), seriestype=:scatter, ls=:dot, ms=1, mc=:purple, msc=:purple, label="LR")
    plot!(p, uvwave, abs.(maskeddata[2,2,1,:]), seriestype=:scatter, ls=:dot, ms=1, mc=:green, msc=:green, label="LL")
    
    # plot the rest of the frequency channels
    if numchan > 1
        plot!(p, uvwave, abs.(maskeddata[1,1,2:end,:]'), seriestype=:scatter, ls=:dot, ms=1, mc=:red, msc=:red, label="")
        plot!(p, uvwave, abs.(maskeddata[1,2,2:end,:]'), seriestype=:scatter, ls=:dot, ms=1, mc=:cyan, msc=:cyan, label="")
        plot!(p, uvwave, abs.(maskeddata[2,1,2:end,:]'), seriestype=:scatter, ls=:dot, ms=1, mc=:purple, msc=:purple, label="")
        plot!(p, uvwave, abs.(maskeddata[2,2,2:end,:]'), seriestype=:scatter, ls=:dot, ms=1, mc=:green, msc=:green, label="")
    end

    plot!(p, xlabel="Projected baseline separation (Gλ)", ylabel="Complex visibility amplitude (Jy)", legend=:outertop, legendcolumns=4)
    savefig(p, saveprefix*"visampvspbs.png")

    # plot visibility phases against projected baseline separation
    if plotphases
        p = plot()
        # plot first frequency channel
        plot!(p, uvwave, angle.(maskeddata[1,1,1,:]), seriestype=:scatter, ls=:dot, ms=1, mc=:red, msc=:red, label="RR")
        plot!(p, uvwave, angle.(maskeddata[1,2,1,:]), seriestype=:scatter, ls=:dot, ms=1, mc=:cyan, msc=:cyan, label="RL")
        plot!(p, uvwave, angle.(maskeddata[2,1,1,:]), seriestype=:scatter, ls=:dot, ms=1, mc=:purple, msc=:purple, label="LR")
        plot!(p, uvwave, angle.(maskeddata[2,2,1,:]), seriestype=:scatter, ls=:dot, ms=1, mc=:green, msc=:green, label="LL")
    
        # plot the rest of the frequency channels
        if numchan > 1
            plot!(p, uvwave, angle.(maskeddata[1,1,2:end,:]'), seriestype=:scatter, ls=:dot, ms=1, mc=:red, msc=:red, label="")
            plot!(p, uvwave, angle.(maskeddata[1,2,2:end,:]'), seriestype=:scatter, ls=:dot, ms=1, mc=:cyan, msc=:cyan, label="")
            plot!(p, uvwave, angle.(maskeddata[2,1,2:end,:]'), seriestype=:scatter, ls=:dot, ms=1, mc=:purple, msc=:purple, label="")
            plot!(p, uvwave, angle.(maskeddata[2,2,2:end,:]'), seriestype=:scatter, ls=:dot, ms=1, mc=:green, msc=:green, label="")
        end

        plot!(p, xlabel="Projected baseline separation (Gλ)", ylabel="Complex visibility phase (rad)", legend=:outertop, legendcolumns=4)
        savefig(p, saveprefix*"visphasevspbs.png")
    end

    # plot visibility amplitudes against time
    p = plot()
    x = times .- times[1]
    # plot first frequency channel
    plot!(p, x, abs.(maskeddata[1,1,1,:]), seriestype=:scatter, ls=:dot, ms=1, mc=:red, msc=:red, label="RR")
    plot!(p, x, abs.(maskeddata[1,2,1,:]), seriestype=:scatter, ls=:dot, ms=1, mc=:cyan, msc=:cyan, label="RL")
    plot!(p, x, abs.(maskeddata[2,1,1,:]), seriestype=:scatter, ls=:dot, ms=1, mc=:purple, msc=:purple, label="LR")
    plot!(p, x, abs.(maskeddata[2,2,1,:]), seriestype=:scatter, ls=:dot, ms=1, mc=:green, msc=:green, label="LL")
    
    # plot the rest of the frequency channels
    if numchan > 1
        plot!(p, x, abs.(maskeddata[1,1,2:end,:]'), seriestype=:scatter, ls=:dot, ms=1, mc=:red, msc=:red, label="")
        plot!(p, x, abs.(maskeddata[1,2,2:end,:]'), seriestype=:scatter, ls=:dot, ms=1, mc=:cyan, msc=:cyan, label="")
        plot!(p, x, abs.(maskeddata[2,1,2:end,:]'), seriestype=:scatter, ls=:dot, ms=1, mc=:purple, msc=:purple, label="")
        plot!(p, x, abs.(maskeddata[2,2,2:end,:]'), seriestype=:scatter, ls=:dot, ms=1, mc=:green, msc=:green, label="")
    end

    plot!(p, xlabel="Time (s)", ylabel="Complex visibility amplitude (Jy)", legend=:outertop, legendcolumns=4)
    savefig(p, saveprefix*"visampvstime.png")

    # plot visibility phases against time
    if plotphases
        p = plot()
        x = times .- times[1]
        # plot first frequency channel
        plot!(p, x, angle.(maskeddata[1,1,1,:]), seriestype=:scatter, ls=:dot, ms=1, mc=:red, msc=:red, label="RR")
        plot!(p, x, angle.(maskeddata[1,2,1,:]), seriestype=:scatter, ls=:dot, ms=1, mc=:cyan, msc=:cyan, label="RL")
        plot!(p, x, angle.(maskeddata[2,1,1,:]), seriestype=:scatter, ls=:dot, ms=1, mc=:purple, msc=:purple, label="LR")
        plot!(p, x, angle.(maskeddata[2,2,1,:]), seriestype=:scatter, ls=:dot, ms=1, mc=:green, msc=:green, label="LL")
        
        # plot the rest of the frequency channels
        if numchan > 1
            plot!(p, x, angle.(maskeddata[1,1,2:end,:]'), seriestype=:scatter, ls=:dot, ms=1, mc=:red, msc=:red, label="")
            plot!(p, x, angle.(maskeddata[1,2,2:end,:]'), seriestype=:scatter, ls=:dot, ms=1, mc=:cyan, msc=:cyan, label="")
            plot!(p, x, angle.(maskeddata[2,1,2:end,:]'), seriestype=:scatter, ls=:dot, ms=1, mc=:purple, msc=:purple, label="")
            plot!(p, x, angle.(maskeddata[2,2,2:end,:]'), seriestype=:scatter, ls=:dot, ms=1, mc=:green, msc=:green, label="")
        end
    
        plot!(p, xlabel="Time (s)", ylabel="Complex visibility phase (rad)", legend=:outertop, legendcolumns=4)
        savefig(p, saveprefix*"visphasevstime.png")
    end

    @info("Done 🙆")
end

"""
    plotstationgains(h5file::String, scanno::Vector{Int32}, times::Vector{Float64}, exposure::Float64, stationnames::Vector{String3})

Plot complex station gains against time.
"""
function plotstationgains(h5file::String, scanno::Vector{Int32}, times::Vector{Float64}, exposure::Float64, stationnames::Vector{String3})
    @info("Plotting station gains against time...")
    fid = h5open(h5file, "r")

    # get unique scan numbers
    uniqscans = unique(scanno)

    # get unique times
    uniqtimes = unique(times)

    p1amp = plot()
    p1phase = plot()
    p2amp = plot()
    p2phase = plot()
    for scan in uniqscans

        # determine indices of missing values
        actualtscanvec = unique(getindex(times, findall(scanno.==scan)))
    	actualtscanveclen = length(actualtscanvec)
	    idealtscanvec = collect(first(actualtscanvec):exposure:last(actualtscanvec))
	    idealtscanveclen = length(idealtscanvec)

        gterms = read(fid["stationgains"]["gjones_scan$(scan)"])
        #indexend = indexend + size(gterms)[3]

        # loop over time/row and apply gjones terms corresponding to each baseline
	    findnearest(A,x) = argmin(abs.(A .- x)) # define function to find nearest neighbour
        indvector = []
        for t in 1:actualtscanveclen
            idealtimeindex = findnearest(idealtscanvec, actualtscanvec[t])
            push!(indvector, idealtimeindex)
        end
            
        for ant in eachindex(stationnames)
            gpol1amp = abs.(gterms[1,1,indvector,ant]) # plot only the indices selected in the previous step
            gpol1phase = rad2deg.(angle.(gterms[1,1,indvector,ant]))
            gpol2amp = abs.(gterms[2,2,indvector,ant]) # plot only the indices selected in the previous step
            gpol2phase = rad2deg.(angle.(gterms[2,2,indvector,ant]))

            xvals = (actualtscanvec .- first(uniqtimes)) ./ 3600.0 # relative time in hours
            if scan == 1
                plot!(p1amp, xvals, gpol1amp, ls=:solid, lw=1, lc=ColorSchemes.mk_15[ant], label=stationnames[ant])
                plot!(p2amp, xvals, gpol2amp, ls=:solid, lw=1, lc=ColorSchemes.mk_15[ant], label=stationnames[ant])

                plot!(p1phase, xvals, gpol1phase, ls=:solid, lw=1, lc=ColorSchemes.mk_15[ant], label=stationnames[ant])
                plot!(p2phase, xvals, gpol2phase, ls=:solid, lw=1, lc=ColorSchemes.mk_15[ant], label=stationnames[ant])
            else
                plot!(p1amp, xvals, gpol1amp, ls=:solid, lw=1, lc=ColorSchemes.mk_15[ant], label="")
                plot!(p2amp, xvals, gpol2amp, ls=:solid, lw=1, lc=ColorSchemes.mk_15[ant], label="")

                plot!(p1phase, xvals, gpol1phase, ls=:solid, lw=1, lc=ColorSchemes.mk_15[ant], label="")
                plot!(p2phase, xvals, gpol2phase, ls=:solid, lw=1, lc=ColorSchemes.mk_15[ant], label="")
            end
        end
        #indexstart = indexend + 1
       
    end
    plot!(p1amp, title="Station gain amplitudes", ylabel="Pol1")
    plot!(p2amp, xlabel="Relative time (hr)", ylabel="Pol2")
    pamp = plot(p1amp, p2amp, layout=(2, 1))
    plot!(pamp, legend=:outertop, legendcolumns=6)

    plot!(p1phase, title="Station gain phases", ylabel="Pol1 (°)")
    plot!(p2phase, xlabel="Relative time (hr)", ylabel="Pol2 (°)")
    pphase = plot(p1phase, p2phase, layout=(2, 1))
    plot!(pphase, legend=:outertop, legendcolumns=6)
    
    savefig(pamp, "gainamplitudes_vs_time.png")
    savefig(pphase, "gainphases_vs_time.png")
    close(fid)
    @info("Done 🙆")
end

"""
    plotbandpass(h5file::String, stationnames::Vector{String3}, chanfreqvec::Vector{Float64})

Plot bandpass gains against time.
"""
function plotbandpass(h5file::String, stationnames::Vector{String3}, chanfreqvec::Vector{Float64})
    @info("Plotting bandpass gains against time...")
    fid = h5open(h5file, "r")

    b = read(fid["bandpass"]["bjonesmatrices"])

    p1 = plot()
    for ant in eachindex(stationnames)
        plot!(p1, chanfreqvec./1e9, abs.(b[1, 1, :, ant]), lw=1, lc=ColorSchemes.mk_15[ant], label="")
    end
    plot!(p1, title="Station bandpass gain amplitudes", ylabel="Pol1 gain amp") #, legend=:outertop, legendcolumns=6)

    p2 = plot()
    for ant in eachindex(stationnames)
        plot!(p2, chanfreqvec./1e9, abs.(b[2, 2, :, ant]), lw=1, lc=ColorSchemes.mk_15[ant], label=stationnames[ant])
    end
    plot!(p2, xlabel="Channel frequency (GHz)", ylabel="Pol2 gain amp") # legend=:outertop, legendcolumns=6)

    p = plot(p1, p2, layout=(2, 1))
    plot!(p, legend=:outertop, legendcolumns=6)
    savefig(p, "bpamplitudes_vs_frequency.png")
    close(fid)
    @info("Done 🙆")
end

"""
    plotpointingerrors(h5file::String, scanno::Vector{Int32}, stationnames::Vector{String3}; save::Bool=true)::Plots.Plot{Plots.GRBackend}

Plot pointing errors.
"""
function plotpointingerrors(h5file::String, scanno::Vector{Int32}, stationnames::Vector{String3}; save::Bool=true)::Plots.Plot{Plots.GRBackend}
    @info("Plotting pointing errors...")
    fid = h5open(h5file, "r")

    # get unique scan numbers
    uniqscans = unique(scanno)

    p_off = plot()
    p_amperr = plot()
    indexstart = 1
    indexend = 0
    for scan in uniqscans
        offsetarr = read(fid["pointingerrors"]["perscanoffsets_scan$(scan)"])
        amperrarr = read(fid["pointingerrors"]["perscanamperrors_scan$(scan)"])
        indexend = indexend + size(offsetarr)[1] # both arrays have same dimensions

        for ant in eachindex(stationnames)
            if scan == 1
                plot!(p_off, 1:indexend, offsetarr[:,ant], lw=1, lc=ColorSchemes.mk_15[ant], label=stationnames[ant])
                plot!(p_amperr, 1:indexend, amperrarr[:,ant], lw=1, lc=ColorSchemes.mk_15[ant], label=stationnames[ant])
            else
                plot!(p_off, indexstart:indexend, offsetarr[:,ant], lw=1, lc=ColorSchemes.mk_15[ant], label="")
                plot!(p_amperr, indexstart:indexend, amperrarr[:,ant], lw=1, lc=ColorSchemes.mk_15[ant], label="")
            end
        end
        indexstart = indexend + 1
       
    end
    plot!(p_off, xlabel="Time (mispointing number)", ylabel="Pointing offsets (arcsec)", legend=:outertop, legendcolumns=6)
    plot!(p_amperr, xlabel="Time (mispointing number)", ylabel="Primary beam response", legend=:outertop, legendcolumns=6)

    if save
        savefig(p_off, "pointingoffsets.png")
        savefig(p_amperr, "pointingamplitudeerrors.png")
    end
    
    close(fid)
    @info("Done 🙆")
    return p_amperr
end

"""
    plotelevationangle(h5file::String, scanno::Vector{Int32}, times::Vector{Float64}, stationnames::Vector{String3}; save::Bool=true)::Plots.Plot{Plots.GRBackend}

Plot evolution of station elevation angles during the course of the observation.
"""
function plotelevationangle(h5file::String, scanno::Vector{Int32}, times::Vector{Float64}, stationnames::Vector{String3}; save::Bool=true)::Plots.Plot{Plots.GRBackend}
    @info("Plotting elevation angles by station...")
    
    # get unique scan numbers
    uniqscans = unique(scanno)

    # get unique times
    uniqtimes = unique(times)
    x = (uniqtimes .- first(uniqtimes)) / 3600.0

    fid = h5open(h5file, "r")

    if "troposphere" in keys(fid)
        elevmat = read(fid["troposphere"]["elevation"])
    elseif "polarization" in keys(fid)
        elevmat = read(fid["polarization"]["elevation"])
    else
        close(fid)
        @error("$h5file does not contain elevation angle information. Not plotting elevation angles 🤷")
        throw(KeyError("elevation"))
    end

    if length(stationnames) != size(elevmat)[2]
        close(fid)
        throw(BoundsError("$h5file does not match the stations in station information. Not plotting elevation angles 🤷"))
    end

    if length(x) != size(elevmat)[1]
        close(fid)
        throw(DimensionMismatch("The number of timestamps in $h5file does not match that in the observation. Not plotting elevation angles 🤷"))
    end

    p = plot()
    for ant in eachindex(stationnames)
        plot!(p, x, rad2deg.(elevmat[:,ant]), seriestype=:scatter, ls=:dot, ms=1, mc=ColorSchemes.mk_15[ant], msc=ColorSchemes.mk_15[ant], label=stationnames[ant])
    end

    plot!(p, xlabel="Relative time (hr)", ylabel="Elevation angle (°)", legend=:outertop, legendcolumns=6)

    if save
    savefig(p, "elevationangle.png")
    end

    close(fid)
    @info("Done 🙆")
    return p
end

"""
    plotparallacticangle(h5file::String, scanno::Vector{Int32}, times::Vector{Float64}, stationnames::Vector{String3}; save::Bool=true)::Plots.Plot{Plots.GRBackend}

Plot evolution of station parallactic angles during the course of the observation.
"""
function plotparallacticangle(h5file::String, scanno::Vector{Int32}, times::Vector{Float64}, stationnames::Vector{String3}; save::Bool=true)::Plots.Plot{Plots.GRBackend}
    @info("Plotting parallactic angles by station...")

    # get unique scan numbers
    uniqscans = unique(scanno)

    # get unique times
    uniqtimes = unique(times)
    x = (uniqtimes .- first(uniqtimes)) / 3600.0

    fid = h5open(h5file, "r")

    if "polarization" in keys(fid)
        parangmat = read(fid["polarization"]["parallacticangle"])
    else
        close(fid)
        @error("$h5file does not contain parallactic angle information. Not plotting parallactic angles 🤷")
        throw(KeyError("parallacticangle"))
    end

    if length(stationnames) != size(parangmat)[2]
        close(fid)
        throw(BoundsError("$h5file does not match the stations in station information. Not plotting parallactic angles 🤷"))
    end

    if length(x) != size(parangmat)[1]
        close(fid)
        throw(DimensionMismatch("The number of timestamps in $h5file does not match that in the observation. Not plotting parallactic angles 🤷"))
    end

    p = plot()
    for ant in eachindex(stationnames)
        plot!(p, x, rad2deg.(parangmat[:,ant]), seriestype=:scatter, ls=:dot, ms=1, mc=ColorSchemes.mk_15[ant], msc=ColorSchemes.mk_15[ant], label=stationnames[ant])
    end

    plot!(p, xlabel="Relative time (hr)", ylabel="Parallactic angle (°)", legend=:outertop, legendcolumns=6)

    if save
        savefig(p, "parallacticangle.png")
    end

    close(fid)
    @info("Done 🙆")
    return p
end

"""
    plotdterms(h5file::String, stationnames::Vector{String3}, chanfreqvec::Vector{Float64}; save::Bool=true)::Plots.Plot{Plots.GRBackend}

Plot frequency-dependent complex instrumental polarization.
"""
function plotdterms(h5file::String, stationnames::Vector{String3}, chanfreqvec::Vector{Float64}; save::Bool=true)::Plots.Plot{Plots.GRBackend}
    @info("Plotting cross-hand instrumental leakage (D-terms) by station...")

    # read in the D-terms
    fid = h5open(h5file, "r")
    d = read(fid["polarization"]["djonesmatrices"])
    close(fid)

    chanfreqvec_ghz = chanfreqvec/1e9 # in GHz

    #=p = plot()
    for ant in eachindex(stationnames)
        plot!(p, real(d[1,2,:,ant]), imag(d[1,2,:,ant]), seriestype=:scatter, markershape=:circle, ms=3, mc=ColorSchemes.mk_15[ant], msc=ColorSchemes.mk_15[ant], label=stationnames[ant]*"-RL")
        plot!(p, real(d[2,1,:,ant]), imag(d[2,1,:,ant]), seriestype=:scatter, markershape=:diamond, ms=3, mc=ColorSchemes.mk_15[ant], msc=ColorSchemes.mk_15[ant], label=stationnames[ant]*"-LR")
    end

    plot!(p, xlabel="Real", ylabel="Imag", title="D-terms", legend=:outertop, legendcolumns=6)=#

    p = plot()
    for ant in eachindex(stationnames)
        plot!(p, chanfreqvec_ghz, abs.(d[1,2,:,ant]), seriestype=:scatter, markershape=:circle, ms=3, mc=ColorSchemes.mk_15[ant], msc=ColorSchemes.mk_15[ant], label=stationnames[ant]*"-RL")
        plot!(p, chanfreqvec_ghz, abs.(d[2,1,:,ant]), seriestype=:scatter, markershape=:diamond, ms=3, mc=ColorSchemes.mk_15[ant], msc=ColorSchemes.mk_15[ant], label=stationnames[ant]*"-LR")
    end
    plot!(p, xlabel="ν (GHz)", ylabel="D-term ampl.", legend=:outertop, legendcolumns=6)

    if save
        savefig(p, "dterms.png")
    end

    @info("Done 🙆")
    return p
end

"""
    plottransmission(h5file::String, stationnames::Vector{String3}, times::Vector{Float64}, chanfreqvec::Vector{Float64}; save::Bool=true)::Plots.Plot{Plots.GRBackend}

Plot tropospheric transmission variation with frequency.
"""
function plottransmission(h5file::String, stationnames::Vector{String3}, times::Vector{Float64}, chanfreqvec::Vector{Float64}; save::Bool=true)::Plots.Plot{Plots.GRBackend}
    @info("Plotting station-based transmission values")

    # get unique times
    uniqtimes = unique(times)
    reltimes = (uniqtimes .- first(uniqtimes))/3600.0 # in units of hours

    chanfreqvec_ghz = chanfreqvec/1e9 # in GHz

    # read in the computed transmission values
    fid = h5open(h5file, "r")
    tr = read(fid["troposphere"]["transmission"])
    close(fid)

    if length(chanfreqvec_ghz) == 1
        p = plot()
        for ant in eachindex(stationnames)
            plot!(p, reltimes, tr[1,:,ant], seriestype=:scatter, ls=:dot, ms=1, mc=ColorSchemes.mk_15[ant], msc=ColorSchemes.mk_15[ant], label=stationnames[ant])
        end
        plot!(p, xlabel="Relative time (hr)", ylabel="Transmission")
        if save
            savefig("transmission.png")
        end
    else
        plotarr = []
        for ant in eachindex(stationnames)
            xtr = tr[:,:,ant]
            p = contour(reltimes, chanfreqvec_ghz, xtr, title=stationnames[ant], xlabel="Rel. time (hr)", ylabel="ν (GHz)")
            push!(plotarr, p)
        end

        # plot all on the same figure
        nplots = length(plotarr)
        ncols = trunc(Int, ceil(sqrt(nplots)))
        nrows = trunc(Int, ceil(nplots/ncols))
        plotsize = (ncols*600, nrows*400)
        p = plot(plotarr...)
        plot!(p, layout=(nrows, ncols), size=plotsize)
        if save
            savefig(p, "transmission.png")
        end
    end

    @info("Done 🙆")
    return p
end

"""
    plotmeandelays(h5file::String, stationnames::Vector{String3}, times::Vector{Float64}, chanfreqvec::Vector{Float64}; save::Bool=true)::Plots.Plot{Plots.GRBackend}

Plot mean delays against time.
"""
function plotmeandelays(h5file::String, stationnames::Vector{String3}, times::Vector{Float64}, chanfreqvec::Vector{Float64}; save::Bool=true)::Plots.Plot{Plots.GRBackend}
    @info("Plotting station-based transmission values")

    # get unique times
    uniqtimes = unique(times)
    reltimes = (uniqtimes .- first(uniqtimes))/3600.0 # in units of hours

    chanfreqvec_ghz = chanfreqvec/1e9 # in GHz

    # read in the computed mean delay values
    fid = h5open(h5file, "r")
    md = read(fid["troposphere"]["meandelays"]) * 1e9 # delays / ns
    close(fid)

    p = plot()
    for ant in eachindex(stationnames)
        plot!(p, reltimes, dropdims(mean(md[:,:,ant], dims=1), dims=1), seriestype=:scatter, ls=:dot, ms=1, mc=ColorSchemes.mk_15[ant], msc=ColorSchemes.mk_15[ant], label=stationnames[ant])
    end
    plot!(p, xlabel="Relative time (hr)", ylabel="Mean delays (ns)")

    if save
        savefig("meandelays.png")
    end

    @info("Done 🙆")
    return p
end