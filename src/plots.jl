export plotvis, plotstationgains, plotbandpass, plotpointingerrors, plotelevationangle, plotparallacticangle, plotdterms,
plottransmission, plotmeandelays

"""
    plotvis(uvw::Matrix{Float64}, chanfreqvec::Array{Float64,1}, flag::Array{Bool,4}, data::Array{Complex{Float32},4},
    numchan::Int64, times::Vector{Float64}; saveprefix="data_")

Plot complex visibilities against time and projected baseline length.
"""
function plotvis(uvw::Matrix{Float64}, chanfreqvec::Array{Float64,1}, flag::Array{Bool,4}, data::Array{Complex{Float32},4},
    numchan::Int64, times::Vector{Float64}; plotphases::Bool=false, saveprefix="data")
    @info("Generating visibility plots...")
    uvwave = sqrt.(uvw[1,:].^2 .+ uvw[2,:].^2) / (299792458.0/mean(chanfreqvec)) / 1e9 # in units of Gλ

    maskindices = findall(isequal(true), flag)
    maskeddata = deepcopy(data)
    maskeddata[maskindices] .= NaN # assign NaN values to flagged visibilities

    xtimes = times .- times[1] # relative time in seconds

    #set pol colours
    colours = [[:red, :cyan], [:purple, :green]]
    labels = [["RR", "RL"], ["LR", "LL"]]

    # plot visibility amplitudes against projected baseline separation
    f = Figure(size=(800, 400))
    ax = Axis(f[1, 1], xlabel="Projected baseline separation (Gλ)", ylabel="Complex visibility amplitude (Jy)")
    for ii in 1:2
        for jj in 1:2
            scatter!(ax, uvwave, abs.(maskeddata[ii,jj,1,:]), color=colours[ii][jj], label=labels[ii][jj], markersize=2)
            if numchan > 1
                for channo in 2:numchan
                    scatter!(ax, uvwave, abs.(maskeddata[ii,jj,channo,:]), color=colours[ii][jj], markersize=2)
                end
            end
        end
    end

    axislegend(ax, position=:rt, merge=true, unique=true, tellheight=true, tellwidth=true)
    save(saveprefix*"_visibilityamplitude_vs_baseline.png", f)

    # plot visibility phases against projected baseline separation
    if plotphases
        f = Figure(size=(800, 400))
        ax = Axis(f[1, 1], xlabel="Projected baseline separation (Gλ)", ylabel="Complex visibility phase (deg.)")
        for ii in 1:2
            for jj in 1:2
                scatter!(ax, uvwave, rad2deg.(angle.(maskeddata[ii,jj,1,:])), color=colours[ii][jj], label=labels[ii][jj], markersize=2)
                if numchan > 1
                    for channo in 2:numchan
                        scatter!(ax, uvwave, rad2deg.(angle.(maskeddata[ii,jj,channo,:])), color=colours[ii][jj], markersize=2)
                    end
                end
            end
        end

        axislegend(ax, position=:rt, merge=true, unique=true, tellheight=true, tellwidth=true)
        save(saveprefix*"_visibilityphase_vs_baseline.png", f)
    end

    # plot visibility amplitudes against time
    f = Figure(size=(800, 400))
    ax = Axis(f[1, 1], xlabel="Projected baseline separation (Gλ)", ylabel="Complex visibility amplitude (Jy)")
    for ii in 1:2
        for jj in 1:2
            scatter!(ax, xtimes, abs.(maskeddata[ii,jj,1,:]), color=colours[ii][jj], label=labels[ii][jj], markersize=2)
            if numchan > 1
                for channo in 2:numchan
                    scatter!(ax, xtimes, abs.(maskeddata[ii,jj,channo,:]), color=colours[ii][jj], label="", markersize=2)
                end
            end
        end
    end

    axislegend(ax, position=:rt, merge=true, unique=true, tellheight=true, tellwidth=true)
    save(saveprefix*"_visibilityamplitude_vs_time.png", f)

    # plot visibility phases against time
    if plotphases
        f = Figure(size=(800, 400))
        ax = Axis(f[1, 1], xlabel="Projected baseline separation (Gλ)", ylabel="Complex visibility phase (deg.)")
        for ii in 1:2
            for jj in 1:2
                scatter!(ax, xtimes, rad2deg.(angle.(maskeddata[ii,jj,1,:])), color=colours[ii][jj], label=labels[ii][jj], markersize=2)
                if numchan > 1
                    for channo in 2:numchan
                        scatter!(ax, xtimes, rad2deg.(angle.(maskeddata[ii,jj,channo,:])), color=colours[ii][jj], label="", markersize=2)
                    end
                end
            end
        end

        axislegend(ax, position=:rt, merge=true, unique=true, tellheight=true, tellwidth=true)
        save(saveprefix*"_visibilityphase_vs_time.png", f)
    end

    @info("Plotted visibilities 🙆")
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

    f = Figure(size=(500, 300))
    axphase1 = Axis(f[1, 1], ylabel="Phases (°)", title="Station gains (Pol1)")
    axamp1 = Axis(f[2, 1], xlabel="Relative time (hr)", ylabel="Amplitudes")
    axphase2 = Axis(f[1, 2], title="Station gains (Pol2)")
    axamp2 = Axis(f[2, 2], xlabel="Relative time (hr)")

    # modify axis attributes
    hidexdecorations!(axphase1, grid=false, ticks=false)
    hidexdecorations!(axphase2, grid=false, ticks=false)
    hideydecorations!(axamp2, grid=false, ticks=false)
    hideydecorations!(axphase2, grid=false, ticks=false)

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
            lines!(axamp1, xvals, gpol1amp, ls=:solid, lw=1, color=ColorSchemes.mk_15[ant], label=stationnames[ant])
            lines!(axphase1, xvals, gpol1phase, ls=:solid, lw=1, color=ColorSchemes.mk_15[ant], label=stationnames[ant])

            lines!(axamp2, xvals, gpol2amp, ls=:solid, lw=1, color=ColorSchemes.mk_15[ant], label=stationnames[ant])
            lines!(axphase2, xvals, gpol2phase, ls=:solid, lw=1, color=ColorSchemes.mk_15[ant], label=stationnames[ant])
        end
        #indexstart = indexend + 1
       
    end
    close(fid) # close HDF5 file

    linkxaxes!(axamp1, axphase1)
    linkxaxes!(axamp2, axphase2)
    linkyaxes!(axamp1, axamp2)
    linkyaxes!(axphase1, axphase2)

    f[1:2, 3] = Legend(f, axamp1, merge=true, unique=true, tellheight=true, tellwidth=true)
    save("StationGains_vs_time.png", f)

    @info("Plotted station gains 🙆")
end

"""
    plotbandpass(h5file::String, stationnames::Vector{String3}, chanfreqvec::Vector{Float64})

Plot bandpass gains against time.
"""
function plotbandpass(h5file::String, stationnames::Vector{String3}, chanfreqvec::Vector{Float64})
    @info("Plotting bandpass gains against time...")
    fid = h5open(h5file, "r")

    chanfreqvec_ghz = chanfreqvec / 1e9 # in GHz

    b = read(fid["bandpass"]["bjonesmatrices"])

    f = Figure(size=(500, 300))
    axphase1 = Axis(f[1, 1], ylabel="Phases (°)", title="Bandpass gains (Pol1)")
    axamp1 = Axis(f[2, 1], xlabel="Frequency (GHz)", ylabel="Amplitudes")
    axphase2 = Axis(f[1, 2], title="Bandpass gains (Pol2)")
    axamp2 = Axis(f[2, 2], xlabel="Frequency (GHz)")

    # modify axis attributes
    hidexdecorations!(axphase1, grid=false, ticks=false)
    hidexdecorations!(axphase2, grid=false, ticks=false)
    hideydecorations!(axamp2, grid=false, ticks=false)
    hideydecorations!(axphase2, grid=false, ticks=false)

    for ant in eachindex(stationnames)
        lines!(axamp1, chanfreqvec_ghz, abs.(b[1, 1, :, ant]), lw=1, color=ColorSchemes.mk_15[ant], label=stationnames[ant])
        lines!(axphase1, chanfreqvec_ghz, rad2deg.(angle.(b[1, 1, :, ant])), lw=1, color=ColorSchemes.mk_15[ant], label=stationnames[ant])

        lines!(axamp2, chanfreqvec_ghz, abs.(b[2, 2, :, ant]), lw=1, color=ColorSchemes.mk_15[ant], label=stationnames[ant])
        lines!(axphase2, chanfreqvec_ghz, rad2deg.(angle.(b[2, 2, :, ant])), lw=1, color=ColorSchemes.mk_15[ant], label=stationnames[ant])
    end
    close(fid) # close HDF5 file

    linkxaxes!(axamp1, axphase1)
    linkxaxes!(axamp2, axphase2)
    linkyaxes!(axamp1, axamp2)
    linkyaxes!(axphase1, axphase2)

    f[1:2, 3] = Legend(f, axamp1, merge=true, unique=true, tellheight=true, tellwidth=true)
    save("BandpassGains_vs_frequency.png", f)

    @info("Plotted bandpass 🙆")
end

"""
    plotpointingerrors(h5file::String, scanno::Vector{Int32}, stationnames::Vector{String3})

Plot pointing errors.
"""
function plotpointingerrors(h5file::String, scanno::Vector{Int32}, stationnames::Vector{String3})
    @info("Plotting pointing errors...")
    fid = h5open(h5file, "r")

    # get unique scan numbers
    uniqscans = unique(scanno)

    f = Figure(size=(500, 300))
    axoff = Axis(f[1, 1], title="Pointing offsets", ylabel="Pointing offsets (arcsec)")
    axamperr = Axis(f[2, 1], title="Pointing amplitude errors", ylabel="Primary beam response", xlabel="Time (mispointing number)")

    # modify axis attributes
    hidexdecorations!(axoff, grid=false, ticks=false)

    indexstart = 1
    indexend = 0
    for scan in uniqscans
        offsetarr = read(fid["pointingerrors"]["perscanoffsets_scan$(scan)"])
        amperrarr = read(fid["pointingerrors"]["perscanamperrors_scan$(scan)"])
        indexend = indexend + size(offsetarr)[1] # both arrays have same dimensions

        for ant in eachindex(stationnames)
            if scan == 1
                lines!(axoff, 1:indexend, offsetarr[:,ant], lw=1, color=ColorSchemes.mk_15[ant], label=stationnames[ant])
                lines!(axamperr, 1:indexend, amperrarr[:,ant], lw=1, color=ColorSchemes.mk_15[ant], label=stationnames[ant])
            else
                lines!(axoff, indexstart:indexend, offsetarr[:,ant], lw=1, color=ColorSchemes.mk_15[ant], label=stationnames[ant])
                lines!(axamperr, indexstart:indexend, amperrarr[:,ant], lw=1, color=ColorSchemes.mk_15[ant], label=stationnames[ant])
            end
        end
        indexstart = indexend + 1
    end
    close(fid) # close HDF5 file

    linkxaxes!(axoff, axamperr)

    f[1:2, 2] = Legend(f, axoff, merge=true, unique=true, tellheight=true, tellwidth=true)
    save("Pointing_offsets_amplitude_errors_vs_time.png", f)

    @info("Plotted pointing offsets and amplitude errors 🙆")
end

"""
    plotelevationangle(h5file::String, scanno::Vector{Int32}, times::Vector{Float64}, stationnames::Vector{String3})

Plot evolution of station elevation angles during the course of the observation.
"""
function plotelevationangle(h5file::String, scanno::Vector{Int32}, times::Vector{Float64}, stationnames::Vector{String3})
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

    f = Figure(size=(500, 300))
    ax = Axis(f[1, 1], xlabel="Relative time (hr)", ylabel="Elevation angle (°)", title="Elevation angles by station")

    for ant in eachindex(stationnames)
        scatter!(ax, x, rad2deg.(elevmat[:,ant]), color=ColorSchemes.mk_15[ant], label=stationnames[ant], markersize=2)
    end
    close(fid) # close HDF5 file

    f[1, 2] = Legend(f, ax, merge=true, unique=true, tellheight=true, tellwidth=true)
    save("Elevation_angle_vs_time.png", f)

    @info("Plotted elevation angles 🙆")
end

"""
    plotparallacticangle(h5file::String, scanno::Vector{Int32}, times::Vector{Float64}, stationnames::Vector{String3})

Plot evolution of station parallactic angles during the course of the observation.
"""
function plotparallacticangle(h5file::String, scanno::Vector{Int32}, times::Vector{Float64}, stationnames::Vector{String3})
    @info("Plotting elevation angles by station...")
    
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

    f = Figure(size=(500, 300))
    ax = Axis(f[1, 1], xlabel="Relative time (hr)", ylabel="Parallactic angle (°)", title="Parallactic angles by station")

    for ant in eachindex(stationnames)
        scatter!(ax, x, rad2deg.(parangmat[:,ant]), color=ColorSchemes.mk_15[ant], label=stationnames[ant], markersize=2)
    end
    close(fid) # close HDF5 file

    f[1, 2] = Legend(f, ax, merge=true, unique=true, tellheight=true, tellwidth=true)
    save("Parallactic_angle_vs_time.png", f)

    @info("Plotted parallactic angles 🙆")
end

"""
    plotdterms(h5file::String, stationnames::Vector{String3}, chanfreqvec::Vector{Float64})

Plot frequency-dependent complex instrumental polarization.
"""
function plotdterms(h5file::String, stationnames::Vector{String3}, chanfreqvec::Vector{Float64})
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

    f = Figure(size=(500, 300))
    axphase1 = Axis(f[1, 1], ylabel="Phases (°)", title="D-terms (CrossPol1)")
    axamp1 = Axis(f[2, 1], xlabel="Frequency (GHz)", ylabel="Amplitudes")
    axphase2 = Axis(f[1, 2], title="D-terms (CrossPol2)")
    axamp2 = Axis(f[2, 2], xlabel="Frequency (GHz)")

    # modify axis attributes
    hidexdecorations!(axphase1, grid=false, ticks=false)
    hidexdecorations!(axphase2, grid=false, ticks=false)
    hideydecorations!(axamp2, grid=false, ticks=false)
    hideydecorations!(axphase2, grid=false, ticks=false)

    for ant in eachindex(stationnames)
        lines!(axamp1, chanfreqvec_ghz, abs.(d[1, 2, :, ant]), lw=1, color=ColorSchemes.mk_15[ant], label=stationnames[ant])
        lines!(axphase1, chanfreqvec_ghz, rad2deg.(angle.(d[1, 2, :, ant])), lw=1, color=ColorSchemes.mk_15[ant], label=stationnames[ant])

        lines!(axamp2, chanfreqvec_ghz, abs.(d[2, 1, :, ant]), lw=1, color=ColorSchemes.mk_15[ant], label=stationnames[ant])
        lines!(axphase2, chanfreqvec_ghz, rad2deg.(angle.(d[2, 1, :, ant])), lw=1, color=ColorSchemes.mk_15[ant], label=stationnames[ant])
    end
    close(fid) # close HDF5 file

    linkxaxes!(axamp1, axphase1)
    linkxaxes!(axamp2, axphase2)
    linkyaxes!(axamp1, axamp2)
    linkyaxes!(axphase1, axphase2)

    f[1:2, 3] = Legend(f, axamp1, merge=true, unique=true, tellheight=true, tellwidth=true)
    save("Dterms_vs_frequency.png", f)

    @info("Plotted D-terms 🙆")
end

"""
    plottransmission(h5file::String, stationnames::Vector{String3}, times::Vector{Float64}, chanfreqvec::Vector{Float64})

Plot tropospheric transmission variation with frequency.
"""
function plottransmission(h5file::String, stationnames::Vector{String3}, times::Vector{Float64}, chanfreqvec::Vector{Float64})
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
        f = Figure(size=(500, 300))
        ax = Axis(f[1, 1], xlabel="Relative time (hr)", ylabel="Transmission")
        for ant in eachindex(stationnames)
            scatter!(ax, reltimes, tr[1,:,ant], markersize=2, color=ColorSchemes.mk_15[ant], label=stationnames[ant])
        end
        f[1, 2] = Legend(f, ax, merge=true, unique=true, tellheight=true, tellwidth=true)
        save("Transmission_vs_time.png", f)
    else
        for ant in eachindex(stationnames)
            xtr = tr[:,:,ant]'
            f = Figure(size=(500, 300))
            ax = Axis(f[1, 1], xlabel="Relative time (hr)", ylabel="Frequency (GHz)", title="Transmission ($(stationnames[ant]))")
            contour!(ax, reltimes, chanfreqvec_ghz, xtr, title=stationnames[ant], levels=10, color=ColorSchemes.mk_15[ant])
            save("Transmission_$(stationnames[ant]).png", f)
        end
    end

    @info("Plotted tropospheric transmission 🙆")
end

"""
    plotmeandelays(h5file::String, stationnames::Vector{String3}, times::Vector{Float64}, chanfreqvec::Vector{Float64})

Plot mean delays against time.
"""
function plotmeandelays(h5file::String, stationnames::Vector{String3}, times::Vector{Float64}, chanfreqvec::Vector{Float64})
    @info("Plotting station-based transmission values")

    # get unique times
    uniqtimes = unique(times)
    reltimes = (uniqtimes .- first(uniqtimes))/3600.0 # in units of hours

    chanfreqvec_ghz = chanfreqvec/1e9 # in GHz

    # read in the computed mean delay values
    fid = h5open(h5file, "r")
    md = read(fid["troposphere"]["meandelays"]) * 1e9 # delays / ns
    close(fid)

    f = Figure(size=(500, 300))
    ax = Axis(f[1, 1], xlabel="Relative time (hr)", ylabel="Mean delays (ns)")

    for ant in eachindex(stationnames)
        scatter!(ax, reltimes, dropdims(mean(md[:,:,ant], dims=1), dims=1), color=ColorSchemes.mk_15[ant], label=stationnames[ant], markersize=2)
    end
    f[1, 2] = Legend(f, ax, merge=true, unique=true, tellheight=true, tellwidth=true)
    save("MeanDelays_vs_time.png", f)
    
    @info("Plotted mean delays 🙆")
end