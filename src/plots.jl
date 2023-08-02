export plotuvcov, plotvis, plotstationgains, plotbandpass, plotpointingerrors, plotelevationangle

"""
    plotuvcov(uvw::Matrix{Float64}, flagrow::Vector{Bool}, chanfreqvec::Vector{Float64}; saveprefix="test_")

Plot uv-coverage of observation
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

Generate various complex visibility plots
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
    plotstationgains(h5file::String, scanno::Vector{Int32}, times::Vector{Float64}, stationnames::Vector{String3}; saveas="stationgainsvstime.png")

Plot station gains against time
"""
function plotstationgains(h5file::String, scanno::Vector{Int32}, times::Vector{Float64}, stationnames::Vector{String3}; saveas="stationgainsvstime.png")
    @info("Plotting station gains against time...")
    fid = h5open(h5file, "r")

    # get unique scan numbers
    uniqscans = unique(scanno)

    # get unique times
    uniqtimes = unique(times)
    x = uniqtimes .- first(uniqtimes)

    p = plot()
    indexstart = 1
    indexend = 0
    for scan in uniqscans
        gterms = read(fid["stationgains"]["gjones_scan$(scan)"])
        indexend = indexend + size(gterms)[3]

        for ant in eachindex(stationnames)
            if scan == 1
                plot!(p, x[indexstart:indexend], abs.(gterms[1,1,:,ant]), lw=1, lc=ColorSchemes.mk_15[ant], label=stationnames[ant])
            else
                plot!(p, x[indexstart:indexend], abs.(gterms[1,1,:,ant]), lw=1, lc=ColorSchemes.mk_15[ant], label="")
            end
        end
        indexstart = indexend + 1
       
    end
    plot!(p, xlabel="Time offset from start of observation (s)", ylabel="Gain amplitudes", legend=:outertop, legendcolumns=6)
    savefig(p, saveas)
    close(fid)
    @info("Done 🙆")
end

"""
    plotbandpass(h5file::String, stationnames::Vector{String3}, chanfreqvec::Vector{Float64}; saveas="bandpassgains.png")

Plot bandpass gains against time
"""
function plotbandpass(h5file::String, stationnames::Vector{String3}, chanfreqvec::Vector{Float64}; saveas="bandpassgains.png")
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
    savefig(p, saveas)
    close(fid)
    @info("Done 🙆")
end

"""
    plotpointingerrors(h5file::String, scanno::Vector{Int32}, stationnames::Vector{String3})

Plot pointing errors
"""
function plotpointingerrors(h5file::String, scanno::Vector{Int32}, stationnames::Vector{String3})
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
    savefig(p_off, "pointingoffsets.png")

    plot!(p_amperr, xlabel="Time (mispointing number)", ylabel="Primary beam response", legend=:outertop, legendcolumns=6)
    savefig(p_amperr, "pointingamplitudeerrors.png")
    
    close(fid)
    @info("Done 🙆")
end

"""
    plotelevationangle(h5file::String, scanno::Vector{Int32}, times::Vector{Float64}, stationnames::Vector{String3})

Plot elevation angle by station
"""
function plotelevationangle(h5file::String, scanno::Vector{Int32}, times::Vector{Float64}, stationnames::Vector{String3})
    
    # get unique scan numbers
    uniqscans = unique(scanno)

    # get unique times
    uniqtimes = unique(times)
    x = uniqtimes .- first(uniqtimes)

    fid = h5open(h5file, "r")

    if "troposphere" in keys(fid)
        elevmat = read(fid["troposphere"]["elevation"])
    elseif "polarization" in keys(fid)
        elevmat = read(fid["polarization"]["elevation"])
    else
        close(fid)
        @error("$h5file does not contain elevation angle information. Not plotting elevation angles!")
        throw(KeyError("elevation"))
    end

    if length(stationnames) != size(elevmat)[2]
        close(fid)
        throw(BoundsError("$h5file does not match the stations in station information. Not plotting elevation angles!"))
    end

    if length(x) != size(elevmat)[1]
        close(fid)
        throw(DimensionMismatch("The number of timestamps in $h5file does not match that in the observation. Not plotting elevation angles!"))
    end

    p = plot()
    for ant in eachindex(stationnames)
        plot!(p, x, rad2deg.(elevmat[:,ant]), seriestype=:scatter, ls=:dot, ms=1, mc=ColorSchemes.mk_15[ant], msc=ColorSchemes.mk_15[ant], label=stationnames[ant])
    end

    plot!(p, xlabel="Relative Time (s)", ylabel="Elevation angle (°)", legend=:outertop, legendcolumns=6)
    savefig(p, "elevationangle.png")

    close(fid)
    @info("Done 🙆")
end

"""
"""
function plotparallacticangle()
end

# plot instrumental polarization
# plot tropospheric quantities