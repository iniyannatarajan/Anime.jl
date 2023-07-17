export plotvisamp_vs_pbs, plotstationgains, plotbandpass, plotpointingerrors

"""
    plotvisamp_vs_pbs(data::Array{Complex{Float32},4}, flag::Array{Bool,4}, uvw::Matrix{Float64}, chanfreqvec::Array{Float64,1}, numchan::Int64)

Plot visibility amplitudes against projected baseline separation
"""
function plotvisamp_vs_pbs(data::Array{Complex{Float32},4}, flag::Array{Bool,4}, uvw::Matrix{Float64}, chanfreqvec::Array{Float64,1}, numchan::Int64; saveas="01_coherencyamplitude_vs_projectedbaseline.png")
    @info("Plotting visibility amplitudes against projected baseline separation...")
    uvwave = sqrt.(uvw[1,:].^2 .+ uvw[2,:].^2) / (299792458.0/mean(chanfreqvec)) / 1e9 # in units of Gλ

    maskindices = findall(isequal(true), flag)
    maskeddata = deepcopy(data)
    maskeddata[maskindices] .= NaN # assign NaN values to flagged visibilities

    # plot first frequency channel
    p = plot(uvwave, abs.(maskeddata[1,1,1,:]), seriestype=:scatter, ls=:dot, ms=1, mc=:red, markerstrokecolor=:red, label="RR")
    plot!(p, uvwave, abs.(maskeddata[1,2,1,:]), seriestype=:scatter, ls=:dot, ms=1, mc=:cyan, markerstrokecolor=:cyan, label="RL")
    plot!(p, uvwave, abs.(maskeddata[2,1,1,:]), seriestype=:scatter, ls=:dot, ms=1, mc=:purple, markerstrokecolor=:purple, label="LR")
    plot!(p, uvwave, abs.(maskeddata[2,2,1,:]), seriestype=:scatter, ls=:dot, ms=1, mc=:green, markerstrokecolor=:green, label="LL")
    
    # plot the rest of the frequency channels
    if numchan > 1
        plot!(p, uvwave, abs.(maskeddata[1,1,2:end,:]'), seriestype=:scatter, ls=:dot, ms=1, mc=:red, markerstrokecolor=:red, label="")
        plot!(p, uvwave, abs.(maskeddata[1,2,2:end,:]'), seriestype=:scatter, ls=:dot, ms=1, mc=:cyan, markerstrokecolor=:cyan, label="")
        plot!(p, uvwave, abs.(maskeddata[2,1,2:end,:]'), seriestype=:scatter, ls=:dot, ms=1, mc=:purple, markerstrokecolor=:purple, label="")
        plot!(p, uvwave, abs.(maskeddata[2,2,2:end,:]'), seriestype=:scatter, ls=:dot, ms=1, mc=:green, markerstrokecolor=:green, label="")
    end

    plot!(p, xlabel="Projected baseline separation (Gλ)", ylabel="Complex visibility amplitude (Jy)", legend=:outertop, legendcolumns=4)
    savefig(p, saveas)
    @info("Done 🙆")
end

"""
    plotstationgains(obs; saveas="stationgainsvstime.png")

Plot station gains against time
"""
function plotstationgains(obs; saveas="stationgainsvstime.png")
    @info("Plotting station gains against time...")
    fid = h5open(obs.yamlconf["hdf5corruptions"], "r")

    # get unique scan numbers
    uniqscans = unique(obs.scanno)

    # get unique times
    uniqtimes = unique(obs.times)
    x = uniqtimes .- first(uniqtimes)

    p = plot()
    indexstart = 1
    indexend = 0
    for scan in uniqscans
        gterms = read(fid["stationgains"]["gjones_scan$(scan)"])
        indexend = indexend + size(gterms)[3]

        for ant in eachindex(obs.stationinfo.station)
            if scan == 1
                plot!(p, x[indexstart:indexend], abs.(gterms[1,1,:,ant]), lw=1, lc=ColorSchemes.mk_15[ant], label=obs.stationinfo.station[ant])
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
    plotbandpass(obs; saveas="bandpassgains.png")

Plot bandpass gains against time
"""
function plotbandpass(obs; saveas="bandpassgains.png")
    @info("Plotting bandpass gains against time...")
    fid = h5open(obs.yamlconf["hdf5corruptions"], "r")

    b = read(fid["bandpass"]["bjonesmatrices"])

    p1 = plot()
    for ant in eachindex(obs.stationinfo.station)
        plot!(p1, obs.chanfreqvec./1e9, abs.(b[1, 1, :, ant]), lw=1, lc=ColorSchemes.mk_15[ant], label="")
    end
    plot!(p1, title="Station bandpass gain amplitudes", ylabel="Pol1 gain amp") #, legend=:outertop, legendcolumns=6)

    p2 = plot()
    for ant in eachindex(obs.stationinfo.station)
        plot!(p2, obs.chanfreqvec./1e9, abs.(b[2, 2, :, ant]), lw=1, lc=ColorSchemes.mk_15[ant], label=obs.stationinfo.station[ant])
    end
    plot!(p2, xlabel="Channel frequency (GHz)", ylabel="Pol2 gain amp") # legend=:outertop, legendcolumns=6)

    p = plot(p1, p2, layout=(2, 1))
    plot!(p, legend=:outertop, legendcolumns=6)
    savefig(p, saveas)
    close(fid)
    @info("Done 🙆")
end

"""
    plotpointingerrors(obs)

Plot pointing errors
"""
function plotpointingerrors(obs)
    @info("Plotting pointing errors...")
    fid = h5open(obs.yamlconf["hdf5corruptions"], "r")

    # get unique scan numbers
    uniqscans = unique(obs.scanno)

    # get unique times
    #uniqtimes = unique(obs.times)
    #x = uniqtimes .- first(uniqtimes)

    p_off = plot()
    p_amperr = plot()
    indexstart = 1
    indexend = 0
    for scan in uniqscans
        offsetarr = read(fid["pointingerrors"]["perscanoffsets_scan$(scan)"])
        amperrarr = read(fid["pointingerrors"]["perscanamperrors_scan$(scan)"])
        indexend = indexend + size(offsetarr)[1] # both arrays have same dimensions

        for ant in eachindex(obs.stationinfo.station)
            if scan == 1
                plot!(p_off, 1:indexend, offsetarr[:,ant], lw=1, lc=ColorSchemes.mk_15[ant], label=obs.stationinfo.station[ant])
                plot!(p_amperr, 1:indexend, amperrarr[:,ant], lw=1, lc=ColorSchemes.mk_15[ant], label=obs.stationinfo.station[ant])
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