export troposphere

include(joinpath("util.jl"))
using LazyGrids
using LinearAlgebra

const Boltzmann = 1.380649e-23
const lightspeed = 2.99792458e8

function run_atm(obs::CjlObservation)::DataFrame
    """
    Run ATM (Pardo et al. 2001) to compute absorption by and dispersive delay in the atmosphere
    """
    # loop through stations
    absdfvec = []
    dispdfvec = []
    for ant in 1:size(obs.stationinfo)[1]
	# absorption
	atmcommand = obs.numchan == 1 ? `absorption --fmin $(obs.chanfreqvec[1]/1e9-1.0) --fmax $(obs.chanfreqvec[1]/1e9) --fstep 1.0 --pwv $(obs.stationinfo.pwv_mm[ant]) 
	--gpress $(obs.stationinfo.gpress_mb[ant]) --gtemp $(obs.stationinfo.gtemp_K[ant])` : `absorption --fmin $((first(obs.chanfreqvec)-obs.chanwidth)/1e9) --fmax $(last(obs.chanfreqvec)/1e9) 
	--fstep $(obs.chanwidth/1e9) --pwv $(obs.stationinfo.pwv_mm[ant]) --gpress $(obs.stationinfo.gpress_mb[ant]) --gtemp $(obs.stationinfo.gtemp_K[ant])`
        run(pipeline(atmcommand, stdout="absorption.csv", append=false))

	# load as dataframe
	push!(absdfvec, CSV.read("absorption.csv", DataFrame; delim=",", ignorerepeated=false, header=["Frequency", "Dry_opacity", "Wet_opacity", "Sky_brightness"], skipto=2))

	# dispersive delay
	atmcommand = obs.numchan == 1 ? `dispersive --fmin $(obs.chanfreqvec[1]/1e9-1.0) --fmax $(obs.chanfreqvec[1]/1e9) --fstep 1.0 --pwv $(obs.stationinfo.pwv_mm[ant]) 
	--gpress $(obs.stationinfo.gpress_mb[ant]) --gtemp $(obs.stationinfo.gtemp_K[ant])` : `dispersive --fmin $((first(obs.chanfreqvec)-obs.chanwidth)/1e9) --fmax $(last(obs.chanfreqvec)/1e9) 
	--fstep $(obs.chanwidth/1e9) --pwv $(obs.stationinfo.pwv_mm[ant]) --gpress $(obs.stationinfo.gpress_mb[ant]) --gtemp $(obs.stationinfo.gtemp_K[ant])`
        run(pipeline(atmcommand, stdout="dispersive.csv", append=false))

	# load as dataframe
	push!(dispdfvec, CSV.read("dispersive.csv", DataFrame; delim=",", ignorerepeated=false, header=["Frequency", "Wet_nondisp", "Wet_disp", "Dry_nondisp"], skipto=2))
    end

    # concatenate all to one dataframe
    absdf = vcat(absdfvec...,source=:Station => obs.stationinfo.station)
    dispdf = vcat(dispdfvec...,source=:Station => obs.stationinfo.station)
    absdf[!, :Wet_nondisp] = dispdf.Wet_nondisp
    absdf[!, :Wet_disp] = dispdf.Wet_disp
    absdf[!, :Dry_nondisp] = dispdf.Dry_nondisp

    CSV.write("atm.csv", absdf) # write all to a single csv file
    run(`rm absorption.csv dispersive.csv`) # remove clutter

    @info("Running ATM to compute absorption by and dispersive delay in the troposphere...")
    return absdf
end

function compute_transmission(obs::CjlObservation, atmdf::DataFrame, elevationmatrix::Array{Float64, 2}, ntimes::Int64, g::HDF5.Group)::Array{Float32, 3}
    """
    Compute transmission matrix
    """
    # compute time and frequency varying transmission matrix for each station
    transmission = zeros(Float32, obs.numchan, ntimes, size(obs.stationinfo)[1])
    for ant in 1:size(obs.stationinfo)[1]
        # get opacity at all frequencies for this station
	opacityvec = obs.yamlconf["troposphere"]["wetonly"] ? atmdf[atmdf.Station .== obs.stationinfo.station[ant],:].Wet_opacity : 
	             atmdf[atmdf.Station .== obs.stationinfo.station[ant],:].Dry_opacity + atmdf[atmdf.Station .== obs.stationinfo.station[ant],:].Wet_opacity
        for t in 1:ntimes
            for chan in 1:obs.numchan
                    transmission[chan, t, ant] = exp(-1*opacityvec[chan]/sin(elevationmatrix[t, ant]))
            end
        end
    end

    if !(haskey(g, "transmission"))
        g["transmission"] = transmission
    end

    return transmission
end

function attenuate(obs::CjlObservation, atmdf::DataFrame, elevationmatrix::Array{Float64, 2}, g::HDF5.Group)
    """
    Compute signal attenuation due to opacity
    """
    # get unique times
    uniqtimes = unique(obs.times)

    transmission = compute_transmission(obs, atmdf, elevationmatrix, size(uniqtimes)[1], g)

    # attenuate visibilities
    row = 1
    for t in 1:size(transmission)[2] # no. of unique times
	# read all baselines present in a given time
	ant1vec = getindex(obs.antenna1, findall(obs.times.==uniqtimes[t]))
	ant2vec = getindex(obs.antenna2, findall(obs.times.==uniqtimes[t]))
	for (ant1,ant2) in zip(ant1vec, ant2vec)
            for chan in 1:obs.numchan
		obs.data[:, :, chan, row] = sqrt(transmission[chan, t, ant1+1]*transmission[chan, t, ant2+1]) .* abs.(obs.data[:, :, chan, row]) .* exp.(angle.(obs.data[:, :, chan, row])*im)
	    end
	    row += 1 # increment the last dimension i.e. row number
        end
    end
    @info("Apply attenuation due to opacity... 🙆")
end

function compute_skynoise(obs::CjlObservation, atmdf::DataFrame, elevation::Array{Float64, 2}, g::HDF5.Group)
    """
    Use ATM (Pardo et al. 2001) to compute atmospheric contribution to temperature
    """
    # get matrix type and size to be created -- data is a 4d array
    elemtype = typeof(obs.data[1])

    # get unique times
    uniqtimes = unique(obs.times)

    transmission = compute_transmission(obs, atmdf, elevation, size(uniqtimes)[1], g)
    sefdarray = zeros(Float64, obs.numchan, size(uniqtimes)[1], size(obs.stationinfo)[1])
    for ant in 1:size(obs.stationinfo)[1]
        skytempvec = atmdf[atmdf.Station .== obs.stationinfo.station[ant],:].Sky_brightness
        for t in 1:size(uniqtimes)[1]
            for chan in 1:obs.numchan
                sefdarray[chan, t, ant] = 2*Boltzmann/(obs.stationinfo.aperture_eff[ant]*pi*((obs.stationinfo.dishdiameter_m[ant]/2.0)^2))*(1e26*skytempvec[chan]*(1.0 - transmission[chan, t, ant]))
            end
        end
    end

    # compute sky noise and add to data
    skynoiserms = zeros(Float64, size(obs.data))
    skynoise = zeros(elemtype, size(obs.data))
    row = 1
    for t in 1:size(transmission)[2] # no. of unique times
        # read all baselines present in a given time
        ant1vec = getindex(obs.antenna1, findall(obs.times.==uniqtimes[t]))
        ant2vec = getindex(obs.antenna2, findall(obs.times.==uniqtimes[t]))
        for (ant1,ant2) in zip(ant1vec, ant2vec)
            for chan in 1:obs.numchan
		skynoiserms[:, :, chan, row] .= (1/obs.yamlconf["correff"]) * sqrt((sefdarray[chan, t, ant1+1]*sefdarray[chan, t, ant2+1])/(2*obs.exposure*obs.chanwidth))
		for jj in 1:2
		    for ii in 1:2
	                skynoise[ii, jj, chan, row] = skynoiserms[ii,jj, chan, row]*randn(obs.rngtrop, elemtype) # sky noise is polarized
		        obs.data[ii, jj, chan, row] += skynoise[ii, jj, chan, row]
		    end
	        end
            end
            row += 1 # increment the last dimension i.e. row number
        end
    end

    # write to h5 file
    if !(haskey(g, "sefdarray"))
        g["sefdarray"] = sefdarray
    end
    if !(haskey(g, "skynoiserms"))
        g["skynoiserms"] = skynoiserms
    end
    if !(haskey(g, "skynoise"))
        g["skynoise"] = skynoise
    end

    @info("Apply tropospheric noise... 🙆")
end

function compute_meandelays(obs::CjlObservation, atmdf::DataFrame, elevationmatrix::Array{Float64, 2}, g::HDF5.Group)
    """
    Use ATM (Pardo et al. 2001) to compute mean delays
    """
    # get element type to be used
    elemtype = typeof(obs.data[1][1])

    # get unique times
    uniqtimes = unique(obs.times)
    ntimes = size(uniqtimes)[1]

    # compute time and frequency varying phase delays for each station
    phasedelays = zeros(elemtype, obs.numchan, ntimes, size(obs.stationinfo)[1])
    for ant in 1:size(obs.stationinfo)[1]
        # get delta path length
	deltapathlengthvec = obs.yamlconf["troposphere"]["wetonly"] ? atmdf[atmdf.Station .== obs.stationinfo.station[ant],:].Wet_disp + 
	                     atmdf[atmdf.Station .== obs.stationinfo.station[ant],:].Wet_nondisp : atmdf[atmdf.Station .== obs.stationinfo.station[ant],:].Wet_disp + 
			     atmdf[atmdf.Station .== obs.stationinfo.station[ant],:].Wet_nondisp + atmdf[atmdf.Station .== obs.stationinfo.station[ant],:].Dry_nondisp
        for t in 1:ntimes
            for chan in 1:obs.numchan
		phasedelays[chan, t, ant] = 2*pi*(deltapathlengthvec[chan]/lightspeed/sin(elevationmatrix[t, ant]))*obs.chanfreqvec[chan]
            end
        end
    end

    if !(haskey(g, "phasedelays"))
        g["mean_phasedelays"] = phasedelays
    end

    # apply phase delays to visibilities
    row = 1
    for t in 1:size(phasedelays)[2] # no. of unique times
        # read all baselines present in a given time
        ant1vec = getindex(obs.antenna1, findall(obs.times.==uniqtimes[t]))
        ant2vec = getindex(obs.antenna2, findall(obs.times.==uniqtimes[t]))
        for (ant1,ant2) in zip(ant1vec, ant2vec)
            for chan in 1:obs.numchan
		    obs.data[:, :, chan, row] .*= exp((phasedelays[chan, t, ant1+1]-phasedelays[chan, t, ant2+1])*im)
            end
            row += 1 # increment the last dimension i.e. row number
        end
    end

    @info("Apply mean delays due to troposphere... 🙆")    

end

function compute_turbulence(obs::CjlObservation, atmdf::DataFrame, elevationmatrix::Array{Float64, 2}, g::HDF5.Group)
    """
    Use ATM (Pardo et al. 2001) to compute turbulent phases
    """
    beta::Float64 = 5/3 # power law index

    # get element type to be used
    elemtype = typeof(obs.data[1][1])

    # get unique scan numbers
    uniqscans = unique(obs.scanno)

    nant::Int64 = size(obs.stationinfo)[1] # get nant

    # loop through each scan
    row = 1 # variable to index obs.data array
    for scan in uniqscans
        # compute ideal ntimes per scan
        actualtscanvec = unique(getindex(obs.times, findall(obs.scanno.==scan)))
        actualtscanveclen = length(actualtscanvec)
        idealtscanvec = collect(first(actualtscanvec):obs.exposure:last(actualtscanvec))
        idealtscanveclen = length(idealtscanvec)

	(xg, yg) = ndgrid(1:idealtscanveclen, 1:idealtscanveclen)
	offset_idealtscanvec = idealtscanvec .- first(idealtscanvec)

        # create 4d array to hold G-Jones terms per time per station
        turbulence_phasedelays = zeros(Float64, obs.numchan, idealtscanveclen, size(obs.stationinfo)[1]) # nchan x ntimes x nant

        for ant in 1:nant
	    structD = (offset_idealtscanvec ./ obs.stationinfo.ctime_sec[ant]) .^ beta # compute structure function
	    autocorrC = abs.(0.5 .* (last(structD) .- structD)) # compute autocorrelation function, clipped at largest mode

	    # compute covariance matrix
	    indices = abs.(xg'.-yg').+1
	    covmatS = zeros(typeof(first(autocorrC)), idealtscanveclen, idealtscanveclen)
	    for ii in 1:idealtscanveclen
		for jj in 1:idealtscanveclen
		    covmatS[ii, jj] = autocorrC[indices[ii,jj]]
	        end
            end
	    L = Array(cholesky(covmatS).L) # Cholesky factorise the covariance matrix
	    timeseries = [sum(L[ii,:] .* randn(obs.rngtrop, Float64, idealtscanveclen)) for ii in 1:idealtscanveclen]

	    # populate turbulence_phasedelays
	    for t in 1:idealtscanveclen
		for chan in 1:obs.numchan
		    if chan == 1
			turbulence_phasedelays[chan, t, ant] = sqrt(1/sin(elevationmatrix[t, ant])) * timeseries[t]
		    else
			turbulence_phasedelays[chan, t, ant] = turbulence_phasedelays[1, t, ant] * obs.chanfreqvec[chan]/first(obs.chanfreqvec)
		    end
		end
	    end
        end

        # loop over time/row and apply gjones terms corresponding to each baseline
        findnearest(A,x) = argmin(abs.(A .- x)) # define function to find nearest neighbour
        for t in 1:actualtscanveclen
            idealtimeindex = findnearest(idealtscanvec, actualtscanvec[t])
            # read all baselines present in a given time
            ant1vec = getindex(obs.antenna1, findall(obs.times.==actualtscanvec[t]))
            ant2vec = getindex(obs.antenna2, findall(obs.times.==actualtscanvec[t]))
            for (ant1,ant2) in zip(ant1vec, ant2vec)
                for chan in 1:obs.numchan
		    obs.data[:,:,chan,row] .*= exp((turbulence_phasedelays[chan,idealtimeindex,ant1+1]-turbulence_phasedelays[chan,idealtimeindex,ant2+1])*im)
                end
                row += 1 # increment obs.data last index i.e. row number
            end
        end	

	g["turbulence_phasedelays_scan$(scan)"] = turbulence_phasedelays

    end

    # add datatype attribute
    attributes(g)["datatype"] = string(typeof(read(g[keys(g)[1]])))

    @info("Introduce turbulence in the troposphere... 🙆")
end

function troposphere(obs::CjlObservation)
    """
    Add tropospheric effects
    """
    @info("Computing tropospheric effects...")

    # open h5 file for writing
    fid = h5open(obs.yamlconf["hdf5corruptions"], "r+")
    g = create_group(fid, "troposphere")
    attributes(g)["desc"] = "all tropospheric signal corruptions"
    attributes(g)["dims"] = "various"
    
    @time atmdf = run_atm(obs)

    elevationmatrix = elevationangle(obs) # compute elevation angle for all stations

    if !(haskey(g, "elevation"))
        g["elevation"] = elevationmatrix
    end

    # attenuate
    obs.yamlconf["troposphere"]["attenuate"] && @time attenuate(obs, atmdf, elevationmatrix, g)

    # skynoise
    obs.yamlconf["troposphere"]["skynoise"] && @time compute_skynoise(obs, atmdf, elevationmatrix, g)

    # meandelays
    obs.yamlconf["troposphere"]["meandelays"] && @time compute_meandelays(obs, atmdf, elevationmatrix, g)

    # turbulence
    obs.yamlconf["troposphere"]["turbulence"] && @time compute_turbulence(obs, atmdf, elevationmatrix, g)

    # close h5 file
    close(fid)

    @info("Finish propagating signal through troposphere... 🙆")
end
