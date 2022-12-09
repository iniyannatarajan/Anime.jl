export troposphere

include(joinpath("util.jl"))

function run_atm(obs::CjlObservation)::DataFrame
    """
    Run ATM (Pardo et al. 2001) to compute absorption by and dispersive delay in the atmosphere
    """
    # loop through stations
    absdfvec = []
    dispdfvec = []
    for ant in 1:size(obs.stationinfo)[1]
	# absorption
	atmcommand = obs.numchan == 1 ? `absorption --freq $(obs.chanfreqvec[1]/1e9) --pwv $(obs.stationinfo.pwv_mm[ant]) --gpress $(obs.stationinfo.gpress_mb[ant]) 
	--gtemp $(obs.stationinfo.gtemp_K[ant])` : `absorption --fmin $((first(obs.chanfreqvec)-obs.chanwidth)/1e9) --fmax $(last(obs.chanfreqvec)/1e9) --fstep $(obs.chanwidth/1e9) 
	    --pwv $(obs.stationinfo.pwv_mm[ant]) --gpress $(obs.stationinfo.gpress_mb[ant]) --gtemp $(obs.stationinfo.gtemp_K[ant])`
        run(pipeline(atmcommand, stdout="absorption.csv", append=false))

	# load as dataframe
	push!(absdfvec, CSV.read("absorption.csv", DataFrame; delim=",", ignorerepeated=false, header=["Frequency", "Dry_opacity", "Wet_opacity", "Sky_brightness"], skipto=2))

	# dispersive delay
	atmcommand = obs.numchan == 1 ? `dispersive --freq $(obs.chanfreqvec[1]/1e9) --pwv $(obs.stationinfo.pwv_mm[ant]) --gpress $(obs.stationinfo.gpress_mb[ant]) 
	--gtemp $(obs.stationinfo.gtemp_K[ant])` : `dispersive --fmin $((first(obs.chanfreqvec)-obs.chanwidth)/1e9) --fmax $(last(obs.chanfreqvec)/1e9) --fstep $(obs.chanwidth/1e9) 
	    --pwv $(obs.stationinfo.pwv_mm[ant]) --gpress $(obs.stationinfo.gpress_mb[ant]) --gtemp $(obs.stationinfo.gtemp_K[ant])`
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

function attenuate(obs::CjlObservation, atmdf::DataFrame, elevationmatrix::Array{Float64, 2})
    """
    Compute signal attenuation due to opacity
    """
    # get unique times
    uniqtimes = unique(obs.times)

    # compute time and frequency varying transmission matrix for each station
    transmissionmatrix = zeros(Float32, obs.numchan, size(uniqtimes)[1], size(obs.stationinfo)[1])
    for ant in 1:size(obs.stationinfo)[1]
	# get opacity at all frequencies for this station
	opacityvec = obs.yamlconf["troposphere"]["wetonly"] ? atmdf[atmdf.Station .== "AA",:].Wet_opacity : opacity = atmdf[atmdf.Station .== "AA",:].Dry_opacity + 
			atmdf[atmdf.Station .== "AA",:].Wet_opacity
        for t in 1:size(uniqtimes)[1]
            for chan in 1:obs.numchan
		    transmissionmatrix[chan, t, ant] = exp(-1*opacityvec[chan]/sin(elevationmatrix[t, ant]))
	    end
        end
    end

    # attenuate visibilities
    row = 1
    for t in 1:size(uniqtimes)[1]
	# read all baselines present in a given time
	ant1vec = getindex(obs.antenna1, findall(obs.times.==uniqtimes[t]))
	ant2vec = getindex(obs.antenna2, findall(obs.times.==uniqtimes[t]))
	for (ant1,ant2) in zip(ant1vec, ant2vec)
            for chan in 1:obs.numchan
		obs.data[:, :, chan, row] = sqrt(transmissionmatrix[chan, t, ant1+1]*transmissionmatrix[chan, t, ant2+1]) .* abs.(obs.data[:, :, chan, row]) .* exp.(angle.(obs.data[:, :, chan, row])*im)
	    end
	    row += 1 # increment the last dimension i.e. row number
        end
    end
    @info("Applying attenuation due to opacity... 🙆")
end

function compute_skynoise()
    """
    Use ATM (Pardo et al. 2001) to compute atmospheric contribution to temperature
    """
    # needs skybrightness, opacity, elevation
end

function compute_meandelays()
    """
    Use ATM (Pardo et al. 2001) to compute mean delays
    """
    # needs elevation
end

function compute_turbulent_phases()
    """
    Use ATM (Pardo et al. 2001) to compute turbulent phases
    """
    # needs additionalpathlength, elevation
end

function troposphere(obs::CjlObservation)
    """
    Add tropospheric effects
    """
    @info("Computing tropospheric effects...")
    
    @time atmdf = run_atm(obs)

    elevationmatrix = elevationangle(obs) # compute elevation angle for all stations

    # attenuate
    obs.yamlconf["troposphere"]["attenuate"] && @time attenuate(obs, atmdf, elevationmatrix)

    @info("Finish propagating signal through troposphere... 🙆")
end
