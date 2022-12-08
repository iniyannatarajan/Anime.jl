export troposphere

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
        run(pipeline(atmcommand, stdout="absorption.csv"))

	# load as dataframe
	push!(absdfvec, CSV.read("absorption.csv", DataFrame; delim=",", ignorerepeated=false, header=["Frequency", "Dry_opacity", "Wet_opacity", "Sky_brightness"], skipto=2))

	# dispersive delay
	atmcommand = obs.numchan == 1 ? `dispersive --freq $(obs.chanfreqvec[1]/1e9) --pwv $(obs.stationinfo.pwv_mm[ant]) --gpress $(obs.stationinfo.gpress_mb[ant]) 
	--gtemp $(obs.stationinfo.gtemp_K[ant])` : `dispersive --fmin $((first(obs.chanfreqvec)-obs.chanwidth)/1e9) --fmax $(last(obs.chanfreqvec)/1e9) --fstep $(obs.chanwidth/1e9) 
	    --pwv $(obs.stationinfo.pwv_mm[ant]) --gpress $(obs.stationinfo.gpress_mb[ant]) --gtemp $(obs.stationinfo.gtemp_K[ant])`
        run(pipeline(atmcommand, stdout="dispersive.csv"))

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
    return absdf
end

function compute_skynoise()
    """
    Use ATM (Pardo et al. 2001) to compute atmospheric contribution to temperature
    """
    # needs skybrightness, opacity, elevation
end

function compute_opacity()
    """
    Use ATM (Pardo et al. 2001) to compute opacity
    """
    # compute and return opacity
end

function compute_additional_pathlength()
    """
    Use ATM (Pardo et al. 2001) to compute opacity
    """
    # compute and return opacity
end

function attenuate(obs::CjlObservation)
    """
    Compute signal attenuation due to opacity
    """
    # needs opacity and elevation 
end

function troposphere(obs::CjlObservation)
    """
    Add tropospheric effects
    """
    @info("Computing tropospheric effects...")
    
    @time atmdf = run_atm(obs)

    # compute opacity and extra path length
    opacity = compute_opacity()
    extrapathlength = compute_additional_pathlength()

    obs.yamlconf["troposphere"]["attenuate"] && @time attenuate()

    @info("Apply tropospheric effects to visibilities... 🙆")
end

