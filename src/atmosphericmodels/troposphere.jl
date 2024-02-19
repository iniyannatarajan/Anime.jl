export troposphere!

using LazyGrids

const Boltzmann = 1.380649e-23
const lightspeed = 2.99792458e8

"""
    run_aatm(obs::CjlObservation; absorptionfile::String="", dispersivefile::String="")::DataFrame

Run AATM (Bjona Nikolic; Pardo et al. 2001) to compute absorption by and dispersive delay in the troposphere. If AATM is not installed,
this function can still accept input absorption and dispersion values in a specific CSV format and populate `atm.csv`.
"""
function run_aatm(obs::CjlObservation; absorptionfile::String="", dispersivefile::String="")::DataFrame
    if absorptionfile == ""
        absorptionfile = "absorption.csv"
    end
    if dispersivefile == ""
        dispersivefile = "dispersive.csv"
    end
    # loop through stations
    absdfvec = []
    dispdfvec = []

    if !isfile(absorptionfile)
        for ant in eachindex(obs.stationinfo.station)
	        # absorption
	        atmcommand = obs.numchan == 1 ? `absorption --fmin $(obs.chanfreqvec[1]/1e9-1.0) --fmax $(obs.chanfreqvec[1]/1e9) --fstep 1.0 --pwv $(obs.stationinfo.pwv_mm[ant]) 
	        --gpress $(obs.stationinfo.gpress_mb[ant]) --gtemp $(obs.stationinfo.gtemp_K[ant])` : `absorption --fmin $((first(obs.chanfreqvec)-obs.chanwidth)/1e9) --fmax $(last(obs.chanfreqvec)/1e9) 
	        --fstep $(obs.chanwidth/1e9) --pwv $(obs.stationinfo.pwv_mm[ant]) --gpress $(obs.stationinfo.gpress_mb[ant]) --gtemp $(obs.stationinfo.gtemp_K[ant])`
            run(pipeline(atmcommand, stdout=absorptionfile, append=false))

	        # load as dataframe
	        push!(absdfvec, CSV.read(absorptionfile, DataFrame; delim=",", ignorerepeated=false, header=["Frequency", "Dry_opacity", "Wet_opacity", "Sky_brightness"], skipto=2))
        end
    else
        for ant in eachindex(obs.stationinfo.station)
            push!(absdfvec, CSV.read(absorptionfile, DataFrame; delim=",", ignorerepeated=false, header=1, skipto=(ant-1)*obs.numchan+2, limit=obs.numchan))
        end
    end

    if !isfile(dispersivefile)
        for ant in eachindex(obs.stationinfo.station)
	        # dispersive delay
	        atmcommand = obs.numchan == 1 ? `dispersive --fmin $(obs.chanfreqvec[1]/1e9-1.0) --fmax $(obs.chanfreqvec[1]/1e9) --fstep 1.0 --pwv $(obs.stationinfo.pwv_mm[ant]) 
	        --gpress $(obs.stationinfo.gpress_mb[ant]) --gtemp $(obs.stationinfo.gtemp_K[ant])` : `dispersive --fmin $((first(obs.chanfreqvec)-obs.chanwidth)/1e9) --fmax $(last(obs.chanfreqvec)/1e9) 
	        --fstep $(obs.chanwidth/1e9) --pwv $(obs.stationinfo.pwv_mm[ant]) --gpress $(obs.stationinfo.gpress_mb[ant]) --gtemp $(obs.stationinfo.gtemp_K[ant])`
            run(pipeline(atmcommand, stdout=dispersivefile, append=false))

	        # load as dataframe
	        push!(dispdfvec, CSV.read(dispersivefile, DataFrame; delim=",", ignorerepeated=false, header=["Frequency", "Wet_nondisp", "Wet_disp", "Dry_nondisp"], skipto=2))
        end
    else
        for ant in eachindex(obs.stationinfo.station)
            push!(dispdfvec, CSV.read(dispersivefile, DataFrame; delim=",", ignorerepeated=false, header=1, skipto=(ant-1)*obs.numchan+2, limit=obs.numchan))
        end
    end

    #= # create input files for bypassing casacore functions for test cases
    absdf = vcat(absdfvec...)
    CSV.write("absorption1.csv", absdf)
    dispdf = vcat(dispdfvec...)
    CSV.write("dispersive1.csv", dispdf)=#

    # concatenate all to one dataframe
    absdf = vcat(absdfvec...,source=:Station => obs.stationinfo.station)
    dispdf = vcat(dispdfvec...,source=:Station => obs.stationinfo.station)
    absdf[!, :Wet_nondisp] = dispdf.Wet_nondisp
    absdf[!, :Wet_disp] = dispdf.Wet_disp
    absdf[!, :Dry_nondisp] = dispdf.Dry_nondisp

    CSV.write("atm.csv", absdf) # write all to a single csv file

    @info("Compute absorption by and dispersive delay in the troposphere using ATM 🙆")
    return absdf
end

"""
    compute_transmission!(transmission::Array{Float64, 3}, obs::CjlObservation, atmdf::DataFrame, elevationmatrix::Array{Float64, 2}, g::HDF5.Group)::Array{Float64, 3}

Compute elevation-dependent (mean) tropospheric transmission given opacity τ and elevation angle θ for each station.
```math
e^{-τ/\\sin{\\theta_{\\rm el}}}
```
"""
function compute_transmission!(transmission::Array{Float64, 3}, obs::CjlObservation, atmdf::DataFrame, elevationmatrix::Array{Float64, 2}, g::HDF5.Group)::Array{Float64, 3}
    # compute time and frequency varying transmission matrix for each station
    for ant in 1:size(transmission)[3]
        # get opacity at all frequencies for this station
	opacityvec = obs.tropwetonly ? atmdf[atmdf.Station .== obs.stationinfo.station[ant],:].Wet_opacity : 
	             atmdf[atmdf.Station .== obs.stationinfo.station[ant],:].Dry_opacity + atmdf[atmdf.Station .== obs.stationinfo.station[ant],:].Wet_opacity
        for t in 1:size(transmission)[2]
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

"""
    attenuate!(obs::CjlObservation, transmission::Array{Float64, 3})

Attenuate the signal as it passes through (mean) troposphere using precomputed transmission values.
```math
I = I_0 e^{-τ/\\sin{\\theta_{\\rm el}}}
```
"""
function attenuate!(obs::CjlObservation, transmission::Array{Float64, 3})
    # get unique times
    uniqtimes = unique(obs.times)
    ntimes = size(uniqtimes)[1]

    # attenuate visibilities
    row = 1
    for t in 1:ntimes # no. of unique times
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
    @info("Apply attenuation due to opacity 🙆")
end

"""
    compute_skynoise!(obs::CjlObservation, atmdf::DataFrame, transmission::Array{Float64, 3}, g::HDF5.Group)

Compute sky contribution to visibility noise using the radiometer equation.
"""
function compute_skynoise!(obs::CjlObservation, atmdf::DataFrame, transmission::Array{Float64, 3}, g::HDF5.Group)
    # get unique times
    uniqtimes = unique(obs.times)
    ntimes = size(uniqtimes)[1]

    # get no of stations
    nant = size(obs.stationinfo)[1]

    sefdarray = zeros(Float64, obs.numchan, ntimes, nant)
    for ant in 1:nant
        skytempvec = atmdf[atmdf.Station .== obs.stationinfo.station[ant],:].Sky_brightness
        for t in 1:ntimes
            for chan in 1:obs.numchan
                sefdarray[chan, t, ant] = 2*Boltzmann/(obs.stationinfo.aperture_eff[ant]*pi*((obs.stationinfo.dishdiameter_m[ant]/2.0)^2))*(1e26*skytempvec[chan]*(1.0 - transmission[chan, t, ant]))
            end
        end
    end

    # compute sky noise and add to data
    skynoiserms = zeros(Float64, size(obs.data))
    skynoise = zeros(eltype(obs.data), size(obs.data))
    row = 1
    for t in 1:ntimes # no. of unique times
        # read all baselines present in a given time
        ant1vec = getindex(obs.antenna1, findall(obs.times.==uniqtimes[t]))
        ant2vec = getindex(obs.antenna2, findall(obs.times.==uniqtimes[t]))
        for (ant1,ant2) in zip(ant1vec, ant2vec)
            for chan in 1:obs.numchan
		            skynoiserms[:, :, chan, row] .= sigmaperchantimebl = (1/obs.correff) * sqrt((sefdarray[chan, t, ant1+1]*sefdarray[chan, t, ant2+1])/(2*obs.exposure*obs.chanwidth))
	                skynoise[:, :, chan, row] = sigmaperchantimebl*randn(obs.rngtrop, eltype(obs.data), 2, 2) # sky noise is polarized
		        obs.data[:, :, chan, row] += skynoise[:, :, chan, row]
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

    @info("Apply tropospheric noise 🙆")
end

"""
    compute_meandelays!(obs::CjlObservation, atmdf::DataFrame, elevationmatrix::Array{Float64, 2}, g::HDF5.Group)

Compute delays due to (mean) troposphere.
"""
function compute_meandelays!(obs::CjlObservation, atmdf::DataFrame, elevationmatrix::Array{Float64, 2}, g::HDF5.Group)
    # get unique times
    uniqtimes = unique(obs.times)
    ntimes = size(uniqtimes)[1]

    # compute time and frequency varying phase delays for each station
    meandelays = zeros(Float32, obs.numchan, ntimes, size(obs.stationinfo)[1])
    phasedelays = zeros(eltype(obs.data), obs.numchan, ntimes, size(obs.stationinfo)[1])
    for ant in 1:size(obs.stationinfo)[1]
        # get delta path length
	    deltapathlengthvec = obs.tropwetonly ? atmdf[atmdf.Station .== obs.stationinfo.station[ant],:].Wet_disp + 
	            atmdf[atmdf.Station .== obs.stationinfo.station[ant],:].Wet_nondisp : atmdf[atmdf.Station .== obs.stationinfo.station[ant],:].Wet_disp + 
			    atmdf[atmdf.Station .== obs.stationinfo.station[ant],:].Wet_nondisp + atmdf[atmdf.Station .== obs.stationinfo.station[ant],:].Dry_nondisp
        
        for t in 1:ntimes
            for chan in 1:obs.numchan
                meandelays[chan, t, ant] = deltapathlengthvec[chan]/lightspeed/sin(elevationmatrix[t, ant])
            end
        end

        for t in 1:ntimes
            for chan in 1:obs.numchan
		        phasedelays[chan, t, ant] = 2*pi*meandelays[chan, t, ant]*obs.chanfreqvec[chan]
            end
        end
    end

    if !(haskey(g, "meandelays"))
        g["meandelays"] = meandelays
    end

    if !(haskey(g, "mean_phasedelays"))
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

    @info("Apply mean delays due to troposphere 🙆")    

end

"""
    compute_turbulence!(obs::CjlObservation, atmdf::DataFrame, elevationmatrix::Array{Float64, 2}, g::HDF5.Group)

Compute phase delays due to tropospheric turbulence. The time series of phase errors for station p is given by
```math
\\{{\\delta \\phi_p(t, \\nu)}\\} = \\frac{1}{\\sqrt{\\sin({\\theta_{\\mathrm{el}}(t)})}} \\{\\delta\\phi^{\\prime}_p(t)\\} \\big(\\frac{\\nu}{\\nu_0}\\big)
```
where ν is the list of channel frequencies and ν_0 is the reference frequency (lowest in the band).
"""
function compute_turbulence!(obs::CjlObservation, atmdf::DataFrame, elevationmatrix::Array{Float64, 2}, g::HDF5.Group)
    beta::Float64 = 5/3 # power law index

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
            timeseries = L * randn(obs.rngtrop, Float64, idealtscanveclen)
            
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
        turbulence_phasedelays = replace(turbulence_phasedelays, Inf => 0) # replace Inf with 0; 1/sin(0) = Inf
        
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
        
        if !haskey(g, "turbulence_phasedelays_scan$(scan)")
            g["turbulence_phasedelays_scan$(scan)"] = turbulence_phasedelays
        end
    end
    
    # add datatype attribute
    if !haskey(HDF5.attributes(g), "datatype")
        HDF5.attributes(g)["datatype"] = string(typeof(read(g[keys(g)[1]])))
    end
    
    @info("Introduce turbulence in the troposphere 🙆")
end

"""
    troposphere!(obs::CjlObservation, h5file::String; absorptionfile::String="", dispersivefile::String="", elevfile::String="")

Main function to compute various components of the tropospheric model. The actual numerical values generated are serialized in HDF5 format.
"""
function troposphere!(obs::CjlObservation, h5file::String; absorptionfile::String="", dispersivefile::String="", elevfile::String="")
    @info("Computing tropospheric effects...")

    # open h5 file for writing
    if !isempty(h5file)
        fid = h5open(h5file, "cw")
        if !haskey(fid, "troposphere")
            g = create_group(fid, "troposphere")
            HDF5.attributes(g)["desc"] = "all tropospheric signal corruptions"
            HDF5.attributes(g)["dims"] = "various"
        else
            g = fid["troposphere"]
        end
    end
    
    atmdf = run_aatm(obs, absorptionfile=absorptionfile, dispersivefile=dispersivefile) # compute necessary atmospheric quantities using atm

    if elevfile != "" && isfile(elevfile)
        fidelev = h5open(elevfile, "r")
        elevationmatrix = read(fidelev["troposphere"]["elevation"])
        close(fidelev)
    else
        elevationmatrix = elevationangle(obs.times, obs.phasedir, obs.stationinfo, obs.pos) # compute elevation angle for all stations
    end

    if !isempty(h5file)
        if !(haskey(g, "elevation"))
            g["elevation"] = elevationmatrix
        end
    end

    transmission = nothing
    if obs.tropattenuate || obs.tropskynoise
        transmission = zeros(Float64, obs.numchan, size(unique(obs.times))[1], size(obs.stationinfo)[1])
        transmission = compute_transmission!(transmission, obs, atmdf, elevationmatrix, g)
    end

    # attenuate
    obs.tropattenuate && attenuate!(obs, transmission)

    # skynoise
    obs.tropskynoise && compute_skynoise!(obs, atmdf, transmission, g)

    # meandelays
    obs.tropmeandelays && compute_meandelays!(obs, atmdf, elevationmatrix, g)

    # turbulence
    obs.tropturbulence && compute_turbulence!(obs, atmdf, elevationmatrix, g)

    # close h5 file
    close(fid)

    @info("Finish propagating signal through troposphere 🙆")
end
