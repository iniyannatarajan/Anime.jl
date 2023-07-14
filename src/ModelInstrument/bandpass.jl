export bandpass

using Interpolations

"""
    bandpass(obs::CjlObservation)

Compute the bandpass model and apply to data. The actual numerical values are serialized as HDF5.
"""
function bandpass(obs::CjlObservation)
    # get element type to be used
    elemtype = typeof(obs.data[1][1])

    # read in the station bandpass file
    bpinfo = CSV.read(obs.yamlconf["bandpass"]["bandpassfile"], DataFrame; delim=",", ignorerepeated=false)

    # open h5 file for writing
    fid = h5open(obs.yamlconf["hdf5corruptions"], "r+")
    g = create_group(fid, "bandpass")
    attributes(g)["desc"] = "Numerical values of frequency-variable per station B-Jones terms applied to data"
    attributes(g)["dims"] = "2 x 2 x nchannels x nant" #for each scan, a 4d array of 2 x 2 x nchan x nant is stored

    # create empty array of dims 2 x 2 x nchannels x nant
    bjonesmatrices = zeros(elemtype, 2, 2, obs.numchan, size(obs.stationinfo)[1]) # 2 x 2 x ntimes x nant

    for station in eachindex(obs.stationinfo.station)
	# read in bandpass info from the input bandpass file
	freqvec = bpinfo[(bpinfo[!,:station].==obs.stationinfo.station[station]),:freq_Hz].*1e9
	pol1amp = bpinfo[(bpinfo[!,:station].==obs.stationinfo.station[station]),:pol1_amp]
	pol2amp = bpinfo[(bpinfo[!,:station].==obs.stationinfo.station[station]),:pol1_amp]
	pol1phaserange = Float32(bpinfo[(bpinfo[!,:station].==obs.stationinfo.station[station]),:pol1_phaserange_deg][1]) # read only the first value
	pol2phaserange = Float32(bpinfo[(bpinfo[!,:station].==obs.stationinfo.station[station]),:pol2_phaserange_deg][1]) # read only the first value

        # throw a warning if extrapolation is needed, by comparing freqvec and obs.chanfreqvec
	if !(first(freqvec) <= first(obs.chanfreqvec) <= last(freqvec)) || !(first(freqvec) <= last(obs.chanfreqvec) <= last(freqvec))
	    @warn("Representative channel frequencies for $(obs.stationinfo.station[station]) bandpass gains do not cover the range of MS channels. Extrapolating in some places...")
	end

	# use linear interpolation for now
	itp1 = linear_interpolation(freqvec, pol1amp, extrapolation_bc=Line())
	itp2 = linear_interpolation(freqvec, pol2amp, extrapolation_bc=Line())

	bjonespol1amp = itp1(obs.chanfreqvec)
	bjonespol1phase = rand(obs.rngcorrupt, Uniform(-pol1phaserange, pol1phaserange), obs.numchan)

	bjonespol2amp = itp2(obs.chanfreqvec)
	bjonespol2phase = rand(obs.rngcorrupt, Uniform(-pol2phaserange, pol2phaserange), obs.numchan)

        # populate the bjones matrices
	for chan in 1:obs.numchan
	    bjonesmatrices[1, 1, chan, station] = bjonespol1amp[chan]*exp(deg2rad(bjonespol1phase[chan])*im)
	    bjonesmatrices[2, 2, chan, station] = bjonespol2amp[chan]*exp(deg2rad(bjonespol2phase[chan])*im)
	end
    end

    # apply bandpass to data
    for row in 1:size(obs.data)[4]
	for chan in 1:obs.numchan
	    obs.data[:,:,chan,row] = bjonesmatrices[:,:,chan,obs.antenna1[row]+1]*obs.data[:,:,chan,row]*adjoint(bjonesmatrices[:,:,chan,obs.antenna2[row]+1])
        end
    end

    # write bandpass gains to HDF5 file
    g["bjonesmatrices"] = bjonesmatrices
    attributes(g)["datatype"] = string(typeof(read(g[keys(g)[1]])))
    # close h5 file
    close(fid)

    @info("Compute and apply bandpass gains 🙆")
end
