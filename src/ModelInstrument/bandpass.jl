export bandpass

using Interpolations

"""
    bandpass(bandpassfile::String, h5file::String, data::Array{Complex{Float32},4}, stationinfo::DataFrame, rngcorrupt::AbstractRNG,
    antenna1::Vector{Int}, antenna2::Vector{Int}, numchan::Int64, chanfreqvec::Vector{Float64})

Compute the bandpass model and apply to data. The actual numerical values are serialized as HDF5.
"""
function bandpass(bandpassfile::String, h5file::String, data::Array{Complex{Float32},4}, stationinfo::DataFrame, rngcorrupt::AbstractRNG,
    antenna1::Vector{Int}, antenna2::Vector{Int}, numchan::Int64, chanfreqvec::Vector{Float64})
    # read in the station bandpass file
    bpinfo = CSV.read(bandpassfile, DataFrame; delim=",", ignorerepeated=false)

    # open h5 file for writing
    fid = h5open(h5file, "r+")
    g = create_group(fid, "bandpass")
    attributes(g)["desc"] = "Numerical values of frequency-variable per station B-Jones terms applied to data"
    attributes(g)["dims"] = "2 x 2 x nchannels x nant" #for each scan, a 4d array of 2 x 2 x nchan x nant is stored

    # create empty array of dims 2 x 2 x nchannels x nant
    bjonesmatrices = zeros(eltype(data), 2, 2, numchan, size(stationinfo)[1]) # 2 x 2 x ntimes x nant

    for station in eachindex(stationinfo.station)
	# read in bandpass info from the input bandpass file
	freqvec = bpinfo[(bpinfo[!,:station].==stationinfo.station[station]),:freq_Hz].*1e9
	pol1amp = bpinfo[(bpinfo[!,:station].==stationinfo.station[station]),:pol1_amp]
	pol2amp = bpinfo[(bpinfo[!,:station].==stationinfo.station[station]),:pol1_amp]
	pol1phaserange = Float32(bpinfo[(bpinfo[!,:station].==stationinfo.station[station]),:pol1_phaserange_deg][1]) # read only the first value
	pol2phaserange = Float32(bpinfo[(bpinfo[!,:station].==stationinfo.station[station]),:pol2_phaserange_deg][1]) # read only the first value

        # throw a warning if extrapolation is needed, by comparing freqvec and chanfreqvec
	if !(first(freqvec) <= first(chanfreqvec) <= last(freqvec)) || !(first(freqvec) <= last(chanfreqvec) <= last(freqvec))
	    @warn("Representative channel frequencies for $(stationinfo.station[station]) bandpass gains do not cover the range of MS channels. Extrapolating in some places...")
	end

	# use linear interpolation for now
	itp1 = linear_interpolation(freqvec, pol1amp, extrapolation_bc=Line())
	itp2 = linear_interpolation(freqvec, pol2amp, extrapolation_bc=Line())

	bjonespol1amp = itp1(chanfreqvec)
	bjonespol1phase = rand(rngcorrupt, Uniform(-pol1phaserange, pol1phaserange), numchan)

	bjonespol2amp = itp2(chanfreqvec)
	bjonespol2phase = rand(rngcorrupt, Uniform(-pol2phaserange, pol2phaserange), numchan)

        # populate the bjones matrices
	for chan in 1:numchan
	    bjonesmatrices[1, 1, chan, station] = bjonespol1amp[chan]*exp(deg2rad(bjonespol1phase[chan])*im)
	    bjonesmatrices[2, 2, chan, station] = bjonespol2amp[chan]*exp(deg2rad(bjonespol2phase[chan])*im)
	end
    end

    # apply bandpass to data
    for row in 1:size(data)[4]
	for chan in 1:numchan
	    data[:,:,chan,row] = bjonesmatrices[:,:,chan,antenna1[row]+1]*data[:,:,chan,row]*adjoint(bjonesmatrices[:,:,chan,antenna2[row]+1])
        end
    end

    # write bandpass gains to HDF5 file
    g["bjonesmatrices"] = bjonesmatrices
    attributes(g)["datatype"] = string(typeof(read(g[keys(g)[1]])))
    # close h5 file
    close(fid)

    @info("Compute and apply bandpass gains 🙆")
end
