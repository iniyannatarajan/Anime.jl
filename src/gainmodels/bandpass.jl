export bandpass!

using Interpolations

"""
    bandpass!(data::Array{Complex{Float32},4}, stationinfo::DataFrame, bandpassinfo::DataFrame, corruptseed::Int,
    antenna1::Vector{Int32}, antenna2::Vector{Int32}, numchan::Int64, chanfreqvec::Vector{Float64}; h5file::String="")

Compute the complex bandpass model and apply to data; write model to HDF5 file.
"""
function bandpass!(data::Array{Complex{Float32},4}, stationinfo::DataFrame, bandpassinfo::DataFrame, corruptseed::Int,
    antenna1::Vector{Int32}, antenna2::Vector{Int32}, numchan::Int64, chanfreqvec::Vector{Float64}; h5file::String="")
    @info("Computing bandpass...")
    #initialize RNG with seed
    rngcorrupt = Xoshiro(corruptseed)

    # create empty array of dims 2 x 2 x nchannels x nant
    bjonesmatrices = zeros(eltype(data), 2, 2, numchan, size(stationinfo)[1]) # 2 x 2 x ntimes x nant

    for station in eachindex(stationinfo.station)
	    # read in bandpass info from the input bandpass file
	    freqvec = bandpassinfo[(bandpassinfo[!,:station].==stationinfo.station[station]),:freq_Hz].*1e9
	    pol1amp = bandpassinfo[(bandpassinfo[!,:station].==stationinfo.station[station]),:pol1_amp]
	    pol2amp = bandpassinfo[(bandpassinfo[!,:station].==stationinfo.station[station]),:pol1_amp]
	    pol1phaserange = Float32(bandpassinfo[(bandpassinfo[!,:station].==stationinfo.station[station]),:pol1_phaserange_deg][1]) # read only the first value
	    pol2phaserange = Float32(bandpassinfo[(bandpassinfo[!,:station].==stationinfo.station[station]),:pol2_phaserange_deg][1]) # read only the first value

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
    if !isempty(h5file)
        fid = h5open(h5file, "cw")
        g = create_group(fid, "bandpass")
        HDF5.attributes(g)["desc"] = "Numerical values of frequency-variable per station B-Jones terms applied to data"
        HDF5.attributes(g)["dims"] = "2 x 2 x nchannels x nant" #for each scan, a 4d array of 2 x 2 x nchan x nant is stored
        g["bjonesmatrices"] = bjonesmatrices
        HDF5.attributes(g)["datatype"] = string(typeof(read(g[keys(g)[1]])))
        close(fid)
    end

    @info("Compute and apply bandpass gains 🙆")
end

"""
    bandpass!(ms::MeasurementSet, stationinfo::DataFrame, bandpassinfo::DataFrame, obsconfig::Dict; h5file::String="")

Alias for use in pipelines.
"""
function bandpass!(ms::MeasurementSet, stationinfo::DataFrame, bandpassinfo::DataFrame, obsconfig::Dict; h5file::String="")
    bandpass!(ms.data, stationinfo, bandpassinfo, obsconfig["corruptseed"], ms.antenna1, ms.antenna2,
    ms.numchan, ms.chanfreqvec, h5file=h5file)
end