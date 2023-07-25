#=@testset "PA and Elevation" begin
    y = YAML.load_file("data/testconfig.yaml", dicttype=Dict{String,Any}) # sample dict to test loadms()

    obs = loadms(y["msname"], y["stations"], Int(y["corruptseed"]), Int(y["troposphere"]["tropseed"]), y["troposphere"]["wetonly"], y["correff"], 
    y["troposphere"]["attenuate"], y["troposphere"]["skynoise"], y["troposphere"]["meandelays"], y["troposphere"]["turbulence"], 
    y["instrumentalpol"]["visibilityframe"], y["instrumentalpol"]["mode"], y["pointing"]["interval"], y["pointing"]["mode"], y["stationgains"]["mode"], 
    y["bandpass"]["bandpassfile"], delim=",", ignorerepeated=false)

    @inferred parallacticangle(obs.times, obs.phasedir, obs.stationinfo, obs.pos)
    @inferred elevationangle(obs.times, obs.phasedir, obs.stationinfo, obs.pos)
end=#

@testset "Time series" begin
    @inferred gentimeseries!(zeros(ComplexF32, 100), "gp", ComplexF32(0.0+0.0*im), Float32(1.0), Float32(3.0), 100, Xoshiro(42))
    @inferred gentimeseries!(zeros(ComplexF32, 100), "normal", ComplexF32(0.0+0.0*im), Float32(1.0), Float32(3.0), 100, Xoshiro(42))

    @inferred gentimeseries!(zeros(Float32, 100), "gp", Float32(0.0), Float32(1.0), Float32(3.0), 100, Xoshiro(42))
    @inferred gentimeseries!(zeros(Float32, 100), "normal", Float32(0.0), Float32(1.0), Float32(3.0), 100, Xoshiro(42))

    @inferred gentimeseries!(zeros(Float64, 100), "gp", 0.0, 1.0, 3.0, 100, Xoshiro(42))
    @inferred gentimeseries!(zeros(Float64, 100), "normal", 0.0, 1.0, 3.0, 100, Xoshiro(42))
end

@testset "Plots" begin
    y = YAML.load_file("data/testconfig.yaml", dicttype=Dict{String,Any}) # sample dict to test loadms()
    h5file = "data/insmodel.h5"

    obs = loadms(y["msname"], y["stations"], Int(y["corruptseed"]), Int(y["troposphere"]["tropseed"]), y["troposphere"]["wetonly"], y["correff"], 
    y["troposphere"]["attenuate"], y["troposphere"]["skynoise"], y["troposphere"]["meandelays"], y["troposphere"]["turbulence"], 
    y["instrumentalpol"]["visibilityframe"], y["instrumentalpol"]["mode"], y["pointing"]["interval"], y["pointing"]["mode"], y["stationgains"]["mode"], 
    y["bandpass"]["bandpassfile"], delim=",", ignorerepeated=false)

    @inferred plotvis(obs.uvw, obs.chanfreqvec, obs.flag, obs.data, obs.numchan, obs.times, saveprefix="test_")
    @inferred plotstationgains(h5file, obs.scanno, obs.times, obs.stationinfo.station)
    @inferred plotbandpass(h5file, obs.stationinfo.station, obs.chanfreqvec)
    @inferred plotpointingerrors(h5file, obs.scanno, obs.stationinfo.station)
end
