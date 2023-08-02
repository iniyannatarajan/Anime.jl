@testset "Plots" begin
    y = YAML.load_file("data/config1.yaml", dicttype=Dict{String,Any}) # sample dict to test loadms()
    h5file = "data/insmodel1.h5"

    obs = loadms(y["msname"], y["stations"], Int(y["corruptseed"]), Int(y["troposphere"]["tropseed"]), y["troposphere"]["wetonly"], y["correff"], 
    y["troposphere"]["attenuate"], y["troposphere"]["skynoise"], y["troposphere"]["meandelays"], y["troposphere"]["turbulence"], 
    y["instrumentalpol"]["visibilityframe"], y["instrumentalpol"]["mode"], y["pointing"]["interval"], y["pointing"]["mode"], y["stationgains"]["mode"], 
    y["bandpass"]["bandpassfile"], delim=",", ignorerepeated=false)

    @inferred plotuvcov(obs.uvw, obs.flagrow, obs.chanfreqvec)
    rm("test_uvcoverage.png")

    @inferred plotvis(obs.uvw, obs.chanfreqvec, obs.flag, obs.data, obs.numchan, obs.times, plotphases=true, saveprefix="test_")
    rm("test_visampvspbs.png")
    rm("test_visphasevspbs.png")
    rm("test_visampvstime.png")
    rm("test_visphasevstime.png")

    @inferred plotstationgains(h5file, obs.scanno, obs.times, obs.stationinfo.station)
    rm("stationgainsvstime.png")

    @inferred plotbandpass(h5file, obs.stationinfo.station, obs.chanfreqvec)
    rm("bandpassgains.png")

    @inferred plotpointingerrors(h5file, obs.scanno, obs.stationinfo.station)
    rm("pointingoffsets.png")
    rm("pointingamplitudeerrors.png")

    @inferred plotelevationangle(h5file, obs.scanno, obs.times, obs.stationinfo.station)
    rm("elevationangle.png")

    fid = h5open("testing.h5", "w")
    close(fid)
    @test_throws KeyError plotelevationangle("testing.h5", obs.scanno, obs.times, obs.stationinfo.station)
    rm("testing.h5")

    st = vcat(obs.stationinfo.station, obs.stationinfo.station)
    @test_throws BoundsError plotelevationangle(h5file, obs.scanno, obs.times, st)

    ts = deepcopy(obs.times)
    push!(ts, ts[end]+(ts[end]-ts[begin]))
    @test_throws DimensionMismatch plotelevationangle(h5file, obs.scanno, ts, obs.stationinfo.station)
end
