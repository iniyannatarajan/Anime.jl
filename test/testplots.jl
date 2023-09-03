@testset "Plots" begin
    y = YAML.load_file("data/config1.yaml", dicttype=Dict{String,Any}) # sample dict to test loadms()
    h5file = "data/insmodel1.h5"

    obs = loadms(y["msname"], y["stations"], Int(y["corruptseed"]), Int(y["troposphere"]["tropseed"]), y["troposphere"]["wetonly"], y["correff"], 
    y["troposphere"]["attenuate"], y["troposphere"]["skynoise"], y["troposphere"]["meandelays"], y["troposphere"]["turbulence"], 
    y["instrumentalpolarization"]["visibilityframe"], y["instrumentalpolarization"]["mode"], y["pointing"]["interval"], y["pointing"]["scale"], y["pointing"]["mode"], y["stationgains"]["mode"], 
    y["bandpass"]["bandpassfile"], delim=",", ignorerepeated=false)

    y2 = YAML.load_file("data/configuvf.yaml", dicttype=Dict{String,Any}) # sample dict to test loadms()
    h5file2 = "data/insmodeluvf.h5"

    obs2 = loadms(y2["msname"], y2["stations"], Int(y2["corruptseed"]), Int(y2["troposphere"]["tropseed"]), y2["troposphere"]["wetonly"], y2["correff"], 
    y2["troposphere"]["attenuate"], y2["troposphere"]["skynoise"], y2["troposphere"]["meandelays"], y2["troposphere"]["turbulence"], 
    y2["instrumentalpolarization"]["visibilityframe"], y2["instrumentalpolarization"]["mode"], y2["pointing"]["interval"], y["pointing"]["scale"], y2["pointing"]["mode"], y2["stationgains"]["mode"], 
    y2["bandpass"]["bandpassfile"], delim=",", ignorerepeated=false)

    @inferred plotuvcov(obs.uvw, obs.flagrow, obs.chanfreqvec)
    rm("test_uvcoverage.png")

    @inferred plotvis(obs.uvw, obs.chanfreqvec, obs.flag, obs.data, obs.numchan, obs.times, plotphases=true, saveprefix="test_")
    rm("test_visampvspbs.png")
    rm("test_visphasevspbs.png")
    rm("test_visampvstime.png")
    rm("test_visphasevstime.png")

    @inferred plotstationgains(h5file, obs.scanno, obs.times, obs.exposure, obs.stationinfo.station)
    rm("gains_vs_time.png")

    @inferred plotstationgains(h5file2, obs2.scanno, obs2.times, obs2.exposure, obs2.stationinfo.station)
    rm("gains_vs_time.png")

    @inferred plotbandpass(h5file, obs.stationinfo.station, obs.chanfreqvec)
    rm("bpamplitudes_vs_frequency.png")

    @inferred plotpointingerrors(h5file, obs.scanno, obs.stationinfo.station)
    rm("pointingoffsets.png")
    rm("pointingamplitudeerrors.png")

    # test elevation angle plotting
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

    # test parallactic angle plotting
    @inferred plotparallacticangle(h5file, obs.scanno, obs.times, obs.stationinfo.station)
    rm("parallacticangle.png")

    fid = h5open("testing.h5", "w")
    close(fid)
    @test_throws KeyError plotparallacticangle("testing.h5", obs.scanno, obs.times, obs.stationinfo.station)
    rm("testing.h5")

    st = vcat(obs.stationinfo.station, obs.stationinfo.station)
    @test_throws BoundsError plotparallacticangle(h5file, obs.scanno, obs.times, st)

    ts = deepcopy(obs.times)
    push!(ts, ts[end]+(ts[end]-ts[begin]))
    @test_throws DimensionMismatch plotparallacticangle(h5file, obs.scanno, ts, obs.stationinfo.station)

    # test d-terms plotting
    @inferred plotdterms(h5file, obs.stationinfo.station, obs.chanfreqvec)
    rm("dterms.png")

    # test transmission plotting
    @inferred plottransmission(h5file, obs.stationinfo.station, obs.times, obs.chanfreqvec)
    rm("transmission.png")

    @inferred plottransmission(h5file2, obs2.stationinfo.station, obs2.times, obs2.chanfreqvec)
    rm("transmission.png")

    @inferred plotmeandelays(h5file, obs.stationinfo.station, obs.times, obs.chanfreqvec)
    rm("meandelays.png")
end