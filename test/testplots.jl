@testset "Plots" begin
    obsconfig1 = readobsconfig("data/config1.yaml")
    h5file = "data/insmodel1.h5"

    ms1 = readms(obsconfig1["msname"])

    obsconfig2 = readobsconfig("data/configuvf.yaml")
    h5file2 = "data/insmodeluvf.h5"

    ms2 = readms(obsconfig2["msname"])

    stationinfo = readstationinfo(obsconfig1["stations"])

    @inferred plotuvcov(ms1.uvw, ms1.flagrow, ms1.chanfreqvec)
    rm("test_uvcoverage.png")

    @inferred plotvis(ms1.uvw, ms1.chanfreqvec, ms1.flag, ms1.data, ms1.numchan, ms1.times, plotphases=true, saveprefix="test")
    rm("test_visibilityphase_vs_baseline.png")
    rm("test_visibilityamplitude_vs_baseline.png")
    rm("test_visibilityphase_vs_time.png")
    rm("test_visibilityamplitude_vs_time.png")

    @inferred plotstationgains(h5file, ms1.scanno, ms1.times, ms1.exposure, stationinfo.station)
    rm("StationGains_vs_time.png")

    @inferred plotstationgains(h5file2, ms2.scanno, ms2.times, ms2.exposure, stationinfo.station)
    rm("StationGains_vs_time.png")

    @inferred plotbandpass(h5file, stationinfo.station, ms1.chanfreqvec)
    rm("BandpassGains_vs_frequency.png")

    @inferred plotpointingerrors(h5file, ms1.scanno, stationinfo.station)
    rm("Pointing_offsets_amplitude_errors_vs_time.png")

    # test elevation angle plotting
    @inferred plotelevationangle(h5file, ms1.scanno, ms1.times, stationinfo.station)
    rm("ElevationAngle_vs_time.png")

    fid = h5open("testing.h5", "w")
    close(fid)
    @test_throws KeyError plotelevationangle("testing.h5", ms1.scanno, ms1.times, stationinfo.station)
    rm("testing.h5")

    st = vcat(stationinfo.station, stationinfo.station)
    @test_throws BoundsError plotelevationangle(h5file, ms1.scanno, ms1.times, st)

    ts = deepcopy(ms1.times)
    push!(ts, ts[end]+(ts[end]-ts[begin]))
    @test_throws DimensionMismatch plotelevationangle(h5file, ms1.scanno, ts, stationinfo.station)

    # test parallactic angle plotting
    @inferred plotparallacticangle(h5file, ms1.scanno, ms1.times, stationinfo.station)
    rm("ParallacticAngle_vs_time.png")

    fid = h5open("testing.h5", "w")
    close(fid)
    @test_throws KeyError plotparallacticangle("testing.h5", ms1.scanno, ms1.times, stationinfo.station)
    rm("testing.h5")

    st = vcat(stationinfo.station, stationinfo.station)
    @test_throws BoundsError plotparallacticangle(h5file, ms1.scanno, ms1.times, st)

    ts = deepcopy(ms1.times)
    push!(ts, ts[end]+(ts[end]-ts[begin]))
    @test_throws DimensionMismatch plotparallacticangle(h5file, ms1.scanno, ts, stationinfo.station)

    # test d-terms plotting
    @inferred plotdterms(h5file, stationinfo.station, ms1.chanfreqvec)
    rm("Dterms_vs_frequency.png")

    # test transmission plotting
    @inferred plottransmission(h5file, stationinfo.station, ms1.times, ms1.chanfreqvec)
    for file in glob("Transmission_*.png")
        rm(file)
    end

    @inferred plottransmission(h5file2, stationinfo.station, ms2.times, ms2.chanfreqvec)
    for file in glob("Transmission_*.png")
        rm(file)
    end

    @inferred plotmeandelays(h5file, stationinfo.station, ms1.times, ms1.chanfreqvec)
    rm("MeanDelays_vs_time.png")
end
