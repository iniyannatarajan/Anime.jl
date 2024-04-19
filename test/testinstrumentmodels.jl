@testset "Troposphere" begin
    obsconfig = readobsconfig("data/config1.yaml")
    stationinfo = readstationinfo(obsconfig["stations"], delim=",", ignorerepeated=false)
    h5file = "tropos.h5"

    ms = readms(obsconfig["msname"])

    @inferred troposphere!(ms, stationinfo, obsconfig, h5file, absorptionfile="data/absorption1.csv", dispersivefile="data/dispersive1.csv", elevfile="data/insmodel1.h5")

    rm("atm.csv")
    rm(h5file)
end

@testset "Polarization" begin
    obsconfig = readobsconfig("data/config1.yaml")
    stationinfo = readstationinfo(obsconfig["stations"], delim=",", ignorerepeated=false)
    h5file = "inspol.h5"

    ms = readms(obsconfig["msname"])

    @inferred instrumentalpolarization!(ms, stationinfo, obsconfig, h5file=h5file)
    rm(h5file)

    @inferred instrumentalpolarization!(ms, stationinfo, obsconfig, h5file=h5file, elevfile="data/insmodel1.h5", parangfile="data/insmodel1.h5")
    rm(h5file)

    # write vis in antenna frame
    ms = readms(obsconfig["msname"])

    @inferred instrumentalpolarization!(ms, stationinfo, obsconfig, h5file=h5file, elevfile="data/insmodel1.h5", parangfile="data/insmodel1.h5")
    rm(h5file)
end

@testset "Primary Beam" begin
    obsconfig = readobsconfig("data/config1.yaml")
    stationinfo = readstationinfo(obsconfig["stations"], delim=",", ignorerepeated=false)
    h5file = "beam.h5"

    ms = readms(obsconfig["msname"])

    @inferred pointing!(ms, stationinfo, obsconfig, h5file=h5file)
    rm(h5file)

    # set pointing interval to 0.0
    obsconfig["pointing"]["interval"] = 0.0

    @inferred pointing!(ms, stationinfo, obsconfig, h5file=h5file)
    rm(h5file)
end

@testset "Station Gains" begin
    obsconfig = readobsconfig("data/config1.yaml")
    stationinfo = readstationinfo(obsconfig["stations"], delim=",", ignorerepeated=false)
    h5file = "gains.h5"

    ms = readms(obsconfig["msname"])

    @inferred stationgains!(ms, stationinfo, obsconfig, h5file=h5file)

    rm(h5file)
end

@testset "Bandpass" begin
    obsconfig = readobsconfig("data/config1.yaml")
    stationinfo = readstationinfo(obsconfig["stations"], delim=",", ignorerepeated=false)
    bandpassinfo = readbandpassinfo(obsconfig["bandpass"]["bandpassfile"], delim=",", ignorerepeated=false)
    h5file = "bandpass.h5"

    ms = readms(obsconfig["msname"])

    @inferred bandpass!(ms, stationinfo, bandpassinfo, obsconfig, h5file=h5file)

    rm(h5file)
end

@testset "Thermal Noise" begin
    obsconfig = readobsconfig("data/config1.yaml")
    stationinfo = readstationinfo(obsconfig["stations"], delim=",", ignorerepeated=false)
    h5file = "noise.h5"

    ms = readms(obsconfig["msname"])

    @inferred thermalnoise!(ms, stationinfo, obsconfig, h5file=h5file)
    rm(h5file)

    @inferred thermalnoise!(ms, stationinfo, obsconfig, h5file=h5file, noisefile="data/insmodel1.h5")
    rm(h5file)

    @inferred thermalnoise!(ms, stationinfo, obsconfig, h5file=h5file, noisefile="data/insmodel2.h5")
    rm(h5file)

end