@testset "Read Station Info" begin
    @inferred readstationinfo("data/eht_2017.stations")
end

@testset "Read Bandpass Info" begin
    @inferred readbandpassinfo("data/eht_2017.bandpass")
end

@testset "Read Alist v6" begin
    @inferred readalistv5("data/alist.v6")
end

@testset "Read MS" begin
    @inferred readms("data/eht1.ms")
end

@testset "Create MS from config" begin
    obsconfig = Dict()
    obsconfig["msname"] = "eht.ms"
    obsconfig["mode"] = "manual"
    obsconfig["stations"] = "data/eht_2017.stations"
    obsconfig["casaanttemplate"] = "data/antenna_table.template"
    obsconfig["spw"] = Dict()
    obsconfig["spw"]["centrefreq"] = [229.0e9]
    obsconfig["spw"]["bandwidth"] = [2.0e9]
    obsconfig["spw"]["channels"] = [32]
    obsconfig["source"] = Dict{String, Any}("M87" => Dict{String, String}("RA"=>"12h30m49.42", "Dec"=>"+12.23.28.04", "epoch"=>"J2000"))
    obsconfig["starttime"] = "UTC,2021/04/28/00:00:00.00"
    obsconfig["exposure"] = 1.0
    obsconfig["scans"] = 2
    obsconfig["scanlengths"] = [900.0, 600.0]
    obsconfig["scanlag"] = 300.0
    obsconfig["autocorr"] = false
    obsconfig["telescopename"] = "VLBA"
    obsconfig["feed"] = "perfect R L"
    obsconfig["shadowlimit"] = 1e-6
    obsconfig["elevationlimit"] = "10deg"
    obsconfig["stokes"] = "RR RL LR LL"
    delim = ","
    ignorerepeated = false

    stationinfo = readstationinfo(obsconfig["stations"])

    @test "somevalue" != Anime.createmsfromconfig(obsconfig, stationinfo)

    rm(obsconfig["msname"], force=true, recursive=true)
    rm("ANTENNA_eht_2017", force=true, recursive=true)
end

@testset "Create MS from UVFITS" begin
    uvfits = "data/hops_lo_3601_M87+zbl-dtcal_selfcal.uvfits"
    msname = "eht.ms"
    mode = "uvfits"

    @test "somevalue" != Anime.createmsfromuvfits(uvfits, msname, mode)

    rm(msname, force=true, recursive=true)
end

@testset "Write MS" begin
    ms = readms("data/eht1.ms")
    h5file = "data/insmodel1.h5"

    @inferred writems(ms, h5file=h5file)
end