@testset "MS from scratch" begin
    msname = "eht.ms"
    mode = "manual"
    stations = "data/eht_2017.stations"
    casaanttemplate = "data/antenna_table.template"
    spw_centrefreq = [229.0e9]
    spw_bw = [2.0e9]
    spw_channels = [32]
    sourcedict = Dict{String, Any}("M87" => Dict{String, String}("RA"=>"12h30m49.42", "Dec"=>"+12.23.28.04", "epoch"=>"J2000"))
    starttime = "UTC,2021/04/28/00:00:00.00"
    exposure = 1.0
    scans = 2
    scanlengths = [900.0, 600.0]
    scanlag = 300.0
    autocorr = false
    telescopename = "VLBA"
    feed = "perfect R L"
    shadowlimit = 1e-6
    elevationlimit = "10deg"
    stokes = "RR RL LR LL"
    delim = ","
    ignorerepeated = false

    @info(pwd())

    @test "somevalue" != Anime.msfromconfig(msname, mode, stations, casaanttemplate, spw_centrefreq, spw_bw, spw_channels, sourcedict, starttime,
    exposure, scans, scanlengths, scanlag, autocorr=autocorr, telescopename=telescopename, feed=feed, shadowlimit=shadowlimit, elevationlimit=elevationlimit,
    stokes=stokes, delim=delim, ignorerepeated=ignorerepeated)

    rm(msname, force=true, recursive=true)
end

@testset "MS from UVFITS" begin
    uvfits = "data/hops_lo_3601_M87+zbl-dtcal_selfcal.uvfits"
    msname = "eht.ms"
    stations = "data/eht_2017.stations"
    mode = "uvfits"
    delim = ","
    ignorerepeated = false

    @test "somevalue" != Anime.msfromuvfits(uvfits, msname, mode, stations, delim=delim, ignorerepeated=ignorerepeated)

    rm(msname, force=true, recursive=true)
end