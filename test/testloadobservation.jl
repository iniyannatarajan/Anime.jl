@testset "Load MS" begin
    y = Dict("x" => "5.0")
    msname = "data/eht.ms"
    stations = "data/eht_2017.stations"
    corruptseed = 4534
    tropseed = 83746

    @inferred CjlObservation loadms(y, msname, stations, corruptseed, tropseed, delim=",", ignorerepeated=false)
end