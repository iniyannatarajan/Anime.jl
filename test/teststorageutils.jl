@testset "Post-processing" begin
    y = YAML.load_file("data/testconfig.yaml", dicttype=Dict{String,Any}) # sample dict to test loadms()
    h5file = "data/insmodel.h5"

    obs = loadms(y["msname"], y["stations"], Int(y["corruptseed"]), Int(y["troposphere"]["tropseed"]), y["troposphere"]["wetonly"], y["correff"], 
    y["troposphere"]["attenuate"], y["troposphere"]["skynoise"], y["troposphere"]["meandelays"], y["troposphere"]["turbulence"], 
    y["instrumentalpol"]["visibilityframe"], y["instrumentalpol"]["mode"], y["pointing"]["interval"], y["pointing"]["mode"], y["stationgains"]["mode"], 
    y["bandpass"]["bandpassfile"], delim=",", ignorerepeated=false)

    postprocessms(obs, h5file=h5file)
end