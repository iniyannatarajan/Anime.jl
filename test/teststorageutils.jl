@testset "Post-processing" begin
    y = YAML.load_file("data/config1.yaml", dicttype=Dict{String,Any}) # sample dict to test loadms()
    h5file = "data/insmodel1.h5"

    obs = loadms(y["msname"], y["stations"], Int(y["corruptseed"]), Int(y["troposphere"]["tropseed"]), y["troposphere"]["wetonly"], y["correff"], 
    y["troposphere"]["attenuate"], y["troposphere"]["skynoise"], y["troposphere"]["meandelays"], y["troposphere"]["turbulence"], 
    y["instrumentalpolarization"]["visibilityframe"], y["instrumentalpolarization"]["mode"], y["pointing"]["interval"], y["pointing"]["scale"], y["pointing"]["mode"], y["stationgains"]["mode"], 
    y["bandpass"]["bandpassfile"], delim=",", ignorerepeated=false)

    postprocessms(obs, h5file=h5file)
end