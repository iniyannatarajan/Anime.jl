@testset "Troposphere" begin
    y = YAML.load_file("data/config1.yaml", dicttype=Dict{String,Any}) # sample dict to test loadms()
    h5file = "tropos.h5"

    obs = loadms(y["msname"], y["stations"], Int(y["corruptseed"]), Int(y["troposphere"]["tropseed"]), y["troposphere"]["wetonly"], y["correff"], 
    y["troposphere"]["attenuate"], y["troposphere"]["skynoise"], y["troposphere"]["meandelays"], y["troposphere"]["turbulence"], 
    y["instrumentalpolarization"]["visibilityframe"], y["instrumentalpolarization"]["mode"], y["pointing"]["interval"], y["pointing"]["mode"], y["stationgains"]["mode"], 
    y["bandpass"]["bandpassfile"], delim=",", ignorerepeated=false)

    @inferred troposphere!(obs, h5file, absorptionfile="data/absorption1.csv", dispersivefile="data/dispersive1.csv", elevfile="data/insmodel1.h5")

    rm("atm.csv")
    rm(h5file)
end

@testset "Polarization" begin
    y = YAML.load_file("data/config1.yaml", dicttype=Dict{String,Any}) # sample dict to test loadms()
    h5file = "inspol.h5"

    obs = loadms(y["msname"], y["stations"], Int(y["corruptseed"]), Int(y["troposphere"]["tropseed"]), y["troposphere"]["wetonly"], y["correff"], 
    y["troposphere"]["attenuate"], y["troposphere"]["skynoise"], y["troposphere"]["meandelays"], y["troposphere"]["turbulence"], 
    y["instrumentalpolarization"]["visibilityframe"], y["instrumentalpolarization"]["mode"], y["pointing"]["interval"], y["pointing"]["mode"], y["stationgains"]["mode"], 
    y["bandpass"]["bandpassfile"], delim=",", ignorerepeated=false)

    @inferred instrumentalpolarization!(obs, h5file=h5file, elevfile="data/insmodel1.h5", parangfile="data/insmodel1.h5")

    rm(h5file)
end

@testset "Primary Beam" begin
    y = YAML.load_file("data/config1.yaml", dicttype=Dict{String,Any}) # sample dict to test loadms()
    h5file = "beam.h5"

    obs = loadms(y["msname"], y["stations"], Int(y["corruptseed"]), Int(y["troposphere"]["tropseed"]), y["troposphere"]["wetonly"], y["correff"], 
    y["troposphere"]["attenuate"], y["troposphere"]["skynoise"], y["troposphere"]["meandelays"], y["troposphere"]["turbulence"], 
    y["instrumentalpolarization"]["visibilityframe"], y["instrumentalpolarization"]["mode"], y["pointing"]["interval"], y["pointing"]["mode"], y["stationgains"]["mode"], 
    y["bandpass"]["bandpassfile"], delim=",", ignorerepeated=false)

    @inferred pointing!(obs, h5file=h5file)

    rm(h5file)
end

@testset "Station Gains" begin
    y = YAML.load_file("data/config1.yaml", dicttype=Dict{String,Any}) # sample dict to test loadms()
    h5file = "gains.h5"

    obs = loadms(y["msname"], y["stations"], Int(y["corruptseed"]), Int(y["troposphere"]["tropseed"]), y["troposphere"]["wetonly"], y["correff"], 
    y["troposphere"]["attenuate"], y["troposphere"]["skynoise"], y["troposphere"]["meandelays"], y["troposphere"]["turbulence"], 
    y["instrumentalpolarization"]["visibilityframe"], y["instrumentalpolarization"]["mode"], y["pointing"]["interval"], y["pointing"]["mode"], y["stationgains"]["mode"], 
    y["bandpass"]["bandpassfile"], delim=",", ignorerepeated=false)

    @inferred stationgains(obs, h5file=h5file)

    rm(h5file)
end

@testset "Bandpass" begin
    y = YAML.load_file("data/config1.yaml", dicttype=Dict{String,Any}) # sample dict to test loadms()
    h5file = "bandpass.h5"

    obs = loadms(y["msname"], y["stations"], Int(y["corruptseed"]), Int(y["troposphere"]["tropseed"]), y["troposphere"]["wetonly"], y["correff"], 
    y["troposphere"]["attenuate"], y["troposphere"]["skynoise"], y["troposphere"]["meandelays"], y["troposphere"]["turbulence"], 
    y["instrumentalpolarization"]["visibilityframe"], y["instrumentalpolarization"]["mode"], y["pointing"]["interval"], y["pointing"]["mode"], y["stationgains"]["mode"], 
    y["bandpass"]["bandpassfile"], delim=",", ignorerepeated=false)

    @inferred bandpass(obs, h5file=h5file)

    rm(h5file)
end

@testset "Thermal Noise" begin
    y = YAML.load_file("data/config1.yaml", dicttype=Dict{String,Any}) # sample dict to test loadms()
    h5file = "noise.h5"

    obs = loadms(y["msname"], y["stations"], Int(y["corruptseed"]), Int(y["troposphere"]["tropseed"]), y["troposphere"]["wetonly"], y["correff"], 
    y["troposphere"]["attenuate"], y["troposphere"]["skynoise"], y["troposphere"]["meandelays"], y["troposphere"]["turbulence"], 
    y["instrumentalpolarization"]["visibilityframe"], y["instrumentalpolarization"]["mode"], y["pointing"]["interval"], y["pointing"]["mode"], y["stationgains"]["mode"], 
    y["bandpass"]["bandpassfile"], delim=",", ignorerepeated=false)

    @inferred thermalnoise(obs, h5file=h5file)
    rm(h5file)

    @inferred thermalnoise(obs, h5file=h5file, noisefile="data/insmodel1.h5")
    rm(h5file)

    @inferred thermalnoise(obs, h5file=h5file, noisefile="data/insmodel2.h5")
    rm(h5file)

end