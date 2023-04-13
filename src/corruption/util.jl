export parallacticangle, elevationangle, gentimeseries!

quanta = pyimport("casatools" => "quanta")
qa = quanta()

me = measures()

"""
    parallacticangle(obs::CjlObservation)

Compute parallactic angle for all stations for all times.
"""
function parallacticangle(obs::CjlObservation)
    # get unique times
    uniqtimes = unique(obs.times)

    ra = qa.quantity(obs.phasedir[1], "rad")
    dec = qa.quantity(obs.phasedir[2], "rad")

    pointing = me.direction("j2000", ra, dec)
    starttime = me.epoch("utc", qa.quantity(uniqtimes[1], "s"))
    me.doframe(starttime)

    nant = size(obs.stationinfo)[1]

    parallacticanglematrix = zeros(size(uniqtimes)[1], nant)

    for ant in 1:nant
        x = qa.quantity(obs.pos[1, ant], "m")
        y = qa.quantity(obs.pos[2, ant], "m")
        z = qa.quantity(obs.pos[3, ant], "m")
        position = me.position("wgs84", x, y, z)
        me.doframe(position)
        sec2rad = 2*pi/(24.0*3600.0)
        hourangle = pyconvert(Float64, me.measure(pointing, "HADEC")["m0"]["value"]) .+ (uniqtimes.-minimum(uniqtimes)).*sec2rad
        earthradius = 6371000.0
        latitude = asin(obs.pos[3, ant]/earthradius)
	parallacticanglematrix[:,ant] = atan.(sin.(hourangle).*cos(latitude), (cos(obs.phasedir[2])*sin(latitude).-cos.(hourangle).*cos(latitude).*sin(obs.phasedir[2])))
    end

    return parallacticanglematrix
end

"""
    elevationangle(obs::CjlObservation)

Compute elevation angle for all stations for all times.

"""
function elevationangle(obs::CjlObservation)
    # get unique times
    uniqtimes = unique(obs.times)

    ra = qa.quantity(obs.phasedir[1], "rad")
    dec = qa.quantity(obs.phasedir[2], "rad")
    
    pointing = me.direction("j2000", ra, dec)
    starttime = me.epoch("utc", qa.quantity(uniqtimes[1], "s"))
    me.doframe(starttime)

    nant = size(obs.stationinfo)[1]
    
    elevationmatrix = zeros(size(uniqtimes)[1], nant)
    
    for ant in 1:nant
        x = qa.quantity(obs.pos[1, ant], "m")
        y = qa.quantity(obs.pos[2, ant], "m")
        z = qa.quantity(obs.pos[3, ant], "m")
        position = me.position("wgs84", x, y, z)
        me.doframe(position)
        sec2rad = 2*pi/(24.0*3600.0)
        hourangle = pyconvert(Float64, me.measure(pointing, "HADEC")["m0"]["value"]) .+ (uniqtimes.-minimum(uniqtimes)).*sec2rad
        earthradius = 6371000.0
        latitude = asin(obs.pos[3, ant]/earthradius)
        elevationmatrix[:,ant] = asin.(sin(latitude)*sin(obs.phasedir[2]).+cos(latitude)*cos(obs.phasedir[2]).*cos.(hourangle))
    end

    elevationmatrix[elevationmatrix .< 0] .= 0.0 # set all -ve values to 0.0 instead of NaN # set all negative values to NaN
    return elevationmatrix
end

"""
    gentimeseries!(series::Vector{ComplexF32}, mode::String, location::ComplexF32, scale::Float64, driftrate::Float64, nsamples::Int64, rng::AbstractRNG)

Generate a complex-valued Gaussian process time-series of length nsamples with the given location, scale, and driftrate parameters.
"""
function gentimeseries!(series::Vector{ComplexF32}, mode::String, location::ComplexF32, scale::Float64, driftrate::Float64, nsamples::Int64, rng::AbstractRNG)
    # TODO this is a crude version of a wiener process -- to be updated
    if mode == "gp"
        sqrtnsamples = sqrt(nsamples)
        series[1] = location + scale*randn(rng, ComplexF32)
        for ii in 2:nsamples
            series[ii] = series[ii-1] + (scale*randn(rng, ComplexF32)/sqrtnsamples) + driftrate*ii
        end
    elseif mode == "normal"
        series = location .+ scale*randn(rng, ComplexF32, nsamples)
    end
    return series # return by convention
end

"""
    gentimeseries!(series::Vector{Float32}, mode::String, location::Float32, scale::Float64, driftrate::Float64, nsamples::Int64, rng::AbstractRNG)

Generate a complex-valued Gaussian process time-series of length nsamples with the given location, scale, and driftrate parameters.
"""
function gentimeseries!(series::Vector{Float32}, mode::String, location::Float32, scale::Float64, driftrate::Float64, nsamples::Int64, rng::AbstractRNG)
    # TODO this is a crude version of a wiener process -- to be updated
    if mode == "gp"
        sqrtnsamples = sqrt(nsamples)
        series[1] = location + scale*randn(rng, Float32)
        for ii in 2:nsamples
            series[ii] = series[ii-1] + (scale*randn(rng, Float32)/sqrtnsamples) + driftrate*ii
        end
    elseif mode == "normal"
        series = location .+ scale*randn(rng, Float32, nsamples)
    end
    return series
end

"""
    gentimeseries!(series::Vector{Float64}, mThis function generates a complex-valued Gaussian process time-series of length nsamples with the given location, scale, and driftrate parameters.ode::String, location::Float64, scale::Float64, driftrate::Float64, nsamples::Int64, rng::AbstractRNG)

Generate a complex-valued Gaussian process time-series of length nsamples with the given location, scale, and driftrate parameters.
"""
function gentimeseries!(series::Vector{Float64}, mode::String, location::Float64, scale::Float64, driftrate::Float64, nsamples::Int64, rng::AbstractRNG)
    # TODO this is a crude version of a wiener process -- to be updated
    if mode == "gp"
        sqrtnsamples = sqrt(nsamples)
        series[1] = location + scale*randn(rng, Float32)
        for ii in 2:nsamples
            series[ii] = series[ii-1] + (scale*randn(rng, Float32)/sqrtnsamples) + driftrate*ii
        end
    elseif mode == "normal"
        series = location .+ scale*randn(rng, Float32, nsamples)
    end
    return series
end
