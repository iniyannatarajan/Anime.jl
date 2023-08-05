export gentimeseries!

"""
    gentimeseries!(series::Vector{ComplexF32}, mode::String, location::ComplexF32, scale::Float32, driftrate::Float32, nsamples::Int64, rng::AbstractRNG)

Generate a complex-valued Gaussian process time-series of length nsamples with the given location, scale, and driftrate parameters.
"""
function gentimeseries!(series::Vector{ComplexF32}, mode::String, location::ComplexF32, scale::Float32, driftrate::Float32, nsamples::Int64, rng::AbstractRNG)
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
    gentimeseries!(series::Vector{Float32}, mode::String, location::Float32, scale::Float32, driftrate::Float32, nsamples::Int64, rng::AbstractRNG)

Generate a Float32-valued Gaussian process time-series of length nsamples with the given location, scale, and driftrate parameters.
"""
function gentimeseries!(series::Vector{Float32}, mode::String, location::Float32, scale::Float32, driftrate::Float32, nsamples::Int64, rng::AbstractRNG)
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
    gentimeseries!(series::Vector{Float64}, mode::String, location::Float64, scale::Float64, driftrate::Float64, nsamples::Int64, rng::AbstractRNG)

Generate a Float64-valued Gaussian process time-series of length nsamples with the given location, scale, and driftrate parameters.
"""
function gentimeseries!(series::Vector{Float64}, mode::String, location::Float64, scale::Float64, driftrate::Float64, nsamples::Int64, rng::AbstractRNG)
    # TODO this is a crude version of a wiener process -- to be updated
    # TODO Look up a squared exponential kernel
    if mode == "gp"
        sqrtnsamples = sqrt(nsamples)
        series[1] = location + scale*randn(rng, Float64)
        for ii in 2:nsamples
            series[ii] = series[ii-1] + (scale*randn(rng, Float64)/sqrtnsamples) + driftrate*ii
        end
    elseif mode == "normal"
        series = location .+ scale*randn(rng, Float64, nsamples)
    end
    return series
end
