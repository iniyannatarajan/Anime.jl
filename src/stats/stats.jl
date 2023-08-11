export genseries1d!

"""
    genseries1d!(series::Vector{ComplexF32}, mode::String, location::ComplexF32, scale::Float32, driftrate::Float32, nsamples::Int64, rng::AbstractRNG)

Generate a complex-valued Gaussian process 1-D series of length nsamples with the given location, scale, and driftrate parameters.
"""
function genseries1d!(series::Vector{ComplexF32}, mode::String, location::ComplexF32, scale::Float32, driftrate::Float32, nsamples::Int64, rng::AbstractRNG)
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
    genseries1d!(series::Vector{Float32}, mode::String, location::Float32, scale::Float32, driftrate::Float32, nsamples::Int64, rng::AbstractRNG)

Generate a Float32-valued Gaussian process 1-D series of length nsamples with the given location, scale, and driftrate parameters.
"""
function genseries1d!(series::Vector{Float32}, mode::String, location::Float32, scale::Float32, driftrate::Float32, nsamples::Int64, rng::AbstractRNG)
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
    genseries1d!(series::Vector{Float64}, mode::String, location::Float64, scale::Float64, driftrate::Float64, nsamples::Int64, rng::AbstractRNG)

Generate a Float64-valued Gaussian process 1-D series of length nsamples with the given location, scale, and driftrate parameters.
"""
function genseries1d!(series::Vector{Float64}, mode::String, location::Float64, scale::Float64, driftrate::Float64, nsamples::Int64, rng::AbstractRNG)
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

"""
    squaredexponentialkernel(x1, x2; σ=1.0, ℓ=1.0)

Generate squared exponential kernel function.
"""
function squaredexponentialkernel(x1, x2; σ=1.0, ℓ=1.0)
    return σ^2 * exp(-0.5 * ((x1 - x2)^2 / ℓ^2))
end

"""
    genseries1d!(series, times, rng::AbstractRNG; μ=0.0, σ=1.0, ℓ=1.0)

Generate a 1-D series using SE kernel.
"""
function genseries1d!(series, times, rng::AbstractRNG; μ=0.0, σ=1.0, ℓ=1.0)
    # Compute covariance matrix
    ntimes = length(times)
    covmat= zeros(ntimes, ntimes)
    for i in 1:ntimes
        for j in 1:ntimes
            covmat[i, j] = squaredexponentialkernel(times[i], times[j], σ=σ, ℓ=ℓ)
        end
    end

    # Add small constant to the diagonal for positive definiteness
    #covmat += 1e-6 * I
    covmat += 1e-6 * I

    # Generate random samples
    meanvector = zeros(ntimes) .+ μ # offset mean to the value of μ

    ndims = length(meanvector)
    stdnormalsample = randn(rng, Float64, ndims)

    # Transform to multivariate normal sample
    series[:] = meanvector .+ cholesky(covmat).L * stdnormalsample
    return series
end
