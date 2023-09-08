export genseries1d!

"""
    genseries1d!(series::Vector{ComplexF32}, mode::String, location::ComplexF32, scale::Float32, driftrate::Float32, nsamples::Int64, rng::AbstractRNG)

Generate a complex-valued Gaussian process 1-D series of length nsamples with the given `location`, `scale`, and `driftrate` parameters. `mode` determines
if "normal" distribution or "gaussian processes" is used to generate the samples.
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

Generate a Float32-valued Gaussian process 1-D series of length nsamples with the given `location`, `scale`, and `driftrate` parameters. `mode` determines
if "normal" distribution or "gaussian processes" is used to generate the samples.
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

Generate a Float64-valued Gaussian process 1-D series of length nsamples with the given `location`, `scale`, and `driftrate` parameters. `mode` determines
if "normal" distribution or "gaussian processes" is used to generate the samples.
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
    squaredexponentialkernel(x1, x2; σ=1.0, ρ=1.0)

Generate squared exponential kernel function of the form
```math
k_{SE}(x-x') = \\sigma^2 e^{-\\frac{(x-x')^2}{2\\rho^2}}
```
where σ^2 is the variance and ρ is the characteristic length.
"""
function squaredexponentialkernel(x1, x2; σ=1.0, ρ=1.0)
    return σ^2 * exp(-0.5 * ((x1 - x2)^2 / ρ^2))
end

"""
    genseries1d!(series, times, rng::AbstractRNG; μ=0.0, σ=1.0, ρ=1.0)

Generate 1-D series using SE kernel.
"""
function genseries1d!(series, times, rng::AbstractRNG; μ=0.0, σ=1.0, ρ=1.0)
    # Compute covariance matrix
    ntimes = length(times)
    covmat= zeros(ntimes, ntimes)
    for i in 1:ntimes
        for j in 1:ntimes
            covmat[i, j] = squaredexponentialkernel(times[i], times[j], σ=σ, ρ=ρ)
        end
    end

    covmat += 1e-6 * I # add small constant to make the covariance matrix positive definite

    # Generate random samples
    meanvector = zeros(ntimes) .+ μ # offset mean to the value of μ

    ndims = length(meanvector)
    stdnormalsample = randn(rng, Float64, ndims)

    # Transform to multivariate normal sample
    series[:] = meanvector .+ cholesky(covmat).L * stdnormalsample
    return series
end

"""
    rationalquadratickernel(x1, x2; σ=1.0, α=1.0, ρ=2.0)

Generate rational quadratic kernel of the form
```math
k_{RQ}(x-x') = \\sigma^2 \\big(1+\\frac{(x-x')^2}{2\\alpha\\rho^2}\\big)^{-\\alpha}
```
"""
function rationalquadratickernel(x1, x2; σ=1.0, α=1.0, ρ=2.0)
    return σ^2 * (1 + ((x1 - x2)^2 / 2*α*ρ^2))^-α
end

"""
    genseries1d!(series, times, rng::AbstractRNG, α::Float64; μ=0.0, σ=1.0, ρ=2.0)

Generate 1-D series using RQ kernel.
"""
function genseries1d!(series, times, rng::AbstractRNG, α::Float64; μ=0.0, σ=1.0, ρ=2.0)
    # Compute covariance matrix
    ntimes = length(times)
    covmat= zeros(ntimes, ntimes)
    for i in 1:ntimes
        for j in 1:ntimes
            covmat[i, j] = rationalquadratickernel(times[i], times[j], σ=σ, α=α, ρ=ρ)
        end
    end

    covmat += 1e-6 * I # add small constant to make the covariance matrix positive definite

    # Generate random samples
    meanvector = zeros(ntimes) .+ μ # offset mean to the value of μ

    ndims = length(meanvector)
    stdnormalsample = randn(rng, Float64, ndims)

    # Transform to multivariate normal sample
    series[:] = meanvector .+ cholesky(covmat).L * stdnormalsample
    return series
end

#="""
    euclideandistance(x1, x2)

Compute Euclidean distance between two vectors.
"""
function euclideandistance(x1, x2)
    return sqrt(sum((x1 .- x2).^2))
end

"""
    maternkernel(r, ν, ρ)

Compute Matern Kernel.
"""
function maternkernel(r, ν, ρ)
    return (2^(1 - ν) / gamma(ν)) * (sqrt(2 * ν) * r / ρ)^ν * exp(-sqrt(2 * ν) * r / ρ)
end

"""
    gen1dseries!(series, samplepoints; ν=2.0, ρ=1.0)

Generate a 1-D series using Matern kernel.
"""
function genseries1d!(series, samplepoints; ν=2.0, ρ=1.0)
    npoints = size(samplepoints, 1)
    distancematrix = zeros(npoints, npoints)
    
    for i in 1:npoints
        for j in 1:npoints
            distancematrix[i, j] = euclideandistance(samplepoints[i, :], samplepoints[j, :])
        end
    end
    
    covmat = maternkernel.(distancematrix, ν, ρ)
    covmat += 1e-6 * I

    series[:] = cholesky(covmat).L
    return series
end=#
