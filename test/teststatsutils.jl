@testset "Time Series" begin
    @inferred gentimeseries!(zeros(ComplexF32, 100), "gp", ComplexF32(0.0+0.0*im), Float32(1.0), Float32(3.0), 100, Xoshiro(42))
    @inferred gentimeseries!(zeros(ComplexF32, 100), "normal", ComplexF32(0.0+0.0*im), Float32(1.0), Float32(3.0), 100, Xoshiro(42))

    @inferred gentimeseries!(zeros(Float32, 100), "gp", Float32(0.0), Float32(1.0), Float32(3.0), 100, Xoshiro(42))
    @inferred gentimeseries!(zeros(Float32, 100), "normal", Float32(0.0), Float32(1.0), Float32(3.0), 100, Xoshiro(42))

    @inferred gentimeseries!(zeros(Float64, 100), "gp", 0.0, 1.0, 3.0, 100, Xoshiro(42))
    @inferred gentimeseries!(zeros(Float64, 100), "normal", 0.0, 1.0, 3.0, 100, Xoshiro(42))

    series = gentimeseries!(zeros(Float64, 100), collect(range(start=1.0, stop=10.0, length=100)), Xoshiro(42))
    @test size(series) == (100,)
end