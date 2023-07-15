export plotvisamp_vs_pbs

"""
    plotvisamp_vs_pbs(data::Array{Complex{Float32},4}, flag::Array{Bool,4}, uvw::Matrix{Float64}, chanfreqvec::Array{Float64,1}, numchan::Int64)

Plot visibility amplitudes against projected baseline separation
"""
function plotvisamp_vs_pbs(data::Array{Complex{Float32},4}, flag::Array{Bool,4}, uvw::Matrix{Float64}, chanfreqvec::Array{Float64,1}, numchan::Int64; saveas="01_coherencyamplitude_vs_projectedbaseline.png")
    uvwave = sqrt.(uvw[1,:].^2 .+ uvw[2,:].^2) / (299792458.0/mean(chanfreqvec)) / 1e9 # in units of Gλ

    maskindices = findall(isequal(true), flag)
    maskeddata = deepcopy(data)
    maskeddata[maskindices] .= NaN # assign NaN values to flagged visibilities

    # plot first frequency channel
    p = plot(uvwave, abs.(maskeddata[1,1,1,:]), seriestype=:scatter, ls=:dot, ms=1, mc=:red, markerstrokecolor=:red, label="RR")
    plot!(p, uvwave, abs.(maskeddata[1,2,1,:]), seriestype=:scatter, ls=:dot, ms=1, mc=:cyan, markerstrokecolor=:cyan, label="RL")
    plot!(p, uvwave, abs.(maskeddata[2,1,1,:]), seriestype=:scatter, ls=:dot, ms=1, mc=:purple, markerstrokecolor=:purple, label="LR")
    plot!(p, uvwave, abs.(maskeddata[2,2,1,:]), seriestype=:scatter, ls=:dot, ms=1, mc=:green, markerstrokecolor=:green, label="LL")
    
    # plot the rest of the frequency channels
    if numchan > 1
        plot!(p, uvwave, abs.(maskeddata[1,1,2:end,:]'), seriestype=:scatter, ls=:dot, ms=1, mc=:red, markerstrokecolor=:red, label="")
        plot!(p, uvwave, abs.(maskeddata[1,2,2:end,:]'), seriestype=:scatter, ls=:dot, ms=1, mc=:cyan, markerstrokecolor=:cyan, label="")
        plot!(p, uvwave, abs.(maskeddata[2,1,2:end,:]'), seriestype=:scatter, ls=:dot, ms=1, mc=:purple, markerstrokecolor=:purple, label="")
        plot!(p, uvwave, abs.(maskeddata[2,2,2:end,:]'), seriestype=:scatter, ls=:dot, ms=1, mc=:green, markerstrokecolor=:green, label="")
    end

    plot!(p, xlabel="Projected baseline separation (Gλ)", ylabel="Source coherency amplitude (Jy)") #, legend=:outerbottom, legendcolumns=4)
    savefig(p, saveas)
end

"""
    plotdata()

"""
function plotdata()
end