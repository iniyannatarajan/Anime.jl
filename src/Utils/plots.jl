export plotcoherencyvis, plotdata

"""
    plotmodel(obs::CjlObservation)

Plot model visibilities before applying instrumental model
"""
function plotcoherencyvis(obs::CjlObservation)
    uvwave = sqrt.(obs.uvw[1,:].^2 .+ obs.uvw[2,:].^2) / (299792458.0/mean(obs.chanfreqvec)) / 1e9 # in units of Gλ

    # plot first frequency channel
    p = plot(uvwave, abs.(obs.data[1,1,1,:]), seriestype=:scatter, ls=:dot, ms=1, mc=:red, markerstrokecolor=:red, label="RR")
    plot!(p, uvwave, abs.(obs.data[1,2,1,:]), seriestype=:scatter, ls=:dot, ms=1, mc=:cyan, markerstrokecolor=:cyan, label="RL")
    plot!(p, uvwave, abs.(obs.data[2,1,1,:]), seriestype=:scatter, ls=:dot, ms=1, mc=:purple, markerstrokecolor=:purple, label="LR")
    plot!(p, uvwave, abs.(obs.data[2,2,1,:]), seriestype=:scatter, ls=:dot, ms=1, mc=:green, markerstrokecolor=:green, label="LL")
    
    # plot the rest of the frequency channels
    if obs.numchan > 1
        plot!(p, uvwave, abs.(obs.data[1,1,2:end,:]'), seriestype=:scatter, ls=:dot, ms=1, mc=:red, markerstrokecolor=:red, label="")
        plot!(p, uvwave, abs.(obs.data[1,2,2:end,:]'), seriestype=:scatter, ls=:dot, ms=1, mc=:cyan, markerstrokecolor=:cyan, label="")
        plot!(p, uvwave, abs.(obs.data[2,1,2:end,:]'), seriestype=:scatter, ls=:dot, ms=1, mc=:purple, markerstrokecolor=:purple, label="")
        plot!(p, uvwave, abs.(obs.data[2,2,2:end,:]'), seriestype=:scatter, ls=:dot, ms=1, mc=:green, markerstrokecolor=:green, label="")
    end

    plot!(p, xlabel="Projected baseline separation (Gλ)", ylabel="Source coherency amplitude (Jy)") #, legend=:outerbottom, legendcolumns=4)
    savefig(p, "01_coherencyamplitude_vs_projectedbaseline.png")
end

"""
    plotdata()

"""
function plotdata()
end