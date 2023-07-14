export plotcoherencyvis, plotdata

"""
    plotmodel(obs::CjlObservation)

Plot model visibilities before applying instrumental model
"""
function plotcoherencyvis(obs::CjlObservation)
    uvwave = sqrt.(obs.uvw[1,:].^2 .+ obs.uvw[2,:].^2) / (299792458.0/mean(obs.chanfreqvec)) / 1e9 # in units of Gλ
    #p = plot(uvwave, abs.(getindex.(obs.data[:],1)), seriestype=:scatter, ls=:dot, ms=1, mc=:red, markerstrokecolor=:red, label="RR")
    #plot(p, uvwave, abs.(getindex.(obs.data[:],2)), seriestype=:scatter, ls=:dot, ms=1, mc=:cyan, markerstrokecolor=:cyan, label="RL")
    #plot(p, uvwave, abs.(getindex.(obs.data[:],3)), seriestype=:scatter, ls=:dot, ms=1, mc=:purple, markerstrokecolor=:purple, label="LR")
    #plot(p, uvwave, abs.(getindex.(obs.data[:],4)), seriestype=:scatter, ls=:dot, ms=1, mc=:green, markerstrokecolor=:green, label="LL")

    p = nothing
    for freqid in 1:obs.numchan
        println(freqid)
        if freqid == 1
            p = plot(uvwave, abs.(obs.data[1,1,freqid,:]), seriestype=:scatter, ls=:dot, ms=1, mc=:red, markerstrokecolor=:red, label="RR")
        else
            plot(p, uvwave, abs.(obs.data[1,1,freqid,:]), seriestype=:scatter, ls=:dot, ms=1, mc=:red, markerstrokecolor=:red, label="RR")
        end
    end
    savefig(p, "01_coherencyvis_amplitude_vs_uvwave.png")
end

"""
    plotdata()

"""
function plotdata()
end