export runwsclean

tb = table()

function copymodeltodata(msname::String)
    # copying model data to data
    tb.open(msname, nomodify=false)
    tb.putcol("DATA", tb.getcol("MODEL_DATA"))
    tb.close()
    tb.clearlocks()
end

function runwsclean(msname::String, fitsdir::String, polarized::Bool, channelgroups::Int64, osfactor::Int64)
    """
    Run wsclean
    """
    fitsfiles = readdir(fitsdir, sort=true)
    nmodels = 0
    if polarized
        length(fitsfiles)%4 == 0 ? nmodels = floor(Int64, length(fitsfiles)/channelgroups/4) : error("polarized=true but polarized models missing!")
    else
        nmodels = floor(Int64, length(fitsfiles)/channelgroups)
    end

    # read in the MS
    tb.open(msname)
    times = tb.getcol("TIME")
    uniqtimes = unique(times)
    tb.close()
    tb.clearlocks()

    # compute timesteps per fits image
    rows_per_modelimg = floor(Int64, length(uniqtimes)/nmodels)

    # initialise start and end rows
    startrow = 0
    endrow = rows_per_modelimg

    # loop through FITS models (polarised or unpolarised)
    @info("Inverting source models to visibilities...")
    for fitsindex in 0:nmodels-1
	infits = "$(fitsdir)/t$(lpad(fitsindex, 4, "0"))"
	polarized ? run(`wsclean -channels-out $channelgroups -predict -name $infits -interval $startrow $endrow -pol I,Q,U,V -no-reorder -oversampling $osfactor -no-small-inversion $msname`) : run(`wsclean -channels-out $channelgroups -predict -name $infits -interval $startrow $endrow -oversampling $osfactor -no-small-inversion $msname`)
        startrow = endrow
        fitsindex == nmodels-2 ? endrow += 2*rows_per_modelimg : endrow += rows_per_modelimg
    end

    # copy MODEL_DATA to DATA -- all corruptions will be added to DATA
    copymodeltodata(msname)

    @info("Uncorrupted visibilities computed 🙆")

end
