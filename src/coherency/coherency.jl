export run_wsclean

tb = table()

"""
    copymodeltodata(msname::String)

Copy ```MODEL_DATA``` to ```DATA``` in MS `msname`.
"""
function copymodeltodata(msname::String)
    # copying model data to data
    tb.open(msname, nomodify=false)
    tb.putcol("DATA", tb.getcol("MODEL_DATA"))
    tb.close()
    tb.clearlocks()
end

"""
    run_wsclean(msname::String, fitsdir::String, polarized::Bool, channelgroups::Int64, osfactor::Int64)

Compute source coherency matrix using `WSClean` and populate MS. `fitsdir` is the directory containing FITS source models,
`polarized` indicates if source model is polarized, and `channelgroups` is the number of images in frequency.
`osfactor` is the oversampling factor which is a `WSClean` argument that controls prediction/imaging accuracy.
"""
function run_wsclean(msname::String, fitsdir::String, polarized::Bool, channelgroups::Int64, osfactor::Int64)
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
    @info("Computing source coherency with WSClean...")
    for fitsindex in 0:nmodels-1
	    infits = "$(fitsdir)/t$(lpad(fitsindex, 4, "0"))"
	    polarized ? run(`wsclean -channels-out $channelgroups -predict -name $infits -interval $startrow $endrow -pol I,Q,U,V -reorder -oversampling $osfactor -no-small-inversion $msname`) : run(`wsclean -channels-out $channelgroups -predict -name $infits -interval $startrow $endrow -oversampling $osfactor -no-small-inversion $msname`)
        startrow = endrow
        fitsindex == nmodels-2 ? endrow += 2*rows_per_modelimg : endrow += rows_per_modelimg
    end

    # copy MODEL_DATA to DATA -- all corruptions will be added to DATA
    copymodeltodata(msname)

    @info("Compute source coherency with WSClean 🙆")

end

#="""
    run_ducc0(msname::String, fitsdir::String, polarized::Bool, channelgroups::Int64, osfactor::Int64)

Predict uncorrupted visibilities using the ducc0 wgridder (experimental)
"""
function run_ducc0(msname::String, fitsdir::String, polarized::Bool, channelgroups::Int64, osfactor::Int64)

end=#