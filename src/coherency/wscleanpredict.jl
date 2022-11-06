using PythonCall

table = pyimport("casatools" => "table")
tb = table()

function runwsclean(msname::String, modelimage::String, toggle_polmodel::Bool, channelgroups::UInt64, start_timestep::UInt64, end_timestep::UInt64, osfactor::UInt64)
    toggle_polmodel ? `wsclean -channels-out $channelgroups -predict -name $modelimage -interval $start_timestep $end_timestep -pol I,Q,U,V -no-reorder -oversampling $osfactor -no-small-inversion $msname` : `wsclean -channels-out $channelgroups -predict -name $modelimage -interval $start_timestep $end_timestep -oversampling $osfactor -no-small-inversion $msname`
end

function slicewscleaninputs(msname::String, fitsdir::String, toggle_polmodel::Bool, channelgroups::UInt64, )
    fitsfiles = readdir(fitsdir, sort=true)
    nfits = 0
    if toggle_polmodel
        length(fitsfiles)%4 ? nfits = length(fitsfiles)/channelgroups/4 : error("toggle_polmodel is set to true but some polarisation images are missing!")
    else
	nfitsfiles = length(fitsfiles)/channelgroups
    end

    # read in the MS
    tb.open(msname)
    times = tb.getcol("TIME")
    uniqtimes = unique(times) 

    # compute timesteps per fits image -- BUT NTIMES MUST BE READ FROM THE MS! use mutable struct somewhere and store values? ALSO, how to read model images for all scans?
    ntimes_per_fits = floor(Int, ntimes/nfitsfiles)

    start_timestep = 0
    endvis = ntimes_per_fits

            for img_ind in range(self.num_images):
                temp_input_fits = '%s/t%04d'%(self.input_fitsimage,img_ind)
                info('Simulating visibilities (corr dumps) from %d to %d using input sky model %s'%(startvis,endvis,temp_input_fits))
                run_wsclean(temp_input_fits, self.input_fitspol, self.input_changroups, startvis, endvis, self.oversampling)
                startvis = endvis
                if img_ind != self.num_images-2:
                    endvis = endvis + self.vis_per_image
                else:
                    endvis = endvis + 2*self.vis_per_image # INI: ensure all vis at the end are accounted for in the next (last) iteration.

            # INI: Copy over data from MODEL_DATA to output_column if output_column is not MODEL_DATA
            if self.output_column != 'MODEL_DATA':
                copy_between_cols(self.output_column, 'MODEL_DATA')
