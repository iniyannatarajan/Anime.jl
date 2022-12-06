export instrumentalpol

include(joinpath("util.jl"))

function instrumentalpol(obs::CjlObservation)
    """
    compute and apply polarization leakage
    """

    # compute D-terms
    if obs.yamlconf["instrumentalpol"]["visibilityframe"] != "antenna"
	# compute in sky frame
	if obs.yamlconf["instrumentalpol"]["visibilityframe"] != "sky"
	    @warn("Visibility frame not recognised! Computing visibilities in sky frame...")
	end
        
	# dfkdfhkdhf
	
    else
        # compute in antenna frame
    end 
end
