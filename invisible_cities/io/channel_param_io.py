import tables as tb
from .. reco     import tbl_functions as tbl


## A dictionary to define the structure of the output table
## At initialization, the function specific entries will
## be defined.
calibration_parameters = {
    'SensorID' : tb.UInt32Col(pos=0)
    }


## Suggested parameter list
## Useful for generic/default fit functions
generic_params = ["normalization", "poisson_mu"    ,
                  "pedestal"     , "pedestal_sigma",
                  "gain"         , "gain_sigma"    ,
                  "fit_limits"   , "n_gaussians_chi2" ]


def create_param_table(h5out, sensor_type, func_name, param_names, covariance=None):
    ## First sort out the column names:
    for i, par in enumerate(param_names):
        calibration_parameters[par] = tb.Float32Col(pos=i + 1, shape=2)

    if covariance:
        # We're saving the covariance matrix too
        calibration_parameters["covariance"] = tb.Float32Col(pos=len(param_names)+1,
                                                             shape=covariance   )
    
    ## If the group 'FITPARAMS' doesn't already exist, create it
    try:                       PARAM_group = getattr(h5out.root, "FITPARAMS")
    except tb.NoSuchNodeError: PARAM_group = h5out.create_group(h5out.root,
                                                                "FITPARAMS")
    ## Define a table for this fitting function
    param_table = h5out.create_table(PARAM_group,
                                     "FIT_"+sensor_type+"_"+func_name,
                                     calibration_parameters,
                                     "Calibration parameters",
                                     tbl.filters("NOCOMPR"))
    return param_table


def store_fit_values(param_table, sensor_id, fit_result):

    channel = param_table.row
    channel["SensorID"] = sensor_id
    for key, param in fit_result.items():
        channel[key] = param
    channel.append()
    param_table.flush()


def channel_param_writer(h5out, *, sensor_type,
                         func_name, param_names,
                         covariance=None):
    """
    Define a group, table and writer function for the ouput
    of fit parameters.
    input:
    h5out : the predefined hdf5 output file
    sensor_type : e.g. pmt, sipm or FE
    func_name : A string with the name of the fitting function being used
    param_names : A list of parameter names
    covariance : None or a tuple with the shape of the covariance matrix
    """
    
    param_table = create_param_table(h5out, sensor_type,
                                     func_name, param_names, covariance)
    def store_channel_fit(sensor_id, fit_result):
        """
        input:
        channel_id : Sensor number
        fit_result : dict with keys as parameter names
                     Fit parameters should be (value, error)
        """
        store_fit_values(param_table, sensor_id, fit_result)
    
    return store_channel_fit
