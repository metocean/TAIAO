import xarray as xr
import pandas as pd
import sys

def read_file(filepath,sensor=None):
    try:
        df = xr.open_dataset(filepath)
    except:
        print("Failed to read ", filepath)
        sys.exit()

    if sensor:
        df=df.sel({'sensor':sensor})

    return df.sel({'site':0}).to_dataframe()

def output_tonc(qcname,df): 


    # Write results to netcdf file
    dset = xr.Dataset.from_dataframe(df)
    dset.time.attrs.update({'standard_name': 'time'})
    dset.elev.attrs.update({'standard_name': 'sea_surface_height_above_reference_ellipsoid',
                            'units': 'm',
                            '_FillValue': -999.})
    dset.ss.attrs.update({'standard_name': 'non_tidal_elevation_of_sea_surface_height',
                            'units': 'm',
                            '_FillValue': -999.})
    dset.tide.attrs.update({'standard_name': 'tidal_elevation_of_sea_surface_height',
                            'units': 'm',
                            '_FillValue': -999.})
    dset.to_netcdf(qcname,
                   encoding={"elev": {"dtype": "short",
                                      "scale_factor": 0.001,
                                      "zlib": True, "complevel": 6},
                             "ss": {"dtype": "short",
                                      "scale_factor": 0.001,
                                      "zlib": True, "complevel": 6},
                             "tide": {"dtype": "short",
                                      "scale_factor": 0.001,
                                      "zlib": True, "complevel": 6},
                             "time": {"zlib": True, "complevel": 6}})

    del df