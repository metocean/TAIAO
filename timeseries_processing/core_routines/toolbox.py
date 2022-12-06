import pandas as pd
import numpy as np
from datetime import datetime, timedelta
from itertools import groupby
from matplotlib.dates import num2date,date2num
from toto.filters.despike_phasespace3d import despike_phasespace3d
from toto.filters.lanczos_filter import lanczos_filter
from toto.filters.detrend import detrend
from storm_surge.timeseries_processing.core_routines import check
from scipy.stats import linregress
import copy


def store_data_as_netcdf(results,
                         lon,
                         lat,
                         outname):
    """
    Stores the results of the storm surge analysis to a netcdf file.
    """
    
    dset_out = results.to_xarray().drop('et')

    encoding_params = {'zlib': True, 'complevel':6,
                       'dtype': 'short',
                       'scale_factor': 0.001,
                       '_FillValue':-999.}


    dset_out.attrs = {'longitude': lon,
                      'latitude': lat}

    dset_out.elev.attrs.update({'long_name': 'water_level'})
    dset_out.trend.attrs.update({'long_name': 'trend'})
    dset_out.tide.attrs.update({'long_name': 'tide'})
    dset_out.msea.attrs.update({'long_name': 'monthly_sea_level_variation'})
    dset_out.ss.attrs.update({'long_name': 'storm_surge'})
    dset_out.res.attrs.update({'long_name': 'residual'})
    dset_out.skew_surge_magnitude.attrs.update({'long_name': 'skew_surge_magnitude'})
    dset_out.skew_surge_lag.attrs.update({'long_name': 'skew_surge_lag'})
    dset_out.tidal_elevation_maximum_over_tidal_cycle.attrs.update({'long_name': 'tidal_elevation_maximum_over_tidal_cycle'})
    dset_out.tidal_elevation_maximum_time_over_tidal_cycle.attrs.update({'long_name': 'tidal_elevation_maximum_time_over_tidal_cycle'})

    
    dset_out.elev.attrs.update({'units': 'm'})
    dset_out.trend.attrs.update({'units': 'm'})
    dset_out.tide.attrs.update({'units': 'm'})
    dset_out.msea.attrs.update({'units': 'm'})
    dset_out.ss.attrs.update({'units': 'm'})
    dset_out.res.attrs.update({'units': 'm'})
    dset_out.skew_surge_magnitude.attrs.update({'units': 'm'})
    dset_out.skew_surge_lag.attrs.update({'units': 'seconds'})
    dset_out.tidal_elevation_maximum_over_tidal_cycle.attrs.update({'units': 'm'})


    dset_out.to_netcdf(outname,
                       encoding={'elev': encoding_params,
                                 'trend': encoding_params,
                                 'tide': encoding_params,
                                 'msea': encoding_params,
                                 'ss': encoding_params,
                                 'res': encoding_params,
                                 'skew_surge_magnitude': encoding_params,
                                 'skew_surge_lag': encoding_params,
                                 'tidal_elevation_maximum_over_tidal_cycle': encoding_params,
                                 'tidal_elevation_maximum_time_over_tidal_cycle': encoding_params,

                       })


def do_analysis(df_in, lat, cutoff=30, dt=1, constit='auto'):
    """
    Analyses the total water level contained in the 'elev' field of the dataframe df.
    The analysis includes extracting trend, tides, monthly mean sea level, storm surge
    and the remaining residuals.
    Inputs:
      - df_in: A Pandas Series containing the total water level data.
      - lat: The latitude of the location where the data were collected.
      - cutoff: The cutoff period to use to extract the storm surge in hours.
      - dt: The time resolution of the timeseries in hours. Default is 30 hours.
      - consist: UTide parameter that select the tidal constituents to extract.
           Default is 'auto' which should mean all.
    Outputs:
      - A Pandas DataFrame that contains all the extracted fields.
    """

    # Turn the input series into a dataframe whose colum 'elev' contains
    # the raw water level data.
    df = copy.deepcopy(df_in).rename('elev').to_frame()
    
    # Detrending but don't think there is much to detrend
    # Before detrending we store the position of all the gaps
    print('\tStoring NaNs position')
    f = np.where(np.isnan(df['elev'].values) == 1)
    
    # Get the detrended time series
    print('\tDetrending')
    df['et'] = df['elev']
    # Step 1: Store non nan data
    not_nan_ind = ~pd.isna(df['et'])
    x =  date2num(df['et'][not_nan_ind].index)
    y = df['et'][not_nan_ind].values
    # Step2: Fit linear curve
    m, b, r_val, p_val, std_err = linregress(x, y)
    # Step3: Store detrended data
    df['et'][not_nan_ind] = y - (m*x + b)
    # Step4: Store the trend
    df['trend'] = df['elev'] - df['et']

    #the tidal analysis
    print('\tTidal analysis')
    constituents = df.TideAnalysis._fit_tides(mag='et',
                                             args={'minimum SNR': 2,
                                                     'trend': False,
                                                     'latitude': lat,
                                                     'constit': constit
                                                })
    df = pd.concat([df, df.TideAnalysis\
                          ._tidal_elevation_from_constituents(constituents=constituents)\
                          .rename(columns={'tidal_elevation': 'tide'})],
                   axis=1)
    
    # Remove the tides
    df['et'] = df['et'] - df['tide']

    print('\tMonthly mean filtering')
    # Now we extract the mean sea level anomaly using a 30-day running average of the non-tidal residual (Haigh et al. 2014)
    # In order not to loose to much data when gaps exist in the data we allow the running average to return a value
    # as long as half of the data is non-nan
    df['msea'] = (df['et']).rolling(window=int(24*30/dt),
                                    min_periods=int(24*30/2),
                                    center=True).mean()

    # Remove the monthly mean sea level variations
    df['et'] = df['et'] - df['msea']

    # We extract the storm surge using a Lanczos lowpass filter
    # here to avoid loosing too much data we fill gaps that contain up to 3 nans using second order Akima spline interpolation.
    print('\tStorm surge filtering')
    df['ss'] = lanczos_filter(df['et'].interpolate(axis=0, limit=3, limit_area='inside', method='akima', order=2),
                              args={'window':30, 'type':'lanczos lowpas 2nd order'})

    # We subtract the storm surge to get the residuals
    df['res'] = df['et'] - df['ss']
    

    # We subtract that component to what is left of the signal and add the tide back (for the skew surge)
    df['et'] = df['et'] - df['msea'] + df['tide']

    # Get skew surge from tide + storm surge + residuals
    print('\tSkew surge filtering')
    df['et'] = df['res'].interpolate(axis=0, limit=3, limit_area='inside', method='akima', order=2) + df['ss'] + df['tide']
    df = pd.concat([df, df.TideAnalysis.skew_surge(mag='et',
                                                   args={'minimum SNR': 2,
                                                         'trend': False,
                                                         'latitude': lat,
                                                         'constituents':constituents,
                                                         'tide_dt': 60*60,
                                                        })],
                    axis=1)

    # Apply back initial mask
    for key in [k for k in df.keys() if not k in ['tide',
                                                  'skew_surge_magnitude',
                                                  'skew_surge_lag',
                                                  'tidal_elevation_maximum_over_tidal_cycle',
                                                  'tidal_elevation_maximum_time_over_tidal_cycle']]:
        df[key].values[f] = np.nan

    return df


def clean_linz(df0,
               dt=60,
               phasespace3d=True,
               despike=True,
               abs_threshold=3,
               detrend_btw_gap=False):
    """
    This routine allows to apply a number of cleaning processes to the data
    contained in the df0 dataframe. Those are targeted towards the cleaning
    of water level data.
    In order the following operations are applied if enabled:
      - Gap filling:
      - Resample of the data to hourly using a nearest neighbout approach
      - Apply a phase-space despiking (more details in
        https://github.com/calypso-science/Toto/blob/master/toto/filters/despike_phasespace3d.py)
      - Apply remove any point whose absolute distance to the mean of the timeseries
        is larger than abs_threshold parameter times the standard deviation of the timeseries.
        This can be seen as threshold-based despiking.
      - "detrend" the timeseries in a piecewise fashion (subtract its mean
        to each segment of data separated by nans).
      - Open an interactive window allowing to select manually points that
        might be left to delete.
    Arguments:
      - df0 (DataFrame): The dataframe to clean
      - dt (int): The sampling frequency in second corresponding to the DataFrame
      - phasespace3d (boolean): Whether to apply the phase-space despiking.
      - despike (boolean): Whether to apply the threshold-based despiking.
      - abs_threshold (float): Parameter to used as part of the threshold despiking
            method. Threshold = abs_threshold*stdv.
      - detrend_btw_gap (boolean): Whether to apply the piecewise detrending
    Outputs:
      A dataframe containing the clean data.
    """

    ### Start filtering
    print('\tFill gap')
    # Turns into an evenly spaced dataframe (interval dt seconds) starting at the first
    # non nan value of df, ending at the last and where missing values have been
    # replaced by the value of their nearest neighbour
    df1 = filled_gap(df0, dt)

    print('\tResample to hourly')
    df1 = df1.resample('1H').nearest()


    df1['filled'] = df1['elev']
    del df1['elev']

    
    if phasespace3d:
        print('\tRemove spike with phasespace3d')
        df1['phasespace'] = despike_phasespace3d(df1['filled'])
    else:
        df1['phasespace']=df1['filled'].copy()


    if despike:
        print('\tRemove spike with abs threshold (if value exceeds)')
        df1['despike'] = df1['phasespace'].copy()
        y = df1['phasespace'].to_numpy(copy=True)
        y = y-np.nanmean(y)
        thresh1 = abs_threshold*np.nanstd(y)
        ind = np.abs(y)> thresh1
        df1.loc[ind, 'despike'] = np.NAN
    else:
        df1['despike']=df1['phasespace'].copy()

    if detrend_btw_gap:
        print('\tDetrending spike')
        df1['detrend'] = f_detrend_between_gap(df1['despike'])
    else:
        df1['detrend'] = df1['despike']

    ## save to column 'elev' just not to keep it under 'detrend' and avoid confusing when 
    ## detrend_btw_gap is False
    df1['elev'] = df1['detrend']
    print('\tFinished automatic cleaning')

    return df1

def clean_manually(df0, df1):

    print('\tClean manually')
    answer=check.plot_graph([df0,df1],['elev','elev'],['Raw','Filter'],save=None,eventname='clean')
    deleted_time={}
    if 'delete' in answer:
        df1['clean']=delete_manually(df1['elev'],answer['delete'])
        
        for i,time in enumerate(answer['delete']):
            deleted_time[str(i+1)]='From %s to %s' % (num2date(time[0]).strftime('%Y-%m-%d %H:00'),num2date(time[1]).strftime('%Y-%m-%d %H:00'))

    else:
        df1['clean']=df1['elev']



    df=df1['clean'].copy().to_frame()
    df.rename(columns={'clean':'elev'},inplace=True)

    return df,df1,deleted_time

def delete_manually(ds,
                    times):
    """
    Inputs:
      - ds (DataFrame): Dataframe from which data needs deleting (replacing by nans)
      - times (list): List of pairs of date indices describing the intervals over
                      which the data has to be delete.
    Output:
      The input DataFrame with a number of data blanked out.
    """
    for time in times:
        ind = np.logical_and(date2num(ds.index) >= time[0],
                             date2num(ds.index) <= time[1])
        ds.values[ind] = np.nan
    return ds


def sort_dataset(df,**args):
    df.sort_index(inplace=True,**args)
    return df


def filled_gap(df,
               dt):
    """
    Returns evenly spaced dataframe (interval dt seconds) starting at the first
    non nan value of df, ending at the last.
    The values of the new dataframe are sourced from the input dataframe
    in a nearest neighbout fashion.
    """

    # Sort and remove nans
    df.sort_index(inplace=True)
    df = df.dropna()

    # Build empty dataframe starting at the first
    # non nan value of df, ending at the last and
    # regularly spaced (interval dt)
    dt = np.round(dt*1000) # *1000 + using ms as units is equivalent to 3 decimal rounding
    idx = pd.period_range(min(df.index),
                          max(df.index),
                          freq='%ims'%dt)\
            .to_timestamp()
    df0 = idx.to_frame(name='time').set_index('time')
    
    
    # Merge original dataframe with df0 to obtain evenly spaced dataframe
    # where missing values have been replaced by the value of their nearest neighbour
    df = pd.merge_asof(df0, df,
                       on='time',
                       direction='nearest',
                       tolerance=pd.Timedelta("%ims"%(int(dt*2))))\
           .set_index('time')
    return df


def f_detrend_between_gap(ds,
                          detrend='mean'):
    """
    Does piecewise detrending of the dataframe.
    Split ds in contiguous data segments (separated by nans)
    and detrend each segment separately ( here by default detrend is really subtracting the mean)
    """

    aa = np.isnan(ds.values).astype(int)

    data_groups = [list(group) for key, group in groupby(enumerate(aa), key=lambda ix: ix[1]) if key==0]

    for data_group in data_groups:
        s = data_group[0][0]
        t = data_group[-1][0]
        if detrend == 'mean':
            ds.values[s:t] = ds.values[s:t]-ds.values[s:t].mean()
        elif detrend == 'linear':
            # This would be a proper piecewise detrending
            ds.values[s:t] = detrend(ds.values[s:t])
        else:
            raise ValueError("f_detrend_between_gap only supports 'mean' and 'linear'"\
                             +" option for its detrend parameter.")
    
    return ds
