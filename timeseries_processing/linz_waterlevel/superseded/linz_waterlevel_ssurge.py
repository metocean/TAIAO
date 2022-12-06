import glob
import os
from datetime import datetime, timedelta

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import xarray as xr
from tqdm import tqdm

pd.options.display.max_rows = 15

################################################################################
#  This is a slightly modified version of the original qc_obs_ssurge.py script
#  from the msl_tools repo ()
#  This scripts read the LINZ water level daily csv files store in /net/datastor1/data/obs/nz qc's them to some extent, does gap filling using tidal data from ttide and extracts the storm surge.
#
#  The script relies on ttide. When "fixing" this script the version of ttide
#  used came from git@github.com:moflaher/ttide_py.git commit 46977c1a9ec09a4db8923180799e4d8585ebe088
#
#  The script also relies on a couple of functions (get_tidal_elevation and apply_filter) that used to
#  live in the hydro branch of the verify library. However, those relied on broken dependencies and were
#  copied to this file and fixed. In consequence, this script relies only on resonnably standard
#  python libraries.
#
#  Last edits were done by Seb Delaux
#
##########################################################################################################


from scipy.signal import butter, filtfilt
from ttide.t_tide import t_tide


def get_tidal_elevation(el, dt=1, lat=30, out_style='pandas', cons=None, verbose=True):
    """
    Returns tidal signal from harmonic analysis [using T_tide]
    """
    if verbose:
        output = 'screen'
    else:
        output = None

    elmean   = el.mean()
    elprime = (el - elmean)

    if cons:
        tideout = t_tide(elprime.to_numpy(), dt=dt, lat=lat, out_style='pandas', constitnames=cons, outfile=output)
    else:
        tideout = t_tide(elprime.to_numpy(), dt=dt, lat=lat, out_style='pandas', outfile=output)

    et = tideout['xout'].squeeze()

    if out_style == 'pandas':
        df =  pd.DataFrame({'et': et}, index=el.index)

    return df

def tidy_up(df, pddt='1H'):
    df = df.dropna()
    df = df.resample(pddt, closed='left').mean().interpolate()
    return df

def butter_lowpass(data, window, dt=1, order=5):
    """
    Inpulse response filter
    """

    fs   = (2*np.pi) / (dt*3600)
    nyq  = 0.5 * fs
    C    = 0.802

    window  = int( window / dt )
    highcut = (2*np.pi) / (window*3600)

    high = (highcut / nyq) / C # to prevent phase lag

    b, a = butter(order, high, btype='low')
    y    = filtfilt(b, a, data)

    return y

def apply_filter(df, window, dt, filter_type):
    df2 = df.copy()
    if filter_type == 'conv':
        df2 = pd.rolling_window(df2, window=window, win_type='hamming')
    elif filter_type == 'iir':
        for key in df2.keys():
            keym = df2[key].mean()
            df2[key] = butter_lowpass(df2[key] - keym, window, dt=dt, order=3)
            df2[key] += keym

    df2 = tidy_up(df2,'1min')

    return df2



def lat2msl(obsdir, df, site):
    # Converts data to mean sea level based on information contained in the station's readme
    f = open(os.path.join(obsdir, "{s}/{s}_readme.txt".format(s=site)),
             encoding='latin-1')
    readme = f.readlines()

    if site == 'CHST':
        df.elev -= 2.422
        return df
    
    if site == 'LOTT':
        df.elev -= 4.59
        return df

    if site == 'MNKT':
        df.elev -= 5.209
        return df

    if site == 'OTAT':
        df.elev -= 7.2
        return df

    if site == 'RBCT':
        df.elev -= 6.002
        return df

    if site == 'SUMT':
        df.elev -= (5.66 - 3.119)
        return df

    for line in readme:
        if 'SUMMARY OF TIDE GAUGE ZERO' in line:
            reference_mark = line.split(' ')[6]

    for line in readme:
        if 'LINZ geodetic code' in line and reference_mark in line:
            ref_datum = float(line.split(',')[-1].split()[0])

    start = None
    for idx, line in enumerate(readme):
        if 'SENSOR 41' in line:
            start = idx
            break

    if not start:
        for idx, line in enumerate(readme):
            if 'SENSOR 40' in line:
                start = idx
                break

    line = 'start'
    ref_gauge = dict()
    count = 1
    line = readme[start + count]
    while line != '\r\n' and line != '\n':
        m1 = datetime.strptime(line.split(' ')[0], '%b').month
        y1 = int(line.split(' ')[1])
        m2 = datetime.strptime(line.split(' ')[3], '%b').month + 1
        y2 = int(line.split(' ')[4][:-1])
        if m2 > 12:
            m2 = 1
            y2 += 1

        key = 'ref{}'.format(count)
        ref_gauge[key] = dict()
        ref_gauge[key]['value'] = float(line.split(' ')[5]) - ref_datum
        ref_gauge[key]['t1'] = datetime(y1, m1, 1)
        ref_gauge[key]['t2'] = datetime(y2, m2, 1)
        count += 1
        line = readme[start + count]

    for c in range(len(ref_gauge)):
        key = 'ref{}'.format(c+1)
        if key == 'ref1': # start
            if df.index[0] < ref_gauge['ref1']['t1']:
                ref_gauge['ref1']['t1'] = df.index[0]


            if ref_gauge['ref{}'.format(c+2)]['t1'] - ref_gauge[key]['t2'] > timedelta(days=1) and len(ref_gauge) > 1:
                ref_gauge[key]['t2'] = ref_gauge['ref{}'.format(c+2)]['t1'] - timedelta(hours=1)

        elif c == len(ref_gauge) - 1: # end
            if ref_gauge[key]['t2'] < df.index[-1]:
                ref_gauge[key]['t2'] = df.index[-1]

        else: # middle
            if ref_gauge['ref{}'.format(c+2)]['t1'] - ref_gauge[key]['t2'] > timedelta(days=1):
                ref_gauge[key]['t2'] = ref_gauge['ref{}'.format(c+2)]['t1'] - timedelta(hours=1)

        df.elev[ref_gauge[key]['t1'] : ref_gauge[key]['t2']] = df.elev[ref_gauge[key]['t1'] : ref_gauge[key]['t2']] - ref_gauge[key]['value']

    return df


def quality_control(df, site):
    if site == 'KAIT': # Kaikoura earthquake correction
        df.elev['2016-11-13 11:00':] =  df.elev['2016-11-13 11:00':] + 0.927
    # elif site == 'GIST':
    #     df.elev[(df.elev > 6) | (df.elev < 2.5)] = np.nan
    #     df.elev['2013-09-26':] =  df.elev['2013-09-26':] + 0.596
    if site == 'NAPT':
        df = df['2008-03-05':'2013-04-01']
        df.elev[(df.elev > 5) | (df.elev < 0)] = np.nan
    elif site == 'LOTT':
        df.elev[(df.elev > 5)] = np.nan
    elif site == 'OTAT':
        df.elev[(df.elev < 2)] = np.nan
    elif site == 'PUYT':
        df.elev[(df.elev < 4)] = np.nan
    elif site == 'AUCT':
        df.elev[(df.elev < 3)] = np.nan
    elif site == 'CHIT':
        df.elev[(df.elev > 4) | (df.elev < 1)] = np.nan
    elif site == 'CPIT':
        df.elev[(df.elev < 2)] = np.nan
    elif site == 'GBIT':
        df.elev[(df.elev < 4)] = np.nan
    elif site == 'GIST':
        df.elev[(df.elev < 2.5)] = np.nan
    
    
        

    return df


project_name = 'nz_ssurge'
obsdir = '/net/datastor1/data/obs/tide/linz/raw'
obsdir = '/home/sebastien/projects/Dev_SS/data/raw'

sites = [site for site in os.listdir(obsdir) if os.path.isdir(site)]
sites = ['AUCT']

def read_file(filepath):
    try:
        df = pd.read_csv(filepath, header=None, names=['station','time','elev'], parse_dates=['time'], usecols=[1,2], index_col='time')
        return xr.Dataset.from_dataframe(df)
    except:
        print("Failed to read ", filepath)
        sys.exit()

for site in tqdm(sites):
    print("Looking at site", site)
    # if site not in ['SUMT']:
    #     continue

    if site in ['RBCT', 'RFRT', 'ROBT', 'GBIT', 'NCPT', 'CHIT']: # outside continental NZ or dodgy data
        continue

    #if site in ['CHST']: # no reference datum in LINZ readme files => datum is 2.422
    #    continue

    if not os.path.isdir(os.path.join(obsdir, site)):
        print("Folder", site, "does not exist. Skipping site")
        continue

    filelist = glob.glob('{o}/{s}/????_??_?????*.zip'.format(o=obsdir, s=site))

    qcname = os.path.join(os.path.dirname(obsdir), 'qc', site + '.nc')#'.csv')

    if not os.path.isdir(os.path.dirname(qcname)):
        os.makedirs(os.path.dirname(qcname))

    # A bit ugly to switch between pandas and xarray twice but worth the speedup
    df = xr.concat([xr.Dataset.from_dataframe(pd.read_csv(filepath, header=None, names=['station','time','elev'], parse_dates=['time'], usecols=[1,2], index_col='time')) for filepath in filelist],dim='time').sortby('time').to_dataframe()

    # Do data qc
    df = quality_control(df, site)
    # Resample to hourly
    df = df.resample('1min').mean()
    # Convert from local datum to mean sea level
    df = lat2msl(obsdir, df, site)
    # Store location of missing data
    f = np.where(np.isnan(df.elev.values) == 1)
    # Fill time series gaps with mean
    df.fillna(df.elev.mean(), inplace=True)
    # Get tidal elevation from ttide for the duration of the time series
    et = get_tidal_elevation(df.elev, dt=1./60., verbose=False)
    # Replace the gap filled data with tidal elevation
    df.elev.values[f] = et.et.values[f]
    # Subtract mean to elevation data
    df.elev = df.elev - df.elev.mean()
    # Use filter to extract storm surge component
    dff = apply_filter(df, 40, 1./60., 'iir')
    df['ss'] = dff.elev
    df['tide'] = et
    # Re-masked values in gaps
    df.elev.values[f] = -999.
    df.ss.values[f] = -999.
    # slide time series
    df.index = df.index - pd.to_timedelta('12H')

    # remove first and last days because of IIR filter edges
    #df = df[24:-24] # assuming hourly interval
    df = df[24*60:-24*60] # assuming 1 minute interval

    # Create plot
 #   plt.figure()
 #   df.elev.plot(linewidth=2, alpha=0.5)
 #   df.ss.plot(linewidth=3)
 #   plt.title(site)
 #   plt.grid()
 #   plt.show()

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


