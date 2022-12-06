import toto
from toto.inputs.nc import NCfile
from distribution_toolbox import (
    get_peaks,
    do_EVA_water,
    plot_kde
)
import numpy as np
import os
import glob
from toto.core.make_table import create_table

root = '/home/remy/projects/019_stormsurge/storm_surge_data/nz_tidal_gauges/'
subfolders = ['linz', 'other', 'uhslc']

# Parameters to find peaks
threshold_type = 'percentile' # or 'value'
thresh = 95.0
min_peak_over_threshold = 3
min_time_interval = 24 # in hours
time_blocking = 'Annual'

# Parameter for the distribution
distribution = 'GPD' # can be 'Weibull','Gumbel','GPD' or 'GEV
method = 'ml' # 'pkd','pwm','mom' or 'ml'
surge_mode = 'Positive' # 'Negative'

def get_distribution(df):
    # prepare the dataframe
    dfout = df.copy()
    dfout['positive_surge'] = df['ss']
    dfout['negative_surge'] = df['ss'] * -1.
    dfout['positive_tide'] = df['tide']
    dfout['negative_tide'] = df['tide'] * -1.
    sint = ( df.index[1] - df.index[0] ).total_seconds()

    # choose which surge we are studying
    if surge_mode == 'Positive':
        surge = 'positive_'
    elif surge_mode == 'Negative':            
        surge = 'negative_'

    pks_opt = {}
    if threshold_type == 'percentile':
        sort_data = np.sort(np.abs(dfout['ss'].values))
        sort_data = sort_data[~np.isnan(sort_data)]
        pks_opt['height'] = sort_data[int(np.round(len(sort_data)*(thresh/100)))]

    else:
        pks_opt['height'] = thresh

    pks_opt['distance']=min_time_interval*(sint/3600)


    peaks_index = get_peaks(dfout,surge+'surge',
                            time_blocking=time_blocking,
                            peaks_options=pks_opt,
                            min_peak=min_peak_over_threshold)

    # in case there is some month with no peaks
    keys=list(peaks_index.keys())
    for key in keys:
        if len(peaks_index[key]) == 0:
            peaks_index.pop(key)

    eva_stats = do_EVA_water(peaks_index,
                             dfout,
                             surge+'surge',
                             surge+'tide',
                             distribution,
                             method,
                             time_blocking,
                             pks_opt['height'],
                             tide_join_prob=False)
    scale = eva_stats['Annual']['Omni'][surge+'surge']['scale']
    shape = eva_stats['Annual']['Omni'][surge+'surge']['shape']
    return scale,shape


for subfolder in subfolders:
    folder = os.path.join(root, subfolder, 'processed')
    files = glob.glob(os.path.join(folder, '*.nc'))
    shape_table = np.empty((len(files)+1,3), dtype="object")
    shape_table[0,0] = 'Filename'
    shape_table[0,1] = 'Scale'
    shape_table[0,2] = 'Shape'
    for i, file in enumerate(files):
        print(file)
        shape_table[i+1, 0] = os.path.split(file)[-1]
        df = NCfile(file)._toDataFrame()[0]
        scale,shape = get_distribution(df)
        shape_table[i+1,1] = '%.8f' % scale
        shape_table[i+1,2] = '%.8f' % shape

    create_table('distribution_shapes.xlsx',
                 subfolder,
                 shape_table)
