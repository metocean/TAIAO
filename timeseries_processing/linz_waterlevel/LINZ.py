import os,sys
from os.path import join
sys.path.append(join(os.path.dirname(__file__),'..')) 
from core_routines.reporting import PDF
from core_routines import check
import copy
from toto.inputs.linz import LINZfile
from core_routines.toolbox import clean_linz,clean_manually,do_analysis,store_data_as_netcdf
import glob
import matplotlib.pyplot as plt

def str2bool(v):
    if isinstance(v, bool):
        return v
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')

def read_linz(file, datum_height=None):
    """Read LINZ data file.
    Args:
        file (str): Full path of LINZ netCDF file.
        datum_height (float|None): LINZ station datum height (m).

    Returns:
        _type_: _description_
    """

    print(f"Reading: {file}")

    # Read the raw the NetCdf file and load both sensors
    dfraw = LINZfile(file, datum_height=datum_height)._toDataFrame()[0]
    lon = dfraw.longitude.drop_duplicates().values
    lat = dfraw.latitude.drop_duplicates().values

    return dfraw, lon, lat

def choose_sensor(dfraw, tmpfolder):

    #check which sensor we want to use, need to plot and click on the graph
    # we want to keep

    df0 = copy.deepcopy(dfraw)

    if 'elev40' in dfraw and 'elev41' in dfraw:

        answer = check.plot_graph([dfraw, dfraw], \
                                    ['elev40', 'elev41'], \
                                    ['Click on this graph to choose sensor: 40','41'],
                                    save=os.path.join(tmpfolder, 'raw_data.png'),
                                    eventname='choose')

        ## switch to right dataframe
        if answer['axes'] == '41':
            sensor = '41'
            del df0['elev40']
        else:
            sensor = '40'
            del df0['elev41']

        print(f'Switching to sensor {sensor}')
        df0.rename(columns={f'elev{sensor}':'elev'}, inplace=True)

    else:
        df0.rename(columns={'elev40':'elev'},inplace=True)

    del df0['time']

    return df0, sensor

def save_pdf(df, df1, pdflog, prefix, filein, fileout, cutoff, sensor, tmpfolder, phasespace3d, despike, \
    abs_threshold, detrend_between_gap, deleted_time, folderout):

    print('Saving log to %s' % pdflog)
    pdf = PDF('cleaning %s' % prefix,'Remy Zyngfogel')


    pdf.add_dict('Processing file',{'File in':filein,
                                'File out':fileout})

    pdf.add_dict('Tidal analysis',{'Function': 'Utide'})
    pdf.add_dict('Mean sea level variation on monthly scale',{'Filter type': 'lanczos lowpas 2nd order','Windows':24*30})
    pdf.add_dict('Storm surge extraction',{'Filter type': 'lanczos lowpas 2nd order','Windows':cutoff})

    pdf.add_page()
    if sensor:
        pdf.add_image('Choosing sensors',os.path.join(tmpfolder,'raw_data.png'),'Sensor # %s got selected' % sensor)
        os.system('rm %s' % os.path.join(tmpfolder,'raw_data.png'))
        pdf.add_page()

    if phasespace3d:
        answer=check.plot_graph([df1],['phasespace'],["Applied phasespace3D filter"],save=os.path.join(tmpfolder,'tmp1.png'),show=False)
        pdf.add_image('Applied phasespace3D filter',os.path.join(tmpfolder,'tmp1.png'),'')
        os.system('rm %s' % os.path.join(tmpfolder,'tmp1.png'))
        pdf.add_page()

    if despike:
        answer=check.plot_graph([df1],['despike'],["Despiked"],save=os.path.join(tmpfolder,'tmp0.png'),show=False)
        pdf.add_image('Applied despike, abs thresh = %i' % (abs_threshold), os.path.join(tmpfolder,'tmp0.png'),'')
        os.system('rm %s' % os.path.join(tmpfolder,'tmp0.png'))
        pdf.add_page()

    if detrend_between_gap:
        answer=check.plot_graph([df1],['detrend'],["Applied detrend between gap"],save=os.path.join(tmpfolder,'tmp2.png'),show=False)
        pdf.add_image('Applied detrend between gap',os.path.join(tmpfolder,'tmp2.png'),None)
        os.system('rm %s' % os.path.join(tmpfolder,'tmp2.png'))
        pdf.add_page()


    if deleted_time:
        pdf.add_dict('These times were deleted mananually', deleted_time)
        pdf.add_page()



    i=int(len(df.index)/2)
    xlim=[df.index[i],df.index[i+24*20]]
    #ylim=[df['res'].min(),df['res'].max()]

    answer=check.plot_graph([df1],['clean'],["Final clean timeserie"],save=os.path.join(tmpfolder,'tmp3.png'),show=False)
    pdf.add_image('Clean timeserie',os.path.join(tmpfolder,'tmp3.png'),None)
    os.system('rm %s' % os.path.join(tmpfolder,'tmp3.png'))

    answer=check.plot_graph([df1],['clean'],["Final clean timeserie zoomed)"],xlim=xlim,save=os.path.join(tmpfolder,'tmp3b.png'),show=False)
    pdf.add_image('',os.path.join(tmpfolder,'tmp3b.png'),None)
    os.system('rm %s' % os.path.join(tmpfolder,'tmp3b.png'))

    pdf.add_page()

    answer=check.plot_graph([df],['elev','tide',],["Tidal analysis"],save=os.path.join(tmpfolder,'tmp4.png'),show=False)
    pdf.add_image('Tidal analysis',os.path.join(tmpfolder,'tmp4.png'),None)
    os.system('rm %s' % os.path.join(tmpfolder,'tmp4.png'))

    answer=check.plot_graph([df],['elev','tide',],["Tidal analysis (zoomed)"],xlim=xlim,save=os.path.join(tmpfolder,'tmp4b.png'),show=False)
    pdf.add_image('',os.path.join(tmpfolder,'tmp4b.png'),None)
    os.system('rm %s' % os.path.join(tmpfolder,'tmp4b.png'))

    pdf.add_page()

    answer=check.plot_graph([df],['elev','trend','msea',],["Trend"],save=os.path.join(tmpfolder,'tmp5.png'),show=False)
    pdf.add_image('Trend analysis',os.path.join(tmpfolder,'tmp5.png'),None)
    os.system('rm %s' % os.path.join(tmpfolder,'tmp5.png'))
    answer=check.plot_graph([df],['elev','trend','msea',],["Trend (zoomed)"],xlim=xlim,save=os.path.join(tmpfolder,'tmp5b.png'),show=False)
    pdf.add_image('',os.path.join(tmpfolder,'tmp5b.png'),None)
    os.system('rm %s' % os.path.join(tmpfolder,'tmp5b.png'))

    pdf.add_page()


    answer=check.plot_graph([df],['elev','ss',],["Storm surge"],save=os.path.join(tmpfolder,'tmp6.png'),show=False)
    pdf.add_image('Storm surge analysis',os.path.join(tmpfolder,'tmp6.png'),None)
    os.system('rm %s' % os.path.join(tmpfolder,'tmp6.png'))
    answer=check.plot_graph([df],['elev','ss',],["Storm surge (zoomed)"],xlim=xlim,save=os.path.join(tmpfolder,'tmp6b.png'),show=False)
    pdf.add_image('',os.path.join(tmpfolder,'tmp6b.png'),None)
    os.system('rm %s' % os.path.join(tmpfolder,'tmp6b.png'))

    pdf.add_page()


    answer=check.plot_graph([df],['elev','skew_surge_magnitude','skew_surge_lag','tidal_elevation_maximum_over_tidal_cycle','total_water_level_maximum_over_tidal_cycle'],["Skew surge"],save=os.path.join(tmpfolder,'tmp7.png'),show=False)
    pdf.add_image('Skew surge analysis',os.path.join(tmpfolder,'tmp7.png'),None)
    os.system('rm %s' % os.path.join(tmpfolder,'tmp7.png'))
    answer=check.plot_graph([df],['elev','skew_surge_magnitude','skew_surge_lag','tidal_elevation_maximum_over_tidal_cycle','total_water_level_maximum_over_tidal_cycle'],["Skew surge (zoomed)"],xlim=xlim,save=os.path.join(tmpfolder,'tmp7b.png'),show=False)
    pdf.add_image('',os.path.join(tmpfolder,'tmp7b.png'),None)
    os.system('rm %s' % os.path.join(tmpfolder,'tmp7b.png'))

    #pdf.add_page()



    # answer=check.plot_graph([df],['elev','tide','res','trend','msea','ss','total_water_level_maximum_over_tidal_cycle'],["SS analysis"],ylim=ylim,xlim=xlim,save=os.path.join(tmpfolder,'tmp5.png'),show=False)
    # pdf.add_image('Clean timeserie zoomed',os.path.join(tmpfolder,'tmp5.png'),None)
    # os.system('rm %s' % os.path.join(tmpfolder,'tmp5.png'))


    pdf.output(os.path.join(folderout,pdflog), 'F')

def process(folderin, folderout, pdfout, tmpfolder, include, exclude, detrend_between_gap,
            phasespace3d, despike, cutoff, datum_height, abs_threshold=3):

    """ Reading and processing LINZ datafile

      Args:
        folderin (str): Path with all the LINZ file to process.
    folderout (str): Path to save all the processed LINZ file.
        pdfout(boolean): default is True. Either to save a pdf with the processing steps.
        tmpfolder(str): default is "/tmp" on linux machine. Temperory folder used during the PDF creation
        exclude (list): List of file that would not be processed. 
        include (list): List of file that would be processed. 
        detrend_between_gap (boolean): default False. Use this option do remove the mean between 2 gaps.
                                      This is usefull if the data are really "jumpy"
        phasespace3d (boolean): Use the phase space filter (recommanded). This will efficiently remove the spike in the timeseries
        despike (boolean): remove the spikes using scipy.signal.find_peaks.
        abs_threshold (float): threshold used in the scipy.signal.find_peaks function. value above this will be deleted
    """


    if '_raw.nc' in folderin:
        files=[folderin]
    else:
        files=glob.glob(os.path.join(folderin,'*.nc'))

    for file in files:
        _,filein=os.path.split(file)
        prefix=filein.replace('_raw.nc','')

        if prefix in exclude:
            continue

        if len(include)>0:
            if prefix not in include:
                continue

        

        fileout=filein.replace('_raw','_processed')
        pdflog=filein.replace('_raw.nc','_log.pdf')

        dfraw, lon, lat = read_linz(file, datum_height)

        df0, sensor = choose_sensor(dfraw, tmpfolder)

        # We have our sensor keep going with this
        # start doing the general clening by doing a phase shift cleaning 
        # and removing the major spikes with a big default value.
        # The user can also clean the timeseries Manually by selecting the period to delete with the mouse
        print('Cleaning ')

        df1 = clean_linz(copy.deepcopy(df0),
                         dt=60,
                         phasespace3d=phasespace3d,
                         despike=despike,
                         abs_threshold=abs_threshold,
                         detrend_btw_gap=detrend_between_gap)

        df, df1, deleted_time = clean_manually(copy.deepcopy(df0), df1)

        # Remove all the different signal (tide, trend, storm surge...)
        print('Do the analysis ')
        df = do_analysis(df['elev'], lat, cutoff=cutoff)

        # Save the new dataset into a NetCDF
        store_data_as_netcdf(df,lon,lat,os.path.join(folderout,fileout))

        # save a pdf with the parameters used and the different steps taken
        if pdfout:
            save_pdf(df, df1, pdflog, prefix, filein, fileout, cutoff, sensor, tmpfolder, \
                phasespace3d, despike, abs_threshold, detrend_between_gap, \
                deleted_time, folderout)

        print('Finished with %s' % filein)

if __name__ == "__main__":
    
    import argparse

    def parse_boolean(value):
        value = value.lower()

        if value in ["true", "yes", "y", "1", "t"]:
            return True
        elif value in ["false", "no", "n", "0", "f"]:
            return False

        return False

    parser = argparse.ArgumentParser(prog='Linz.py', usage='%(prog)s raw_folder processed_folder [options]')
    ## main arguments
    parser.add_argument('raw_folder', type=str,help='Raw file folder')
    parser.add_argument('processed_folder', type=str,help='Processed file folder')

    ## options  

    parser.add_argument('-l', '--log', type=parse_boolean, help='Export log', default=True)
    parser.add_argument('-t', '--temp', type=str, help='Tempory file folder', default='/tmp/')
    parser.add_argument('-e', '--exclude', type=str, nargs='+', help='list of file to exclude', default=[])
    parser.add_argument('-i', '--include', type=str, nargs='+', help='list of file to include', default=[])
    parser.add_argument('-d', '--detrend', type=parse_boolean, help='detrend betwwen gaps', default=False)
    parser.add_argument('-p', '--phasespace3d', type=parse_boolean, help='despike using phasespace3d', default=True)
    parser.add_argument('-de', '--despike', type=parse_boolean, help='despike using threshold', default=True)
    parser.add_argument('-c', '--cutoff', type=int, help='lanczos filter cut-off period (in hours) for storm surge', default=30)
    parser.add_argument('-da', '--datum_height', type=float, help='LINZ station datum height (m). \
        If defined, the given value will be used to apply sensor vertical drifting corrections and \
        reference water level around 0. If omitted, no sensor vertical drifting corrections will be \
        applied and the time average of sea level will be used to reference water level around 0. \
        Defaults to None.', default=None)

    args = parser.parse_args()
    process(args.raw_folder,\
        args.processed_folder,\
        args.log, \
        args.temp,\
        args.include,
        args.exclude,
        args.detrend,
        args.phasespace3d,
        args.despike,
        args.cutoff,
        args.datum_height)
