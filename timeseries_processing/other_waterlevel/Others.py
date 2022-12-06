import os,sys,glob
sys.path.append(os.path.join(os.path.dirname(__file__),'..'))
from core.reporting import PDF
from core import check

import copy
from toto.inputs.nc import NCfile
from core.toolbox import clean_linz,do_analysis,store_data_as_netcdf

def str2bool(v):
    if isinstance(v, bool):
        return v
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')

def add_pos(file):
    lon=0
    lat=0

    if 'Kawhia' in file:
        lon=174.82271804
        lat=-38.06586322
    elif 'LittleKaiteriteri' in file:
        lon=173.020234
        lat=-41.038576
    elif 'PortTaranaki' in file:
        lon=174.037381
        lat=-39.053697
    elif 'PoutoPoint' in file:
        lon=174.10540
        lat=-36.39637
    elif 'QueensWharf' in file:
        lon=174.7669504
        lat=-36.8405635
    elif 'Tarakohe' in file:
        lon=172.891660
        lat=-40.818808
    elif 'Thames' in file:
        lon=175.52079872
        lat=-37.12728521
    elif 'Whitianga' in file:
        lon=175.70882796
        lat=-36.83278872

    return lon,lat

def process(folderin,folderout,pdfout,tmpfolder,exclude,include,
    detrend_between_gap,phasespace3d,despike,cutoff,abs_threshold=3,
    diff_threshold=3):


    """ Reading and processing Other datafile

      Args:
        folderin (str): Path with all the Other file to process.
        folderout (str): Path to save all the processed LINZ file.
        pdfout(boolean): default is True. Either to save a pdf with the processing steps.
        tmpfolder(str): default is "/tmp" on linux machine. Temperory folder used during the PDF creation
        exclude (list): List of file that would not be processed. 
        detrend_between_gap (boolean): default False. Use this option do remove the mean between 2 gaps.
                                      This is usefull if the data are really "jumpy"
        phasespace3d (boolean): Use the phase space filter (recommanded). This will efficiently remove the spike in the timeseries
        despike (boolean): remove the spikes using scipy.signal.find_peaks.
        abs_threshold (float): threshold used in the scipy.signal.find_peaks function. value above this will be deleted
        diff_threshold (float): threshold used between 2 consecutive data point. value above this will be deleted
    """


    for file in glob.glob(os.path.join(folderin,'*.nc')):
        _,filein=os.path.split(file)
        prefix=filein.replace('_raw.nc','')

        if len(include)>0 and prefix not in include:
            continue


        if prefix in exclude:
            continue


        

        fileout=filein.replace('_all','_processed')
        pdflog=filein.replace('_all.nc','_log.pdf')

        print('Reading: %s' % filein)
        
        dfraw=NCfile(file)._toDataFrame()[0]
        lon,lat=add_pos(file)


        sens=False
        df0=copy.deepcopy(dfraw)
        del df0['time']
        # We have our sensor keep going with this
        # start doing the general clening by doing a phase shift cleaning 
        # and removing the major spikes with a big default value.
        # The user can also clean the timeseries Manually by selecting the period to delete with the mouse


        
        print('Cleaning ')

        df,df1,deleted_time=clean_linz(copy.deepcopy(df0),dt=60,
                                        phasespace3d=phasespace3d,
                                        despike=despike,
                                        abs_threshold=abs_threshold,diff_threshold=diff_threshold,
                                        detrend_btw_gap=detrend_between_gap)

        print('Do the analysis ')
        df = do_analysis(df, lat, cutoff=cutoff)


        store_data_as_netcdf(df,lon,lat,os.path.join(folderout,fileout))



        if pdfout:
            print('Saving log to %s' % pdflog)
            pdf = PDF('cleaning %s' % prefix,'Remy Zyngfogel')


            pdf.add_dict('Processing file',{'File in':filein,
                                        'File out':fileout})

            answer=check.plot_graph([dfraw],['elev'],["Raw data"],save=os.path.join(tmpfolder,'raw.png'),show=False)
            pdf.add_image('Raw data',os.path.join(tmpfolder,'raw.png'),'')
            os.system('rm %s' % os.path.join(tmpfolder,'raw.png'))

            if phasespace3d:
                answer=check.plot_graph([df1],['phasespace'],["Applied phasespace3D filter"],save=os.path.join(tmpfolder,'tmp1.png'),show=False)
                pdf.add_image('Applied phasespace3D filter',os.path.join(tmpfolder,'tmp1.png'),'')
                os.system('rm %s' % os.path.join(tmpfolder,'tmp1.png'))



            if despike:
                answer=check.plot_graph([df1],['despike'],["Despiked"],save=os.path.join(tmpfolder,'tmp0.png'),show=False)
                pdf.add_image('Applied despike, abs thresh = %i and diff thesh = %i' % (abs_threshold,diff_threshold),os.path.join(tmpfolder,'tmp0.png'),'')
                os.system('rm %s' % os.path.join(tmpfolder,'tmp0.png'))

            if detrend_between_gap:
                answer=check.plot_graph([df1],['detrend'],["Applied detrend between gap"],save=os.path.join(tmpfolder,'tmp2.png'),show=False)
                pdf.add_image('Applied detrend between gap',os.path.join(tmpfolder,'tmp2.png'),None)
                os.system('rm %s' % os.path.join(tmpfolder,'tmp2.png'))


            if deleted_time:
                pdf.add_dict('This time were deleted mananually',deleted_time)


            answer=check.plot_graph([df1],['clean'],["Final clean timeserie"],save=os.path.join(tmpfolder,'tmp3.png'),show=False)
            pdf.add_image('Clean timeserie',os.path.join(tmpfolder,'tmp3.png'),None)
            os.system('rm %s' % os.path.join(tmpfolder,'tmp3.png'))

            pdf.add_dict('Tidal analysis',{'Function': 'Utide'})
            pdf.add_dict('Mean sea level variation on monthly scale',{'Filter type': 'lanczos lowpas 2nd order','Windows':24*30})
            pdf.add_dict('Storm surge extraction',{'Filter type': 'lanczos lowpas 2nd order','Windows': cutoff})


            answer=check.plot_graph([df],['elev','tide','res','trend','msea','ss'],["SS analysis"],save=os.path.join(tmpfolder,'tmp4.png'),show=False)
            pdf.add_image('Clean timeserie',os.path.join(tmpfolder,'tmp4.png'),None)
            os.system('rm %s' % os.path.join(tmpfolder,'tmp4.png'))

            i=int(len(df.index)/2)
            xlim=[df.index[i],df.index[i+240]]
            ylim=[df['res'].min(),df['res'].max()]
            answer=check.plot_graph([df],['elev','tide','res','trend','msea','ss'],["SS analysis"],ylim=ylim,xlim=xlim,save=os.path.join(tmpfolder,'tmp5.png'),show=False)
            pdf.add_image('Clean timeserie zoomed',os.path.join(tmpfolder,'tmp5.png'),None)
            os.system('rm %s' % os.path.join(tmpfolder,'tmp5.png'))


            pdf.output(os.path.join(folderout,pdflog), 'F')

        print('Finished with %s' % filein)

if __name__ == "__main__":
    
    import argparse
    parser = argparse.ArgumentParser(prog='Linz.py', usage='%(prog)s raw_folder processed_folder [options]')
    ## main arguments
    parser.add_argument('raw_folder', type=str,help='Raw file folder')
    parser.add_argument('processed_folder', type=str,help='Processed file folder')

    ## options  
    parser.add_argument('-l','--log', type=str2bool,help='Export log',default=True,)
    parser.add_argument('-t','--temp', type=str,help='Tempory file folder',default='/tmp/',)
    parser.add_argument('-i','--include', type=str,nargs='+',help='list of file to include',default=[])
    parser.add_argument('-e','--exclude', type=str,nargs='+',help='list of file to exclude',default=[])
    parser.add_argument('-d','--detrend',type=str2bool,help='detrend betwwen gaps',default=False)
    parser.add_argument('-p','--phasespace3d',type=str2bool,help='despike using phasespace3d')#,default=True)
    parser.add_argument('-de','--despike',type=str2bool,help='despike using threshold',default=True)
    parser.add_argument('-c', '--cutoff', type=int, help='lanczos filter cut-off period (in hours) for storm surge', default=30)


    args = parser.parse_args()

    process(args.raw_folder,\
        args.processed_folder,\
        args.log, \
        args.temp,\
        args.exclude,
        args.include,
        args.detrend,
        args.phasespace3d,
        args.despike,
        args.cutoff)
