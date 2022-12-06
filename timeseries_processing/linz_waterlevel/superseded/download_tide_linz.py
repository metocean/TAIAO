import os
import re
import urllib.request as urllib

from bs4 import BeautifulSoup

#import pandas as pd

#needs to run as metocean user to be able to write into datastor!!!
rootdir = '/net/datastor1/data/obs/tide/linz/raw'
rootdir = '/home/sebastien/projects/Dev_SS/data/raw'
toplevel = 'http://apps.linz.govt.nz/'
ftp_url  = os.path.join( toplevel, 'ftp/sea_level_data/' )


sitesdict = {
'AUCT': {'x':174.783, 'y':-36.833, 'id':'auckland'},                            # http://apps.linz.govt.nz/ftp/sea_level_data/AUCT/AUCT_readme.txt
'RBCT': {'x':177.900, 'y':-29.283, 'id':'boat_cove-raoul_island'},              # http://apps.linz.govt.nz/ftp/sea_level_data/RBCT/RBCT_readme.txt
'ROBT': {'x':163.183, 'y':-77.033, 'id':'cape-roberts_antarctica'},             # http://apps.linz.govt.nz/ftp/sea_level_data/ROBT/ROBT_readme.txt
'CPIT': {'x':176.217, 'y':-40.917, 'id':'castlepoint'},                         # http://apps.linz.govt.nz/ftp/sea_level_data/CPIT/CPIT_readme.txt
'CHST': {'x':171.433, 'y':-41.900, 'id':'charleston'},                          # http://apps.linz.govt.nz/ftp/sea_level_data/CHST/CHST_readme.txt
'RFRT': {'x':177.900, 'y':-29.250, 'id':'fishing_rock-raoul_island'},           # http://apps.linz.govt.nz/ftp/sea_level_data/RFRT/RFRT_readme.txt
'GIST': {'x':178.033, 'y':-38.667, 'id':'gisborne'},                            # http://apps.linz.govt.nz/ftp/sea_level_data/GIST/GIST_readme.txt
'KAIT': {'x':173.700, 'y':-42.417, 'id':'kaikoura'},                            # http://apps.linz.govt.nz/ftp/sea_level_data/KAIT/KAIT_readme.txt
'GBIT': {'x':175.483, 'y':-36.183, 'id':'korotiti_bay-great_barrier_island'},   # http://apps.linz.govt.nz/ftp/sea_level_data/GBIT/GBIT_readme.txt
'LOTT': {'x':178.167, 'y':-37.550, 'id':'lottin_point'},                        # http://apps.linz.govt.nz/ftp/sea_level_data/LOTT/LOTT_readme.txt
'MNKT': {'x':174.517, 'y':-37.050, 'id':'manukau'},                             # http://apps.linz.govt.nz/ftp/sea_level_data/MNKT/MNKT_readme.txt
'NAPT': {'x':176.917, 'y':-39.483, 'id':'napier'},                              # http://apps.linz.govt.nz/ftp/sea_level_data/NAPT/NAPT_readme.txt
'NCPT': {'x':173.050, 'y':-34.417, 'id':'north_cape'},                          # http://apps.linz.govt.nz/ftp/sea_level_data/NCPT/NCPT_readme.txt
'CHIT': {'x':176.367, 'y':-44.033, 'id':'owenga-chatham_islands'},              # http://apps.linz.govt.nz/ftp/sea_level_data/CHIT/CHIT_readme.txt
'OTAT': {'x':170.650, 'y':-45.817, 'id':'port_chalmers'},                       # http://apps.linz.govt.nz/ftp/sea_level_data/OTAT/OTAT_readme.txt
'PUYT': {'x':166.583, 'y':-46.083, 'id':'puysegur'},                            # http://apps.linz.govt.nz/ftp/sea_level_data/PUYT/PUYT_readme.txt
'SUMT': {'x':172.567, 'y':-43.567, 'id':'sumner'},                              # http://apps.linz.govt.nz/ftp/sea_level_data/SUMT/SUMT_readme.txt
'TAUT': {'x':176.183, 'y':-37.650, 'id':'tauranga'},                            # http://apps.linz.govt.nz/ftp/sea_level_data/TAUT/TAUT_readme.txt
'WLGT': {'x':174.783, 'y':-41.283, 'id':'wellington'}                           # http://apps.linz.govt.nz/ftp/sea_level_data/WLGT/WLGT_readme.txt
}


sites = sitesdict.keys()
#sites = ['WLGT'] # For now focus on this site
sensor = '41' # '40'
sites = ['ROBT']
sensor = '00'

to_remove = []

for site in sites:
    sitedir = os.path.join(rootdir, site)
    gauge_url = os.path.join( ftp_url, site )

    # Start by downloading readme in case it was updated
    readmefile = urllib.URLopener()
    readmefile.retrieve(os.path.join(gauge_url, "{0}_readme.txt".format(site)),
                        os.path.join(sitedir, "{0}_readme.txt".format(site)))

    # get years of available data
    html_page = urllib.urlopen(gauge_url)
    years = []
    for link in BeautifulSoup(html_page, "html.parser").findAll('a'):
        year = re.findall('\d+', link.get('href'))
        if len(year):
            years.append(year[0])

    for year in years:
        filedir = os.path.join( gauge_url, year + '/'+sensor+'/')
        # get file URL level
        html_page = urllib.urlopen(filedir)
        zipfiles = []
        for link in BeautifulSoup(html_page, "html.parser").findAll('a'):
            fileurl = link.get('href')
            if '.zip' in os.path.basename(fileurl):
                zipfiles.append(toplevel+str(fileurl)[1:])

        for zipfile in zipfiles:
            # Local directory to save files to...
            if not os.path.isdir(sitedir):
                os.system('mkdir -p %s'%sitedir)

            outfile = os.path.join(sitedir, os.path.basename(zipfile))
            csvfile = os.path.basename(zipfile).replace('.zip', '.csv')

            # ensuring ww download daily files only
            if re.match(r".*_\d{7}.zip", os.path.basename(zipfile)):

                # Check if data has been downloaded, skip if so...
                if not os.path.isfile(outfile) and not os.path.isfile(outfile.replace('.zip', '.csv')):
                    print("Downloading", zipfile)
                    #sys.exit()
                    # Download tidal gauges data
                    linzfile = urllib.URLopener()
                    linzfile.retrieve(zipfile, outfile)

                    # Unzip and get rid of compressed files
                    #os.system( "unzip %s -d %s" %(outfile, sitedir))
                    #to_remove.append(outfile)

                    #csvfile = outfile.replace('.zip', '.csv')
                    #df = pd.read_csv(csvfile,header=None,index_col=1)
                    #df = pd.read_csv(zipfile,header=None,index_col=1)
                    #tstr = pd.to_datetime(df.index[0]).strftime('%Y%m%d')
                    #newcsvfile = csvfile.replace( re.findall('_(\d+).csv',csvfile)[0], tstr)
                    #os.rename(csvfile, newcsvfile)

# for zipfile in to_remove:
#     os.remove(zipfile)
