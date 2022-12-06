import os
from datetime import datetime
import pandas as pd
import requests
import xarray as xr
import tempfile
import sys


linz_gauges_metadata = {
    'AUCT': {'location': 'Auckland', 'tstart': datetime(2009, 3, 26), 'lon': 174.783, 'lat': -36.833 , 'sensors': [40, 41]},
    'CHIT': {'location': 'Owenga', 'tstart': datetime(2007, 12, 16), 'lon': -176.366, 'lat': -44.033, 'sensors': [40, 41]},
    'CHST': {'location': 'Charleston', 'tstart': datetime(2015, 7, 14), 'lon': 171.433, 'lat': -41.9, 'sensors': [40, 41]},
    'CPIT': {'location': 'Castlepoint', 'tstart': datetime(2009, 10, 6), 'lon': 176.216, 'lat': -40.91, 'sensors': [40, 41]},
    'GBIT': {'location': 'Korotiti Bay', 'tstart': datetime(2010, 7, 30), 'lon': 175.483, 'lat': -36.183, 'sensors': [40, 41]},
    'GIST': {'location': 'Gisborne', 'tstart': datetime(2008, 3, 10), 'lon': 178.033, 'lat': -38.666, 'sensors': [40, 41]},
    'KAIT': {'location': 'Kaikoura', 'tstart': datetime(2010, 5, 26), 'lon': 173.7, 'lat': -42.416, 'sensors': [40, 41]},
    'LOTT': {'location': 'Lottin Point', 'tstart': datetime(2008, 10, 10), 'lon': 178.166, 'lat': -37.55, 'sensors': [40, 41]},
    'MNKT': {'location': 'Manukau', 'tstart': datetime(2010, 7, 28), 'lon': 174.516, 'lat': -37.05, 'sensors': [40, 41]},
    'NAPT': {'location': 'Napier', 'tstart': datetime(2007, 9, 19), 'lon': 176.916, 'lat': -39.483, 'sensors': [40, 41]},
    'NCPT': {'location': 'North Cape', 'tstart': datetime(2008, 12, 24), 'lon': 173.05, 'lat': -34.416, 'sensors': [40, 41]},
    'OTAT': {'location': 'Port Chalmers', 'tstart': datetime(2010, 3, 1), 'lon': 170.64, 'lat': -45.816, 'sensors': [40, 41]},
    'PUYT': {'location': 'Puysegur Point', 'tstart': datetime(2009, 12, 14), 'lon': 166.583, 'lat': -46.083, 'sensors': [40, 41]},
    'RBCT': {'location': 'Boat Cove', 'tstart': datetime(2009, 5, 29), 'lon': -177.9, 'lat': -29.283, 'sensors': [40, 41]},
    'RFRT': {'location': 'Fishing Rock', 'tstart': datetime(2009, 5, 29), 'lon': 177.9, 'lat': -29.25, 'sensors': [40, 41]},
    'ROBT': {'location': 'Cape Roberts', 'tstart': datetime(2007, 11, 8), 'lon': 163.183, 'lat': -77.033, 'sensors': [0]},
    'SUMT': {'location': 'Sumner', 'tstart': datetime(2010, 8, 11), 'lon': 172.566, 'lat': -43.566, 'sensors': [40, 41]},
    'TAUT': {'location': 'Tauranga', 'tstart': datetime(2008, 7, 6), 'lon': 176.183, 'lat': -37.65, 'sensors': [40, 41]},
    'WLGT': {'location': 'Wellington', 'tstart': datetime(2007, 3, 23), 'lon': 174.783, 'lat': -41.283, 'sensors': [40, 41]}
}

def _read_file(filepath):
    try:
        df = pd.read_csv(filepath,
                         header=None,
                         usecols=[1,2],
                         names=['time','elev'],
                         parse_dates=['time'],
                         index_col='time')
        return xr.Dataset.from_dataframe(df)
    except:
        print("Failed to read ", filepath)
        sys.exit()


def _concat_files_to_xarray(files, sensor):
    dset = xr.concat([_read_file(filepath) for filepath in files],
                 dim='time').sortby('time')

    dset.time.attrs.update({'standard_name': 'time'})
    dset.elev.attrs.update({'standard_name': 'sea_surface_height_above_reference_ellipsoid',
                            'units': 'm',
                            '_FillValue': -999.})

    return dset.isel(time=~dset.to_dataframe().index.duplicated(keep='first')).expand_dims('sensor').assign_coords({'sensor': [int(sensor)]})


def _download_file(file_url, output_path, verbose=False):
    """
    Downloads a binary file from input url to the file system.
    
    Inputs:
      file_url    (string): The URL of the file to download
      output_path (string): The full path, including filename to download the file to.
      verbose     (string): Whether to print informations while downloading

    Return:
       (string): The full path where the downloaded data is store. None if no the
                 file was not found.

    """
    
    if verbose:
        print("Downloading ", file_url, "to", output_path)

    r = requests.get(file_url)
    with open(output_path, 'wb') as f:
        f.write(r.content)

    if r.status_code == 200:
        if verbose:
            print("Dowloaded "+output_path+" dowloaded successfully")
        return output_path
    elif (r.text == "Not found"):
        print(file_url+" not found")
        os.remove(output_path)
        return None
    else:
        print("Failed downloading on "+output_path)
        print(r.text, r.status_code)
        os.remove(output_path)
        raise


class LINZ_Site(object):

    def __init__(self,
                 site_name,
                 filename=None,
                 tstart=None,
                 tend=datetime.now(),
                 verbose=False):
        """
          Instantiates a LINZ_Site object that allow to access water elevation data from the
          LINZ website. By default the object will allow to download the full dataset for a 
          given location. However, it is also possible to provide a filename that points to
          a previously downloaded version of the timeseries. The class will allow to update
          that data rather than redownloading everything.

          Inputs:
            site_name (string): A four letter name corresponding to the site of interest.
                                Available names can be found in the linz_gauges_metadata
                                at the top of this file or on the LINZ website
                                https://sealevel-data.linz.govt.nz/
            filename (string):  The name of file that already contains previously downloaded
                                data for the site that are to be either used of updated.
            tstart (datetime):  Date when to start to download the data from. Default is either
                                the start of the available record or the end of the record
                                already available on file.
            tend (datetime):    The end date to consider when downloading data.
                                Default is now.

          Return:
            The instantiated object
        """
        

        self.filename = filename
        self.site_name = site_name
        self.download_folder = tempfile.mkdtemp()
        self.verbose = verbose

        if site_name not in linz_gauges_metadata:
            print('Invalid site name.')
            print('Supported site names are: ', ','.join(linz_gauges_metadata.keys()))
                  
        site_info = linz_gauges_metadata[site_name]
                  
        self.location = site_info['location']
        self.longitude = site_info['lon']
        self.latitude = site_info['lat']
        self.tstart = site_info['tstart']
        self.sensors = site_info['sensors']

        if tstart:
            self.tstart = tstart
        self.tend = tend

        self.dset = None
        if filename and os.path.isfile(filename):
            self.dset = xr.open_dataset(filename)
            self.sensors = self.dset.sensor.values
            self.tstart = pd.to_datetime(self.dset.sortby('time').time.values[-1]).to_pydatetime()


    def download_root_url(self, sensor, year):
        """
          Returns the base dowload url for a given year and sensor

          Input:
            sensor (int): The id of the sensor
            year   (int): The year

           Return:
            (string) the url to the data
        """
        
        root_url_template = 'https://sealevel-data.linz.govt.nz/tidegauge/%s/%04d/%02d/'

        return root_url_template % (self.site_name,
                                    year,
                                    sensor)
        
    def download_filename(self, sensor, year, day):
        """
          Return the name of the file to download based on sensor, year and day.

          Input:
            sensor (int): The id of the sensor
            year   (int): The year
            day    (int): The day

           Return:
            (string) the name of the data file
        """
        
        filename_template = '%s_%02d_%04d%03d'

        return filename_template%(self.site_name,
                                  sensor,
                                  year, day)
        
    def _download_data(self):
        """
          Download data from the LINZ website from tstart and defined when
          instantiating the object to tend.

          Return:
            xarray dataset: The downloaded data in the shape of a xarray dataset.
        """

        tstart = self.tstart
        tend = self.tend

        data = {}
        for sensor in self.sensors:

            print(f'Downloading data from sensor {sensor}')
            download_folder = self.download_folder if not self.download_folder is None\
                                                   else tempfile.TemporaryDirectory().name

            downloaded_files = []

            # Loop over all years the site was acive for
            for year in range(tstart.year, tend.year+1):
                print("Year: ", year)
        
                root_url = self.download_root_url(sensor=sensor,
                                                  year=year)

                start_day = tstart.timetuple().tm_yday if year == tstart.year else 1
                # Last day of the year we expect data for 365/366 or less if it is this year
                end_day = tend.timetuple().tm_yday if tend.year == year\
                                                             else datetime(year,12,31).timetuple().tm_yday

                # Loop over days
                for day in range(start_day, end_day+1):
                    # Expected name for file to download
                    fname = self.download_filename(sensor=sensor,
                                                   year=year,
                                                   day=day)
                    
                    # Check if already downloaded
                    if os.path.exists(os.path.join(download_folder, fname+".zip")):
                        downloaded_files.append(os.path.join(download_folder, fname+".zip"))
                    elif os.path.exists(os.path.join(download_folder, fname+".csv")):
                        downloaded_files.append(os.path.join(download_folder, fname+".csv"))
                    else:
                        try:
                            # Try to download
                            downloaded_file = _download_file(os.path.join(root_url, fname)+".zip",
                                                             os.path.join(download_folder, fname+".zip"),
                                                             verbose=self.verbose)

                            if downloaded_file is not None:
                                downloaded_files.append(downloaded_file)
                            
                        except:
                            if self.verbose:
                                print("Download Failed")
                                print(os.path.join(root_url, fname)+".zip")
                            pass
                        
            if len(downloaded_files) > 0:
                data[sensor] = _concat_files_to_xarray(files=downloaded_files,
                                                      sensor=sensor)

        return xr.concat([data[sensor] for sensor in self.sensors],
                         dim='sensor')

    def get_data(self):
        """
          Return the downloaded/updated timeseries of water level data.

          Return:
            The downloaded/updated timeseries of water level data in the shape of
            an xarray dataset
        """

        if self.dset:
            dset = xr.concat([self.dset.drop_vars(['longitude',
                                                   'latitude']),
                              self._download_data()],
                             dim='time')
        else:
            dset = self._download_data()

        self.dset = dset.assign({"longitude": (('site',),[self.longitude]),
                                 "latitude":(('site',),[self.latitude])})\
                        .assign_attrs({'site_name': self.site_name,
                                       'location_name': self.location})

        return self.dset

    def to_netcdf(self, filename):
        """
          Store the data in a netcdf file

          Input:
            filename (string): Full path where to store the data
        """

        self.dset.to_netcdf(filename,
                            encoding={'elev': {"dtype": "short",
                                               "scale_factor": 0.001}})
    
    def get_readme(self, filename):
        """Download readme/description ascii file for a specific site.

        Args:
            filename (str): Full path to store ascii file
        """

        url = f"http://sealevel-data.linz.govt.nz/api/getGaugeData/{self.site_name}"
        site_info = requests.get(url)
        site_desc = site_info.json()['description']
        output = open(filename, "w")
        output.write(site_desc)
        output.close()


if __name__ == "__main__":

    site_name = 'AUCT'
    input_file = '/home/seb/metocean/ssurge/storm_surge_data/nz_tidal_gauges/linz/raw/AUCT_raw.nc'
    output_file = '/home/seb/metocean/ssurge/storm_surge_data/nz_tidal_gauges/linz/raw/AUCT_raw_updated.nc'
    readme_file = '/home/seb/metocean/ssurge/storm_surge_data/nz_tidal_gauges/linz/raw/AUCT_readme.txt'

    site = LINZ_Site(site_name=site_name, filename=input_file)

    site.get_data()
    site.to_netcdf(output_file)
    site.get_readme(readme_file)
