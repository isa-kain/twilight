import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from astropy.io import fits
from datetime import date, datetime, timedelta, timezone
from astroquery.jplhorizons import Horizons
import astropy.units as u
import pandas as pd
import xarray as xr
from tqdm import tqdm
import calendar
import time
from matplotlib.patches import Rectangle
from astropy.coordinates import AltAz, EarthLocation, SkyCoord, get_body, get_sun
from astropy.time import Time
import multiprocessing

from utils import *

ephempath = '/Users/isabelkain/Desktop/Twilight_/twilight-observing-tool/ephemeris-tables/'
toipath = '/Users/isabelkain/Desktop/Twilight_/twilight-observing-tool'

# Read in TOI list, ephemeris for this year
year = '2027'
toi_list = pd.read_csv(f'{toipath}/exofop_dec-30_starsep1500mas_2026aug02.csv', skiprows=2, on_bad_lines='skip')
# toi_list = pd.read_csv(f'{toipath}/toi_list_short.csv', skiprows=0, on_bad_lines='skip')
ephem = pd.read_csv(f'{ephempath}/twilight_ephemerides_keckII_{year}.csv')
fname = '2026aug02'


if __name__ == '__main__':

    # Set observatory location
    keck = EarthLocation.of_site('Keck')

    # Isolate unique TOI identifiers
    unique_TOIs = toi_list['TOI'].unique()
    print(f'{len(unique_TOIs)} unique TOIs.')
    np.save(f'{toipath}/vis_matrix_TOIs_{fname}.npy', unique_TOIs)

    # Set up empty matrix to catch boolean visibility results
    vis_matrix = np.full((len(unique_TOIs), len(ephem)), False)

    # Pull lists of TOI parameters to iterate over using multiprocess (shape: N_TOI x 365)
    toi_loc_list = [np.where(toi_list['TOI']==toi_num)[0][0] for toi_num in unique_TOIs]
    toi_coord_list = [SkyCoord(toi_list['RA'][loc], toi_list['Dec'][loc], unit=(u.hourangle, u.deg)) for loc in toi_loc_list]
    keck_site_list = [keck]*len(toi_loc_list)

    sunrise_temp = [datetime.strptime(ephem['sunrise_UTC'][j], '%Y-%b-%d %H:%M') for j in range(len(ephem))]
    sunrise_list = [sunrise_temp]*len(toi_loc_list) # shape: 423x365
    
    sunaz_temp = [ephem['sunrise_az'][j] for j in range(len(ephem))]
    sunaz_list = [sunaz_temp]*len(toi_loc_list) # shape: 423x365


    ########################
    ## Do multiprocessing ##
    ########################

    start_time = time.time()

    # Initialize pool
    pool = multiprocessing.Pool(4)

    # Do multiprocessing for checking if this target is up for each of 365 days of the year
    for i, result in enumerate(pool.starmap(multiprocess_toi_visibility, zip(toi_coord_list, sunrise_list, sunaz_list, keck_site_list))):
        vis_matrix[i, :] = result

    pool.close() # are these necessary?
    pool.join()

    # Calculate execution time
    end_time = time.time()
    execution_time = end_time - start_time
    print(f'Execution time: {execution_time:0.3f} s.') # typical time for first execution: 2.5s +- 0.5s.


    # Save result
    np.save(f'{toipath}/vis_matrix_{fname}.npy', vis_matrix)

    print('Done.')
