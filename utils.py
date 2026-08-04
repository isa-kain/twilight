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
from matplotlib.patches import Rectangle
from astropy.coordinates import AltAz, EarthLocation, SkyCoord, get_body, get_sun
from astropy.time import Time


ephempath = '/Users/isabelkain/Desktop/Twilight_/twilight-observing-tool/ephemeris-tables'
toipath = '/Users/isabelkain/Desktop/Twilight_/twilight-observing-tool'


def parse_datestring(datestr=None):
    '''
    Ingests string of unknown datetime format and checks it against a number of datetime formats:
    
    ** Year, month and day specified: **
    YYY-MM-DD > %Y-%m-%d
    YY-MM-DD > %y-%m-%d
    MM-DD-YYYY > %m-%d-%Y
    MM-DD-YY > %m-%d-%y
    mnth DD YYYY > %b %d %Y
    mnth DD YY > %b %d %y
    mnth-DD-YYYY > %b-%d-%Y
    mnth-DD-YY > %b-%d-%y
    Month DD YYYY > %B %d %Y
    Month DD YY > %B %d %y
    YYYY-mnth-DD > %Y-%b-%d
    YYYY mnth DD > %Y %b %d
    YY-mnth-DD > %y-%b-%d
    YY mnth DD > %y %b %d
    
    
    ** Month and day only specified: **
    MM-DD > %m-%d
    mnth DD > %b %d
    month DD > %B %d
    
    ** User specifies current day: **
    Today
    Tonight
    Now
    
    Returns:
    date :: datetime object
    year :: int, e.g. 2026
    '''
    
#     print(':: DEBUG ::', datestr, type(datestr), datestr==None)
    
    # Write out datetime formats to check user input against
    list_of_dtformats_withyear = ['%Y-%m-%d', '%y-%m-%d', '%m-%d-%Y', '%m-%d-%y', 
                                  '%b %d %Y', '%b %d %y', '%b-%d-%Y', '%b-%d-%y',
                                  '%B %d %Y', '%B %d %y', '%Y-%b-%d', '%Y %b %d', 
                                  '%y-%b-%d', '%y %b %d']
    list_of_dtformats_woutyear = ['%m-%d', '%b %d', '%B %d']
    list_of_dtformats_today = ['today', 'tonight', 'now']
    
    
    # Create datetime object for current day
    now = datetime.now()
    
    # If datestr not specified, take today's date
    if datestr==None:
        date = now # contains current hour, minute, second, ect – shouldn't matter
        year = now.strftime('%Y')
        return date, year
    
    # Check entered datestring against strptime formats (year user-specified)
    for fmt in list_of_dtformats_withyear:    
        try:
            date = datetime.strptime(datestr, fmt)
            year = date.strftime('%Y') 
            print('User-entered observation date:', date)
        except ValueError:
            pass

        
    # Check entered datestring against strptime formats (year NOT user-specified)
    for fmt in list_of_dtformats_woutyear:    
        try:
            date = datetime.strptime(datestr, fmt)
            year = now.strftime('%Y') 
            date = date.replace(year=int(year))
            print('User-entered observation date:', date)
        except ValueError:
            pass
    
    # Maybe your user is overeager and specifies that they TODAY
    if any( fmt.casefold() in datestr.casefold() for fmt in list_of_dtformats_today ):
        date = now # contains current hour, minute, second, ect – shouldn't matter
        year = now.strftime('%Y')
    
    # If no formats matched, return ValueError; else, return datetime object and year string
    try:
        return date, year
    except NameError:
        raise ValueError('User-entered date format not recognized. Try again!')


def read_solar_ephemeris(datestr=None):
    '''
    Read pre-generated lookup table of solar ephemerides from Keck II. 
    User can input an optional UTC date string (YYY-MM-DD) to pull ephemerides
    from this date; otherwise, the code pulls ephemerides from today.
    '''
    
    # Parse date if given by user; otherwise take today's date
    dt, year = parse_datestring(datestr)
    date = dt.strftime('%Y-%m-%d')
#     print(':: DEBUG ::', date, type(date), year, type(year))
        
        
    # Read in solar ephemeris table from matching year
    ephem = pd.read_csv(f'{ephempath}/twilight_ephemerides_keckII_{year}.csv')
    
    
    # Identify row of ephemeris table that matches date
    try: 
        
        ix = ephem.index[ephem['Date'] == date].tolist()

        try: 
            assert len(ix)==1
            ix = ix[0]
        except: 
            print(f'Multiple rows for the same date. Something went wrong for the {date} ephemeris!')
            return 0

    except: 
        print(f'Solar ephemeris not available for {date}. Go make one.')
        return 0
    
    
    # Parse sunrise/sunset parameters from matching date
    sunset_UTC = ephem['sunset_UTC'][ix]
    sunset_az = ephem['sunset_az'][ix]
    sunset_el = ephem['sunset_el'][ix]
    sunset_ra = ephem['sunset_ra'][ix]
    sunset_dec = ephem['sunset_dec'][ix]
    sunset_sidereal = ephem['sunset_sidereal'][ix]
    
    sunrise_UTC = ephem['sunrise_UTC'][ix]
    sunrise_az = ephem['sunrise_az'][ix]
    sunrise_el = ephem['sunrise_el'][ix]
    sunrise_ra = ephem['sunrise_ra'][ix]
    sunrise_dec = ephem['sunrise_dec'][ix]
    sunrise_sidereal = ephem['sunrise_sidereal'][ix]
    
    print('Sunset time:', sunset_UTC)
    print('Sunrise time:', sunrise_UTC)
    
    return sunset_UTC, sunset_az, sunset_el, sunset_ra, sunset_dec, sunset_sidereal, \
            sunrise_UTC, sunrise_az, sunrise_el, sunrise_ra, sunrise_dec, sunrise_sidereal
    


def keckII_pointing_limits(az, el):
    '''
    Returns if given az, el is within Keck II telescope pointing constraints.
    See documentation: https://www2.keck.hawaii.edu/inst/common/TelLimits.html
    '''
    
    # Check against nasmyth deck pointing limit
    if (az >= 185.3) and (az <= 332.8):
        if el >= 36.8:
            return True
        else:
            return False
        
    # Check against shutter vignetting
    else:
        if el >= 18.0:
            return True
        else:
            return False
        
        
def twilight_pointing_limits(obj_az, sun_az, az_sep=90.):
    '''
    Check if object az is >=90˚ from the sun. (az, sun_az units in degrees)
    '''

    # Find difference between sun, telescope az angles
    az_diff = np.abs(obj_az-sun_az) % 360.
#     print(':: DEBUG :: ', az_diff)
#     print(':: DEBUG :: twilight_pointing_limits() sun_az / obj_az / az_diff:', sun_az, obj_az, az_diff)
                
    # Is object az > minimum separation angle from the rising sun?
    if az_diff <= az_sep: 
        return False
    else: 
        return True
    
    
    

def is_target_up_morn_twi(sunrise_dt, sunrise_az, obj_coord, site):
    '''
    sunrise_dt :: datetime object, datetime of sunrise at site
    sunrise_az :: az angle of sunrise (deg)
    obj_coord ::  SkyCoord (ICRS) of target object
    site :: EarthLocation object of observing site
    '''

    # Sample target visibility at 5 10-minute increments during hour before sunrise 
    # (twilight_samptimes increases in time, e.g. sunrise is last position in list)
    twilight_samptimes = [sunrise_dt - timedelta(minutes=10*x) for x in range(4,-1,-1)]
#     print(':: DEBUG :: ', twilight_samptimes)


    # Check if target is within Keck II pointing limits at each sample time
    sampled_visibility = np.full(len(twilight_samptimes), False)

    for j, time in enumerate(twilight_samptimes):

        # What is az/el location of TOI at given time? 
        obj_coord_azel = obj_coord.transform_to(AltAz(obstime=Time(time), location=site))
        obj_az = obj_coord_azel.az.value
        obj_el = obj_coord_azel.alt.value
#         print(':: DEBUG :: is_target_up_morn_twi() object azel:', obj_az, obj_el)

        # Is it within telescope pointing constraints?
        kvis = keckII_pointing_limits(obj_az, obj_el) # returns True/False
        
        # Is it within twilight pointing constraints?
        tvis = twilight_pointing_limits(obj_az, sunrise_az) # default az_sep=90.
        
        # At this sampled time, is it within both Keck & twilight constraints?
        sampled_visibility[j] = kvis and tvis


    # Is object up for at least 3 out of 5 10-min sampled periods (i.e. up for >60% of the 40 mins before sunrise)
    # AND is the object up at sunrise (to prevent looking at a target and having it set early)
#     is_up = (np.sum(sampled_visibility) >= 3.) and (sampled_visibility[-1]==True) # True/False
    is_up = (np.sum(sampled_visibility) >= 3.)  # True/False -- FIXME is this the right litmus?
    
#     if np.sum(sampled_visibility)==0:
#         print(':: DEBUG :: this TOI is not up at all:', obj_coord, np.sum(sampled_visibility), np.shape(sampled_visibility))
    
    return is_up


def multiprocess_toi_visibility(toi_coord, sunrise_list, sunaz_list, site): 
    '''
    [thing it does]
    
    Inputs:
    toi_coord :: SkyCoord object specifying TOI (ra, dec)
    sunrise_list :: 
    site :: EarthLocation object specifying observing site.
    
    Returns:
    vis_matrix_row :: (1 x N_TOIs) boolean array. If True, target is visible during morning twilight.
    '''
    
    vis_matrix_row = np.full(len(sunrise_list), False)
    
    for j in range(len(sunrise_list)): # for each day in given year
        
        sunrise_dt = sunrise_list[j]
        sunrise_az = sunaz_list[j]
        is_up = is_target_up_morn_twi(sunrise_dt, sunrise_az, toi_coord, site=site) # sunrise_dt, sunrise_az, obj_coord, site

        vis_matrix_row[j] = is_up
        
    return vis_matrix_row

    

    
###################################################################
## Solar System stuff
###################################################################

    
def query_JPL_horizons(start, end):
    '''
    start: datetime object, year month day 03:00:00
    start: datetime object, year month day 18:00:00
    '''

    # Set epoch to query
    epoch =  {'start':start.strftime('%Y-%m-%d %H:%M:%S'), 'stop':end.strftime('%Y-%m-%d %H:%M:%S'), 'step':'1m'}

    
    # Send query
    obj = Horizons(id='sun', location='T17', epochs=epoch) # t17 = Keck II
    result = obj.ephemerides(skip_daylight=True)

    if not (result['solar_presence'][0]=='C') & (result['solar_presence'][-1]=='C'):
        print('FUCKED UP!', epoch)

        
    # Extract sunset params
    sunset_time = result['datetime_str'][0] # datetime string, e.g. '2025-Jan-01 04:04'
    sunset_az = result['AZ'][0] # deg
    sunset_el = result['EL'][0] # deg
    sunset_ra = result['RA_app'][0] # deg
    sunset_dec = result['DEC_app'][0] # deg
    sunset_sidereal = result['siderealtime'][0] # hours

    
    # Extract sunrise params
    sunrise_time = result['datetime_str'][-1] # datetime string, e.g. '2025-Jan-01 04:04'
    sunrise_az = result['AZ'][-1] # deg
    sunrise_el = result['EL'][-1] # deg
    sunrise_ra = result['RA_app'][-1] # deg
    sunrise_dec = result['DEC_app'][-1] # deg
    sunrise_sidereal = result['siderealtime'][-1] # hours

    
    # Return all params as single row
    return start.strftime('%Y-%m-%d'), sunset_time, sunset_az, sunset_el, sunset_ra, sunset_dec, \
            sunset_sidereal, sunrise_time, sunrise_az, sunrise_el, sunrise_ra, sunrise_dec, sunrise_sidereal



def grab_ss_trace(objname, site, sunset_UTC, sunrise_UTC):
    '''
    objname :: string name of major SS body (e.g. Jupiter, Neptune, NOT any of the moons)
    site :: EarthLocation object
    sunset_UTC :: datetime object
    '''
    
    # Generate list of times throughout night to check where SS object is in the sky
    sample_times = Time( pd.date_range(start=sunset_UTC, end=sunrise_UTC, periods=30).to_numpy() )
    obj_az = np.zeros(len(sample_times))
    obj_el = np.zeros(len(sample_times)) 


    # Grab object az,el location at each timestamp
    for i, time in enumerate(sample_times):

        obj_coord = get_body(objname, time=time, location=site) 
        obj_coord_altaz = obj_coord.transform_to(AltAz(obstime=time, location=keck))

        obj_az[i] = obj_coord_altaz.az.value
        obj_el[i] = obj_coord_altaz.alt.value
        
        # If object not effectively visible with Keck II, mask values with NaNs
        if not keckII_pointing_limits(obj_az[i], obj_el[i]):
            obj_el[i] = np.nan

#     # Mask timestamps with NaNs where el < 18deg
#     obj_el[obj_el < 18.] = np.nan
    
    # Return location
    return obj_az, obj_el, sample_times




def make_azel_plot(sunset_UTC, sunset_az, sunrise_UTC, sunrise_az, obj_az, obj_el, trace_times, obj_label='Target'):
    
    # Initialize plot
    fig, ax = plt.subplots(1, 1, figsize=(5, 7), subplot_kw={'projection': 'polar'}, layout='constrained')
    plt.tight_layout()
    
    ax.set_rmax(0.)
    ax.set_rmin(90.)
    ax.grid(True)
    ax.set_theta_zero_location('N')

    ax.annotate('N', (np.deg2rad(10.), -6.), xycoords='data', fontsize=18, annotation_clip=False, fontweight='bold')
    ax.annotate('S', (np.deg2rad(185.), -12.), xycoords='data', fontsize=18, annotation_clip=False, fontweight='bold')
    ax.annotate('E', (np.deg2rad(85.), -14.), xycoords='data', fontsize=18, annotation_clip=False, fontweight='bold')
    ax.annotate('W', (np.deg2rad(275.), -6.), xycoords='data', fontsize=18, annotation_clip=False, fontweight='bold')

    # Calculate Keck II pointing limits
    linspace_az = np.linspace(0., 2.*np.pi, 100)

    tel_constraints_el = np.full_like(linspace_az, 18.0)
    tel_constraints_el[(linspace_az >= np.deg2rad(185.3)) & (linspace_az <= np.deg2rad(332.8))] = 36.8

    horizon_el = np.zeros_like(linspace_az)
    
    # Plot Keck II pointing limits
    ax.plot(linspace_az, tel_constraints_el, ls='--', color='gray')
    ax.fill_between(linspace_az, horizon_el, tel_constraints_el, color='gray', alpha=0.3, label='Keck II limits')
    
    # Calculate twilight pointing constraints
    
    az_lolim = sunrise_az - 45.
    az_uplim = sunrise_az + 45.

    if az_lolim < 0:
        az_lolim += 360.
    if az_uplim > 360.:
        az_uplim %= 360.

    twi_constraints_el = np.full_like(linspace_az, 0.0)
    twi_constraints_el[(linspace_az > np.deg2rad(az_lolim)) & (linspace_az < np.deg2rad(az_uplim))] = 90.


    # Plot twilight pointing constraints
    ax.scatter(np.deg2rad(sunrise_az), 3., marker='*', color='r', edgecolor='k', s=200, label=f'Sunrise location ({sunrise_UTC} UTC)')

    ax.fill_between(linspace_az, horizon_el, twi_constraints_el, color='red', alpha=0.3, 
                label=f'Sunrise limits ({az_lolim:0.0f} < az < {az_uplim:0.0f})')
    
    # Find timestamps of when object rises & sets
    start_vis = np.where(np.isfinite(obj_el))[0][0] # index where object rises/observing window starts
    end_vis = np.where(np.isfinite(obj_el))[0][-1] # index where object sets/observing window ends
    print(start_vis, end_vis)
    
    # Plot trace of object across sky
    plt.plot(np.deg2rad(obj_az), obj_el, label=obj_label)
    plt.scatter(np.deg2rad(obj_az[start_vis]), obj_el[start_vis], marker='o', fc=None, ec='C0')
    plt.scatter(np.deg2rad(obj_az[end_vis]), obj_el[end_vis], marker='o', color='C0')
    
    # Annotate rise, set times of object
    t1 = trace_times[start_vis].strftime('%H:%M')
    t2 = trace_times[end_vis].strftime('%H:%M')
    
    ax.annotate(t1, (np.deg2rad(nep_az[start_vis])*1.02, nep_el[start_vis]*1.02), xycoords='data')
    ax.annotate(t2, (np.deg2rad(nep_az[end_vis])*1.02, nep_el[end_vis]*1.02), xycoords='data')

    # Set legend
    datestr = sunrise_UTC.split(' ')[0]
    fig.legend(title=f'Ephemeris {datestr}', frameon=True, fancybox=True, loc='upper center',
               title_fontproperties={'weight':'semibold', 'size':'large'})
    
    # Display
    return fig




# def visperday(toi_num, year, site):
#     '''
#     toi_num :: 
#     year :: str, e.g. '2025'
#     site :: EarthLocation object
#     '''
    
#     print('New row:', toi_num)
    
#     # Read in toi_list, ephem
#     toi_list = pd.read_csv(f'{toipath}/exofop_tess(1).csv', skiprows=2, on_bad_lines='skip')
#     ephem = pd.read_csv(f'{ephempath}/twilight_ephemerides_keckII_{year}.csv')
    
#     # Manually find ra/dec of test TOI
#     loc = np.where(toi_list['TOI']==toi_num)[0][0] # Find first instance of TOI in table
#     toi_coord = SkyCoord(toi_list['RA'][loc], toi_list['Dec'][loc], unit=(u.hourangle, u.deg))

#     # Iterate through 365 days per year, check twilight visibility
#     vis_matrix_row = np.full(len(ephem), False)

#     for j in range(len(ephem)): # for each day in given year

#         # Grab
#         sunrise_dt = datetime.strptime(ephem['sunrise_UTC'][j], '%Y-%b-%d %H:%M')
#         is_up = is_target_up_morn_twi(sunrise_dt, toi_coord, site=site)

#         vis_matrix_row[j] = is_up
        
#     return vis_matrix_row



def get_sidereal_time(t, site):
    '''
    t :: datetime object
    site :: EarthLocation object specifying observing site
    
    Returns: sidereal time (degrees) at specified time and location
    '''
    
    tt = Time(t, scale='utc') # make given time into astropy.Time object
    stime = tt.sidereal_time(kind='apparent', longitude=site)
    
    return stime.degree


def write_starlist(df, tois=None, ra=None, diff=60., Hmag=False, filename='toi_list', savepath='~'):
    '''
    Write .txt starlist file with columns [TOI, RA, DEC].
    tois     :: list of TOI numbers (floats). If not None, write only these targets to starlist.
    ra       :: right ascension (deg). if not None, only write targets within +-XXX degrees of this ra
    diff     :: [ra-diff, ra+diff] range to return if ra is not None
    Hmag     :: if True, write .txt file with additional Hmag column
    filename :: specify custom filename if you'd like
    '''
    
    # If list of specific TOIs is given, limit to only these objects
    if tois is not None:
        df = df[df['TOI'].isin(tois)]
    
    
    # Make new df with stripped-down columns
    df_new = df[['TOI', 'RA', 'Dec']]
    
    
    # Re-format TOI column: XXX.X > TOIXXX.X
    new_TOI_col = ['TOI'+str(df_new['TOI'].values[i]) for i in range(len(df_new))]
    df_new['TOI'] = new_TOI_col
    
    
    # Add in Hmag if desired
    if Hmag:
        df_new['Hmag'] = df['Hmag'].values
        
        
    # If specified, throw out targets outside given RA range
    if ra is not None:
    
        target_ra = df['RA_deg'].values
        angle_diff = np.abs(target_ra - ra) % 360.

        ix = (angle_diff <= diff)
        df_new = df_new[ix]
        
    # Write .txt file
    np.savetxt(f'{savepath}/{filename}.txt', df_new.values, fmt='%s', delimiter=',')
    
    return df_new



