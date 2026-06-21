import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
# from astropy.io import fits
from datetime import date, datetime, timedelta, timezone
# from astroquery.jplhorizons import Horizons
import astropy.units as u
import pandas as pd
# import xarray as xr
# from tqdm import tqdm
import calendar
from matplotlib.patches import Rectangle
from astropy.coordinates import AltAz, EarthLocation, SkyCoord, get_body, get_sun
from astropy.time import Time
# import multiprocessing

from utils import *

# %load_ext autoreload
# %autoreload 2

ephempath = '/Users/isabelkain/Desktop/Twilight_/twilight-observing-tool'
toipath = '/Users/isabelkain/Desktop/Twilight_/twilight-observing-tool'


def parse_datestring(datestr):
    '''
    Ingests string of unknown datetime format and checks it against a number of datetime formats:
    
    ** Year, month and day specified: **
    YYY-MM-DD > %Y-%m-%d
    YY-MM-DD > %y-%m-%d
    MM-DD-YYYY > %m-%d-%Y
    MM-DD-YY > %m-%d-%y
    mnth DD YYYY > %b %d %Y
    mnth DD YY > %b %d %y
    Month DD YYYY > %B %d %Y
    Month DD YY > %B %d %y
    
    ** Month and day only specified: **
    MM-DD > %m-%d
    mnth DD > %b %d
    month DD > %B %d
    
    ** User specifies current day: **
    Today
    Tonight
    Now
    '''
    
    # Write out datetime formats to check user input against
    list_of_dtformats_withyear = ['%Y-%m-%d', '%y-%m-%d', '%m-%d-%Y', '%m-%d-%y', '%b %d %Y', '%b %d %y', '%B %d %Y', '%B %d %y']
    list_of_dtformats_woutyear = ['%m-%d', '%b %d', '%B %d']
    list_of_dtformats_today = ['today', 'tonight', 'now']
    
    
    # Create datetime object for current day
    now = datetime.now()
    
    
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
    User can input an optional date string (YYY-MM-DD) to pull ephemerides
    from this date; otherwise, the code pulls ephemerides from today.
    '''
    
    # Parse date if given by user
    if datestr is not None:
    	date, year = parse_datestring(datestr)
    	datestr = date.strftime('%Y-%m-%d') # overwrite with standard format, just in case
    
    # Otherwise, use today's date
    else:
        date = datetime.now()
        datestr = date.strftime('%Y-%m-%d')
        year = now.strftime('%Y')
        
        
    # Read in solar ephemeris table from matching year
    ephem = pd.read_csv(f'{ephempath}/twilight_ephemerides_keckII_{year}.csv')
    
    
    # Identify row of ephemeris table that matches date
    try: 
        
        ix = ephem.index[ephem['Date'] == datestr].tolist()

        try: 
            assert len(ix)==1
            ix = ix[0]
        except: 
            print(f'Multiple rows for the same date. Something went wrong for the {datestr} ephemeris!')
            return 0

    except: 
        print(f'Solar ephemeris not available for {datestr}. Go make one.')
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
        
        
def twilight_pointing_limits(az, sun_az, az_sep=45.):
    '''
    Check if telescope az is >=45˚ from the sun. (az, sun_az units in degrees)
    Since telescope shutters are vertical, el angle does NOT matter.
    '''

    # Find difference between sun, telescope az angles
    if az > 180.:
        az = az - 360.

    if sun_az > 180.:
        sun_az = sun_az - 360.
    
    az_diff = np.abs(az-sun_az)
    print(az_diff)
                
    
    # Check if input az obeys twilight pointing constraints
    if az_diff < az_sep: 
        return False
    else: 
        return True


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


if __name__ == "__main__":

	###############################
	# Examples of running script
	###############################

	# For a same-day observation of Neptune:
	# > python generate_ss_azel_plot.py Neptune
	# > python generate_ss_azel_plot.py neptune tonight

	# For an observation of Titan on a specified night (if year is not specified, current year is assumed):
	# > python generate_ss_azel_plot.py Saturn jan 31
	# > python generate_ss_azel_plot.py Saturn January 31 2026
	# > python generate_ss_azel_plot.py saturn 1-31
	# > python generate_ss_azel_plot.py saturn 01-31
	# > python generate_ss_azel_plot.py saturn 2026-01-31


	# Take user input for Solar System target
	try:
		target_name = str(sys.argv[1])
	except:
		print('Specify which Solar System target you would like to generate an az-el pointing limits plot for.')
		print('Options: Neptune, Uranus, Titan, etc (case insensitive).')


	# Set observatory location
	keck = EarthLocation.of_site('Keck')


	# Set date of observation being planned. If no date specified, the current date is assumed.
	try:
		datestr = str(sys.argv[2]) 
        date, year = parse_datestring(datestr)
        datestr = date.strftime('%Y-%m-%d') # overwrite datestr format for convenience
        print('User-entered observation date:', datestr)
	except:
        date = datetime.now()
        year = date.strftime('%Y')
        datestr = date.strftime('%Y-%m-%d')
		print('No observation date specified, taking today\'s date:', datestr)


	# For given observation date, find UTC times, az angle of sunset and sunrise
	sunset_UTC, sunset_az, _, _, _, _, sunrise_UTC, sunrise_az, _, _, _, _, = read_solar_ephemeris(datestr)


	# Pull target's path across the sky from pre-generated ephemeris tables
	targ_az, targ_el, trace_times = grab_ss_trace(target_name, keck, sunset_UTC, sunrise_UTC)


	# Make figure
	fig = make_azel_plot(sunset_UTC, sunset_az, sunrise_UTC, sunrise_az, 
	                     targ_az, targ_el, trace_times, obj_label=target_name.capitalize())

	plt.show()


