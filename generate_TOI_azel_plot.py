import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import sys, os
from datetime import date, datetime, timedelta, timezone
import astropy.units as u
import pandas as pd
import calendar
from matplotlib.patches import Rectangle
from astropy.coordinates import AltAz, EarthLocation, SkyCoord, get_body, get_sun
from astropy.time import Time

from generate_ss_azel_plot import *


def grab_TOI_trace(target_name, site, sunset_UTC, sunrise_UTC):
	'''
	target_name :: string name of TOI target
	site        :: EarthLocation object
	sunset_UTC  :: datetime object
	sunrise_UTC :: datetime object
	'''

	# Read in (ra, dec) of one TOI host
	toi_coord = SkyCoord.from_name(target_name)


	# Generate list of times throughout night to check where SS object is in the sky
	sample_times = Time( pd.date_range(start=sunset_UTC, end=sunrise_UTC, periods=30).to_numpy() )
	targ_az = np.zeros(len(sample_times))
	targ_el = np.zeros(len(sample_times)) 


	# Grab object az,el location at each timestamp
	for i, time in enumerate(sample_times):

		toi_coord_azel = toi_coord.transform_to(AltAz(obstime=Time(datetime.strptime(sunrise_UTC, '%Y-%b-%d %H:%M')), location=keck))
		targ_az[i] = toi_coord_azel.az.value
		targ_el[i] = toi_coord_azel.alt.value

		# If object not effectively visible with Keck II, mask values with NaNs
		if not keckII_pointing_limits(targ_az[i], targ_el[i]):
			targ_el[i] = np.nan

	return targ_az, targ_el, sample_times


if __name__ == "__main__":


	# For now, hardcode TOI target name -- later, replace with user input FIXME
	target_name = 'TOI 103.01'

	# Set observatory location
	keck = EarthLocation.of_site('Keck')

	# Read in solar ephemeris for today's date FIXME later do user selection of date
	sunset_UTC, sunset_az, _, _, _, _, sunrise_UTC, sunrise_az, _, _, _, _,  = read_solar_ephemeris()

	# FIXME need a check about whether this TOI is remotely visible

	# Pull target's path across the sky from pre-generated ephemeris tables
	targ_az, targ_el, trace_times = grab_TOI_trace(target_name, keck, sunset_UTC, sunrise_UTC)

	# print('TARGET AZ:', targ_az)
	# print('TARGET AZ % NAN:', 100*np.sum(~np.isfinite(targ_az))/np.size(targ_az), '%')
	# print('TARGET EL:', targ_el)
	# print('TARGET EL % NAN:', 100*np.sum(~np.isfinite(targ_el))/np.size(targ_el), '%')

	frac_vis = np.sum(np.isfinite(targ_el))/np.size(targ_el)

	if frac_vis < 0.75:
		raise AssertionError(f'Your chosen target is visible for <75% of the twilight window ({100*frac_vis:0.1f}% to be exact). Pick a better one!')


	# Make figure
	fig = make_azel_plot(sunset_UTC, sunset_az, sunrise_UTC, sunrise_az, 
						 targ_az, targ_el, trace_times, obj_label=target_name.upper())

	plt.show()

