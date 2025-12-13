import numpy as np
import pandas as pd
from tqdm import tqdm
from astroquery.jplhorizons import Horizons
import calendar
from datetime import datetime
from datetime import timedelta

'''
This script generates a lookup table .csv of twilight-relevant solar ephemerides for Keck II. It takes ~4 minutes to run. 
To generate a lookup table for a give year, change the value of the 'year' parameter on Line 16. If you want to generate
a lookup table for another observatory, change the 'location' input of the JPL Horizons query on Line 36.
'''

#################################################################################
year = 2026
savepath = '/Users/isabelkain/Desktop/Twilight_/twilight-observing-tool'

print( 'Querying JPL Horizons to build Keck II solar ephemeris lookup table.' )
print( 'Ephemeris year: ' + str(year) )
print( 'File location:  ' + savepath )

# If you're having problems (or need to re-query results), uncomment this line and run the script again:
# Horizons.clear_cache()
#################################################################################

def query_JPL_horizons(start, end):
    '''
    start: datetime object, year month day 03:00:00
    start: datetime object, year month day 18:00:00
    '''

    # Set epoch to query
    epoch =  {'start':start.strftime('%Y-%m-%d %H:%M:%S'), 'stop':end.strftime('%Y-%m-%d %H:%M:%S'), 'step':'1m'}

    
    # Send query
    obj = Horizons(id='sun', location='T17', epochs=epoch) # T17 = Keck II, T16 = Keck I, T15 = Gemini N
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
    sunrise_time = result['datetime_str'][-1] # datetime string, e.g. '2025-Jan-01 16:46'
    sunrise_az = result['AZ'][-1] # deg
    sunrise_el = result['EL'][-1] # deg
    sunrise_ra = result['RA_app'][-1] # deg
    sunrise_dec = result['DEC_app'][-1] # deg
    sunrise_sidereal = result['siderealtime'][-1] # hours

    
    # Return all params as single row
    return start.strftime('%Y-%m-%d'), sunset_time, sunset_az, sunset_el, sunset_ra, sunset_dec, \
            sunset_sidereal, sunrise_time, sunrise_az, sunrise_el, sunrise_ra, sunrise_dec, sunrise_sidereal


#################################################################################
#################################################################################


# Set range of dates to query

start_time = '03:00:00' # start search time UTC (evening, 5pm HT)
end_time = '18:00:00' # end search time UTC (morning, 8am HT)

numdays = 365 + calendar.isleap(year) # year must be int


# Make lists of datetime inputs

start_date = datetime.strptime(f'{year}-01-01 {start_time}', '%Y-%m-%d %H:%M:%S')
start_datetimes = [start_date + timedelta(days=x) for x in range(numdays)]

end_date = datetime.strptime(f'{year}-01-01 {end_time}', '%Y-%m-%d %H:%M:%S')
end_datetimes = [end_date + timedelta(days=x) for x in range(numdays)]


# Make pandas dataframe to catch results

columns = ['Date', 'sunset_UTC', 'sunset_az', 'sunset_el', 'sunset_ra', 'sunset_dec', 'sunset_sidereal', \
           'sunrise_UTC', 'sunrise_az', 'sunrise_el', 'sunrise_ra', 'sunrise_dec', 'sunrise_sidereal']

df = pd.DataFrame(columns=columns, index=np.arange(numdays))


# Iterate through all 365 days

assert len(start_datetimes)==len(end_datetimes)

for i in tqdm(range(len(start_datetimes))):
    
    start = start_datetimes[i]
    end = end_datetimes[i]
    
    results = query_JPL_horizons(start, end) # type(row)=tuple
    
    # Save new row to DataFrame
    newrow = dict(zip(columns, results))
    df.loc[i] = newrow
    
    
# Save lookup table

df.to_csv(f'{savepath}/twilight_ephemerides_keckII_{year}.csv', index=False)

