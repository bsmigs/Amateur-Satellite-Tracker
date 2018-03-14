import datetime
import numpy as np
import pytz
import sys
import constants as c


def ConvertLocalTimeToUTC(local_time, used_format='%Y %m %d %H %M %S'):
    local_time = datetime.datetime.strptime(local_time, used_format)
    local = pytz.timezone("America/New_York")
    local_dt = local.localize(local_time, is_dst=False)
    utc_dt = local_dt.astimezone(pytz.utc)
    utc = utc_dt.strftime(used_format)
    
    return utc

def GenerateTimeVec(utc_start_time, utc_end_time, tle_epoch_year, tle_epoch_days):
    # start_time = a string of 'year month day hour min sec'
    # end_time = a string of 'year month day hour min sec'
    # tle_epoch_year = 2 digit year
    # tle_epoch_days = Nth day of the year (with fractional component)
    
    if (tle_epoch_year < 57): # Sputnik launched in 1957
        tle_epoch_year += 2000
    else:
        tle_epoch_year += 1900

    # convert string year portion to float
    future_start_year = float(utc_start_time[0:4])
    future_end_year = float(utc_end_time[0:4])
    if (future_start_year > future_end_year):
        future_start_year = tle_epoch_year
        future_end_year = tle_epoch_year
        print "Forcing entered start and end year to be same as TLE epoch year"
    
    if (future_start_year < tle_epoch_year):
        future_start_year = tle_epoch_year
        print "Forcing entered start year to be same as TLE epoch year"

    if (future_end_year < tle_epoch_year):
        future_end_year = tle_epoch_year
        print "Max future end year forced to be same as TLE epoch year"

    # convert future date to nth date of year
    future_start_days = Date_to_nth_day(utc_start_time)
    future_end_days = Date_to_nth_day(utc_end_time)
    
    if (future_start_days > future_end_days):
        print "Cannot choose future start time to be greater than future end time"
        sys.exit()

    if (future_start_days < tle_epoch_days):
        future_start_time = tle_epoch_days
        print "Entered future start day of year set equal to epoch from TLE: " + str(tle_epoch_days)

    # if the future time is actually in the future
    # then make a vector of times that we can evaluate
    time_vec = np.linspace(future_start_days, future_end_days, num=c.num_time_pts, endpoint=True)

	# this is in days
    return time_vec, tle_epoch_year


#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
def Date_to_nth_day(date, used_format='%Y %m %d %H %M %S'):
    date = datetime.datetime.strptime(date, used_format)
    new_year_day = datetime.datetime(year=date.year, month=1, day=1, hour=0, minute=0, second=0)
    delta = (date - new_year_day)
    num_days = ( delta.days + 1 ) + ( delta.seconds / (24.0 * 3600.0) )
    #print "num_days=",num_days
    
    return num_days

#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
def Nth_day_to_date(year, ndays):
	# it's encouraged that ndays is a float
	# so that we can encode more info about
	# hours, minutes, and seconds
	year_array_len = np.size(year)
	days_array_len = np.size(ndays)
	
	results = np.zeros((days_array_len, 6), dtype=int)

	if (year_array_len == 1 and days_array_len > 1):
		# it's assumed that the year applies
		# to all the days
		year = year*np.ones((days_array_len,), dtype=int)
		
		for ii in np.arange(0, days_array_len):
			this_date = datetime.datetime(year[ii], 1, 1, 0, 0, 0) + datetime.timedelta(ndays[ii] - 1.0)
			this_date = this_date.strftime('%Y %m %d %H %M %S')    
			this_date = np.fromstring(this_date, dtype=int, sep=' ')
			results[ii,:] = this_date
			
	elif (year_array_len == 1 and days_array_len == 1):
		#year = np.array([year], dtype=int)
		#ndays = np.array([ndays], dtype=float)

		# year is just a number, not an array. Not sure
		# if this is good coding practice
		this_date = datetime.datetime(year, 1, 1, 0, 0, 0) + datetime.timedelta(ndays[0] - 1.0)
		this_date = this_date.strftime('%Y %m %d %H %M %S')    
		this_date = np.fromstring(this_date, dtype=int, sep=' ')
		results[0,:] = np.array(this_date, dtype=int)
	
		
	return results

#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
'''
/* -----------------------------------------------------------------------------
 *
 *                           procedure days2mdhms
 *
 *  this procedure converts the day of the year, days, to the equivalent month
 *    day, hour, minute and second.
 *
 *  algorithm     : set up array for the number of days per month
 *                  find leap year - use 1900 because 2000 is a leap year
 *                  loop through a temp value while the value is < the days
 *                  perform int conversions to the correct day and month
 *                  convert remainder into h m s using type conversions
 *
 *  author        : david vallado                  719-573-2600    1 mar 2001
 *
 *  inputs          description                    range / units
 *    year        - year                           1900 .. 2100
 *    days        - julian day of the year         0.0  .. 366.0
 *
 *  outputs       :
 *    mon         - month                          1 .. 12
 *    day         - day                            1 .. 28,29,30,31
 *    hr          - hour                           0 .. 23
 *    min         - minute                         0 .. 59
 *    sec         - second                         0.0 .. 59.999
 *
 *  locals        :
 *    dayofyr     - day of year
 *    temp        - temporary extended values
 *    inttemp     - temporary int value
 *    i           - index
 *    lmonth[12]  - int array containing the number of days per month
 *
 *  coupling      :
 *    none.
 * --------------------------------------------------------------------------- */

def Days2ymdhms(year, days):
    lmonth = [31, (year % 4) == 0 ? 29 : 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31];
    dayofyr = np.floor(days);

    #  ----------------- find month and day of month ----------------
    i = 1;
    inttemp = 0;
    while ((dayofyr > (inttemp + lmonth[i - 1])) && i < 12):
        inttemp += lmonth[i - 1];
        i += 1;

    mon = i;
    day = dayofyr - inttemp;

    #  ----------------- find hours minutes and seconds -------------
    temp = (days - dayofyr) * 24.0;
    hr = np.floor(temp);
    temp = (temp - hr) * 60.0;
    minute = np.floor(temp);
    sec = (temp - minute) * 60.0;

    return np.array([year, mon, day, hr, minute, sec])
'''

#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
'''
/* -----------------------------------------------------------------------------
 *
 *                           procedure jday
 *
 *  this procedure finds the julian date given the year, month, day, and time.
 *    the julian date is defined by each elapsed day since noon, jan 1, 4713 bc.
 *
 *  algorithm     : calculate the answer in one step for efficiency
 *
 *  author        : david vallado                  719-573-2600    1 mar 2001
 *
 *  inputs          description                    range / units
 *    year        - year                           1900 .. 2100
 *    mon         - month                          1 .. 12
 *    day         - day                            1 .. 28,29,30,31
 *    hr          - universal time hour            0 .. 23
 *    min         - universal time min             0 .. 59
 *    sec         - universal time sec             0.0 .. 59.999
 *
 *  outputs       :
 *    jd          - julian date                    days from 4713 bc
 *
 *  locals        :
 *    none.
 *
 *  coupling      :
 *    none.
 *
 *  references    :
 *    vallado       2007, 189, alg 14, ex 3-14
 *
 * --------------------------------------------------------------------------- */
'''
def JdayInternal(ymdhms):

	year = ymdhms[:,0]
	mon = ymdhms[:,1]
	day = ymdhms[:,2]
	hr = ymdhms[:,3]
	min = ymdhms[:,4]
	sec = ymdhms[:,5]

	jday = (367.0 * year) - np.floor(7.0 * (year + np.floor((mon + 9.0) / 12.0)) * 0.25) + \
      np.floor(275.0 * mon / 9.0) + day + 1721013.5 # ut in days
		
	jdfrac = (sec + min*60.0 + hr*3600.0) / 86400.0		
	jday += jdfrac
	
	return jday

#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
'''
/* -----------------------------------------------------------------------------
 *
 *                           function gstime
 *
 *  this function finds the greenwich sidereal time.
 *
 *  author        : david vallado                  719-573-2600    1 mar 2001
 *
 *  inputs          description                    range / units
 *    jdut1       - julian date in ut1             days from 4713 bc
 *
 *  outputs       :
 *    gstime      - greenwich sidereal time        0 to 2pi rad
 *
 *  locals        :
 *    temp        - temporary variable for doubles   rad
 *    tut1        - julian centuries from the
 *                  jan 1, 2000 12 h epoch (ut1)
 *
 *  coupling      :
 *    none
 *
 *  references    :
 *    vallado       2004, 191, eq 3-45
 * --------------------------------------------------------------------------- */
'''
def CalculateGMSTFromJD(jdut1, time_vec):
	
	gmst = np.zeros((jdut1.size,), dtype=float)
	for ii in np.arange(0, jdut1.size):
	
		# need to figure out previous midnight
		# call it JD0
		JDmin = np.floor(jdut1[ii]) - 0.5
		JDmax = np.floor(jdut1[ii]) + 0.5
		if (jdut1[ii] > JDmin):
			JD0 = JDmin

		if (jdut1[ii] > JDmax):
			JD0 = JDmax
			
		#JD0 = JdayInternal(Nth_day_to_date(2018, np.floor(time_vec[ii])))

		# get GMST at previous midnight
		T = (JD0 - 2451545.0) / 36525.0
		#gmst = (-6.2e-6 * T * T * T) + (0.093104 * T * T) + (876600.0 * 3600 + 8640184.812866) * T + 67310.54841 # sec
		gmst00 = (-6.2e-6 * T * T * T) + (0.093104 * T * T) + (8640184.812866 * T) + 24110.548416 # sec
		gmst00 *= (360.0 / 86400.0) * c.deg2rad # 360/86400 = 1/240, to deg, to rad 
	
		# take GMST from previous midnight results
		# and add the angular displacement up to current UT vals
		# to do this, get fraction part of time vec since this
		# gives us the additional (hours,min,sec) the Earth
		# has rotated since midnight of the same day
		fractional_part = ( time_vec[ii] - np.floor(time_vec[ii]) ) * (24.0 * 3600.0) # (convert to sec)
		gmst[ii] = np.mod( gmst00 + c.omega_earth * fractional_part, c.twoPi )
    

	'''
    # Another way. Agrees REALLY well with above
    gmst = np.zeros((jdut1.size,), dtype=float)
    for ii in np.arange(0, jdut1.size):
        JDmin = np.floor(jdut1[ii]) - 0.5
        JDmax = np.floor(jdut1[ii]) + 0.5
        if (jdut1[ii] > JDmin):
            JD0 = JDmin

        if (jdut1[ii] > JDmax):
            JD0 = JDmax
        
        H = (jdut1[ii] - JD0)*24.       #Time in hours past previous midnight
        D = jdut1[ii] - 2451545.0       #Compute the number of days since J2000
        D0 = JD0 - 2451545.0            #Compute the number of days since J2000
        T = D / 36525.                  #Compute the number of centuries since J2000
        #Calculate GMST in hours (0h to 24h) ... then convert to degrees
        gmst[ii] = np.mod(6.697374558 + 0.06570982441908*D0  + 1.00273790935*H + 0.000026*(T*T), 24) * 15.0 * c.deg2rad
	'''
  
	return gmst
