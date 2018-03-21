import numpy as np
from tle_to_kep import *
from TimeRoutines import *
from CoordinateConversions import *
from datetime import datetime


def ConvertKepToStateVectors(tle_dict):
#def ConvertKepToStateVectors(tle_dict, utc_start_time, utc_end_time):
	# parse the current TLE file
	#tle_dict = ParseTwoLineElementFile()

	# construct orbital/keplerian elements at requested time
	utc_start_time = datetime.utcnow()
	utc_start_time = utc_start_time.strftime('%Y %m %d %H %M %S')
	utc_end_time = utc_start_time

	kep_elem_dict, time_vec, epoch_year = ConvertTLEToKepElem(tle_dict, utc_start_time, utc_end_time)
	time_array = Nth_day_to_date(epoch_year, time_vec)

	# get the julian day
	jday = JdayInternal(time_array)
	#print "jday=",jday

	# get the greenwich mean sidereal time
	gmst = CalculateGMSTFromJD(jday, time_vec)	
	#print "GMST=",gmst*c.rad2deg

	latslons_dict = {}
	for key in kep_elem_dict:
		values = kep_elem_dict[key]
		
		a = values[:,0] # semi major axis 
		e = values[:,1] # eccentricity
		i = values[:,2] # inclination
		Omega = values[:,3] # Right ascension of ascending node (big omega)
		w = values[:,4] # argument of perigee (little omega)
		nu = values[:,5] # true anomaly
		epoch_days = values[:,8] # days from epoch
		
		# Convert orbital elements to ECI frame
		delta_time_vec = time_vec - epoch_days
		X_eci, Y_eci, Z_eci, Xdot_eci, Ydot_eci, Zdot_eci = ConvertKeplerToECI(a, e, i, Omega, w, nu, delta_time_vec)

		# convert ECI to ECEF
		X_ecef, Y_ecef, Z_ecef = ConvertECIToECEF(X_eci, Y_eci, Z_eci, gmst)
		
		# convert ECEF to geodetic (only keep lat/lon)
		lons = ComputeGeodeticLon(X_ecef, Y_ecef)
		#lats = np.array([ scipy.optimize.newton(ComputeGeodeticLat, c.lat0, fprime=DComputeGeodeticLat, args=(X_ecef[ii], Y_ecef[ii], Z_ecef[ii], a[ii], e[ii])) for ii in range(0, a.size) ])
		lats = ComputeGeodeticLat2(X_ecef, Y_ecef, Z_ecef, a, e)

		# convert lats/lons to degrees
		lats *= c.rad2deg
		lons *= c.rad2deg

		# store answer in dictionary
		n_rows = len(lats)
		results = np.zeros((n_rows, 2), dtype=float)
		results[:,0] = lons
		results[:,1] = lats
		latslons_dict[key] = results

	return latslons_dict
