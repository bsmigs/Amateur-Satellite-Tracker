import numpy as np
import scipy
from scipy.optimize import newton
from keplerian_parser import *
from TimeRoutines import *
import datetime
import sys
import constants as c

def KeplerEquation(eccentric_anomalies, mean_anomalies, ecc):
	return ( mean_anomalies - (eccentric_anomalies - ecc * np.sin(eccentric_anomalies)) )
	
def DKeplerEquation(eccentric_anomalies, mean_anomalies, ecc):
	return (-1.0 + ecc * np.cos(eccentric_anomalies))
	
def GetTrueAnomaly(eccentric_anomalies, ecc):	
    sinnu = np.sqrt(1.0 - ecc*ecc) * np.sin(eccentric_anomalies) # leave out denominator since it cancels out anyway
    cosnu = np.cos(eccentric_anomalies) - ecc                    # leave out denominator since it cancels out anyway
    nu = np.arctan2(sinnu, cosnu)
    
    return nu
	

def ConvertTLEToKepElem(tle_dict, utc_start_time, utc_end_time):    
	for key in tle_dict:
		tle_elem = tle_dict[key]

		epoch_year = int(tle_elem[0])
		epoch_days = tle_elem[1]                          # units of days
		inclination = tle_elem[2]*c.deg2rad               # radians
		raan = tle_elem[3]*c.deg2rad                      # radians
		ecc = tle_elem[4]                                 # dimensionless
		arg_perigee = tle_elem[5]*c.deg2rad               # radians
		mean_anomaly = tle_elem[6]*c.deg2rad              # radians
		mean_motion = tle_elem[7]*c.twoPi*c.day2sec       # radians per s
		ftdmm = tle_elem[8]*c.twoPi*c.day2sec*c.day2sec   # radians per s^2
		
        # semi-major axis
		#a = np.power(c.GM / np.power(mean_motion, 2.0), 1.0/3.0)

        # if the future time is actually in the future
        # then make a vector of times that we can evaluate (units = days)
		time_vec, epoch_year = GenerateTimeVec(utc_start_time, utc_end_time, epoch_year, epoch_days)
		
		#print "Satellite", key, "has date from TLE = ",Nth_day_to_date(epoch_year, epoch_days)
		
        # create vector of mean anomalies based on time vec
		delta_time_vec = ( (time_vec - epoch_days) * (24.0 * 3600.0) )
		current_mean_motions = (mean_motion + 0.5 * ftdmm * delta_time_vec)
		mean_anomalies = mean_anomaly + current_mean_motions * delta_time_vec
		mean_anomalies = np.mod(mean_anomalies, c.twoPi)
		
		# get current semi-major axes based on current times
		a = np.power( np.divide(c.GM, np.power(current_mean_motions, 2.0)), 1.0/3.0)

        # now figure out the eccentric anomaly
		# by solving Kepler's equation
		# M = E - e sin(E)
		initial_guess = mean_anomaly
		eccentric_anomalies = [ scipy.optimize.newton(KeplerEquation, initial_guess, fprime=DKeplerEquation, args=(x0, ecc)) for x0 in mean_anomalies]
		eccentric_anomalies = np.array(eccentric_anomalies)
        
		# derive true anomaly
		true_anomalies = GetTrueAnomaly(eccentric_anomalies, ecc)

		n_rows = true_anomalies.size        
		tmp = np.zeros((n_rows, 9), dtype=float)
		tmp[:,0] = a
		tmp[:,1] = ecc
		tmp[:,2] = inclination
		tmp[:,3] = raan
		tmp[:,4] = arg_perigee
		tmp[:,5] = true_anomalies
		tmp[:,6] = eccentric_anomalies
		tmp[:,7] = epoch_year
		tmp[:,8] = epoch_days

		tle_dict[key] = tmp
        
	return tle_dict, time_vec, epoch_year
	

'''
# parse the keplerian file
tle_dict = ParseTwoLineElementFile()

# create the dictionary of orbital elements
time = '2018 02 15 06 07 00'
tle_dict = ConvertTLEToKepElem(tle_dict, time)
print tle_dict["AO-85"]

# get the julian day
time = np.fromstring(time, dtype=int, sep=' ')
jday = JdayInternal(time)
print jday

# get the greenwich mean sidereal time
gmst = CalculateGMSTFromJD(jday)
print gmst
'''
