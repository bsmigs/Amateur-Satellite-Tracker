import numpy as np
from tle_to_kep import *
from CoordinateConversions import *
import scipy
from scipy.optimize import newton
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
from datetime import datetime
import constants as c
import Tkinter as tk
import tkFont
from matplotlib import animation

fig = plt.figure()

# miller projection
map = Basemap(projection = 'mill', llcrnrlat=-90, urcrnrlat=90, \
                      llcrnrlon=-180, urcrnrlon=180, resolution='c')

# plot coastlines, draw label meridians and parallels.
map.drawparallels(np.arange(-90,90,30), labels=[1,0,0,0])
map.drawmeridians(np.arange(map.lonmin, map.lonmax+30,60), labels=[0,0,0,1])
map.bluemarble()



def animate(ii, selected_tle_dict):    
    latslons_dict = ConvertKepToStateVectors(selected_tle_dict)			
    for key in latslons_dict:
        latslons = latslons_dict[key]
        lons = latslons[:,0]
        lats = latslons[:,1]
        x, y = map(lons, lats)
        map.plot(x, y, 'ro', markersize=2, latlon=False)
        #map.plot(x[0], y[0], 'ks', markersize=8, latlon=False, label='Start')
        #map.plot(x[-1], y[-1], 'bs', markersize=8, latlon=False, label='End')

    # shade the night areas, with alpha transparency so the
    # map shows through. Use current time in UTC.
    date = datetime.utcnow()
    CS = map.nightshade(date)
    #plt.title('Satellite Tracks for %s (UTC)' % date.strftime("%d %b %Y %H:%M:%S"))
    #plt.legend()

    
    
    

	
#def runPredictionTool(start_time_list, end_time_list, checkbox_dict, tle_dict):
def runPredictionTool(checkbox_dict, tle_dict):
    '''
    utc_start_time = ''
    for elem in start_time_list:
        utc_start_time += elem.get()
        utc_start_time += ' '
	
    utc_end_time = ''
    for elem in end_time_list:
        utc_end_time += elem.get()
        utc_end_time += ' '
		
    utc_start_time = utc_start_time.rstrip()
    utc_end_time = utc_end_time.rstrip()
	'''
	
    selected_tle_dict = {}
    for key in checkbox_dict:
        if (checkbox_dict[key].get() == 1):
            selected_tle_dict[key] = tle_dict[key]

	# get lat/lons for all selected sats
	#latslons_dict = ConvertKepToStateVectors(selected_tle_dict, utc_start_time, utc_end_time)	
    ani = animation.FuncAnimation(fig, animate, fargs=(selected_tle_dict,), interval=1000, \
                                      blit=True, repeat=False)

    plt.show()

    



def SetupWindow(root):
	root.title("Amateur Radio Satellite Tracking") # set title for main window
	helv36 = tkFont.Font(size=12, weight=tkFont.BOLD)

    # parse the current TLE file
	tle_dict = ParseTwoLineElementFile()

    # number of desired columns
	num_cols = 4
	num_rows = len(tle_dict) / num_cols

    # define a frame for the checkboxes
	frame_1 = tk.Frame(root)
	frame_1.grid(row=0, column=1)

	# populate the checkboxes with satellite names
	row_counter = 0
	col_counter = 0
	checkbox_dict = {}
	for key in tle_dict:
		var = tk.IntVar()
		cb = tk.Checkbutton(frame_1, text=key, font=("Arial", 14), justify=tk.LEFT, variable=var)
		cb.grid(row=row_counter, column=col_counter, sticky=tk.W)
		row_counter += 1
	
		if (np.mod(row_counter, num_rows) == 0):
			col_counter += 1
			row_counter = 0
			
		checkbox_dict[key] = var

	# define a frame for the entered start/end times
	'''
	frame_2 = tk.Frame(root)
	frame_2.grid(row=1, column=1, columnspan=12)

	# set up labels
	start_year_label = tk.Label(frame_2, text="UTC Start Time: Year", font=helv36)
	start_year_label.grid(row=0, column=0, pady=2, sticky=tk.W)
	end_year_label = tk.Label(frame_2, text="UTC End Time: Year", font=helv36)
	end_year_label.grid(row=1, column=0, pady=2, sticky=tk.W)

	start_month_label = tk.Label(frame_2, text="Month", font=helv36)
	start_month_label.grid(row=0, column=2, pady=2, sticky=tk.W)
	end_month_label = tk.Label(frame_2, text="Month", font=helv36)
	end_month_label.grid(row=1, column=2, pady=2, sticky=tk.W)

	start_day_label = tk.Label(frame_2, text="Day", font=helv36)
	start_day_label.grid(row=0, column=4, pady=2, sticky=tk.W)
	end_day_label = tk.Label(frame_2, text="Day", font=helv36)
	end_day_label.grid(row=1, column=4, pady=2, sticky=tk.W)

	start_hour_label = tk.Label(frame_2, text="Hour", font=helv36)
	start_hour_label.grid(row=0, column=6, pady=2, sticky=tk.W)
	end_hour_label = tk.Label(frame_2, text="Hour", font=helv36)
	end_hour_label.grid(row=1, column=6, pady=2, sticky=tk.W)

	start_minute_label = tk.Label(frame_2, text="Minute", font=helv36)
	start_minute_label.grid(row=0, column=8, pady=2, sticky=tk.W)
	end_minute_label = tk.Label(frame_2, text="Minute", font=helv36)
	end_minute_label.grid(row=1, column=8, pady=2, sticky=tk.W)
	
	start_second_label = tk.Label(frame_2, text="Second", font=helv36)
	start_second_label.grid(row=0, column=10, pady=2, sticky=tk.W)
	end_second_label = tk.Label(frame_2, text="Second", font=helv36)
	end_second_label.grid(row=1, column=10, pady=2, sticky=tk.W)

	# set up entry fields
	start_year = tk.Entry(frame_2, width=5)
	start_year.grid(row=0, column=1, pady=2, sticky=tk.W)
	end_year = tk.Entry(frame_2, width=5)
	end_year.grid(row=1, column=1, pady=2, sticky=tk.W)

	start_month = tk.Entry(frame_2, width=5)
	start_month.grid(row=0, column=3, pady=2, sticky=tk.W)
	end_month = tk.Entry(frame_2, width=5)
	end_month.grid(row=1, column=3, pady=2, sticky=tk.W)

	start_day = tk.Entry(frame_2, width=5)
	start_day.grid(row=0, column=5, pady=2, sticky=tk.W)
	end_day = tk.Entry(frame_2, width=5)
	end_day.grid(row=1, column=5, pady=2, sticky=tk.W)

	start_hour = tk.Entry(frame_2, width=5)
	start_hour.grid(row=0, column=7, pady=2, sticky=tk.W)
	end_hour = tk.Entry(frame_2, width=5)
	end_hour.grid(row=1, column=7, pady=2, sticky=tk.W)
	
	start_minute = tk.Entry(frame_2, width=5)
	start_minute.grid(row=0, column=9, pady=2, sticky=tk.W)
	end_minute = tk.Entry(frame_2, width=5)
	end_minute.grid(row=1, column=9, pady=2, sticky=tk.W)

	start_second = tk.Entry(frame_2, width=5)
	start_second.grid(row=0, column=11, pady=2, sticky=tk.W)
	end_second = tk.Entry(frame_2, width=5)
	end_second.grid(row=1, column=11, pady=2, sticky=tk.W)

	start_time_list = [start_year, start_month, start_day, start_hour, start_minute, start_second]
	end_time_list = [end_year, end_month, end_day, end_hour, end_minute, end_second]
	'''

	# define a frame for the run/quit buttons
	frame_3 = tk.Frame(root)
	frame_3.grid(row=2, column=1, columnspan=12)

	# run button
	#rb = tk.Button(frame_3, text='Run Prediction', font=helv36, bg='green', command=(lambda a=start_time_list, b=end_time_list, c=checkbox_dict, d=tle_dict: runPredictionTool(a, b, c, d)))
	rb = tk.Button(frame_3, text='Run Prediction', font=helv36, bg='green', command=(lambda c=checkbox_dict, d=tle_dict: runPredictionTool(c, d)))
	rb.grid(row=2, column=1, padx=15, pady=15, sticky=tk.S)

	# quit button
	qb = tk.Button(frame_3, text='Quit Program', font=helv36, bg='red', command=root.quit)
	qb.grid(row=2, column=2, padx=15, pady=15, sticky=tk.S)
	
	
	
def ConstructSatelliteTracks(latslons_dict):
	# miller projection
	map = Basemap(projection = 'mill', llcrnrlat=-90,urcrnrlat=90,\
					  llcrnrlon=-180,urcrnrlon=180,resolution='c')

	# plot coastlines, draw label meridians and parallels.
	#map.drawcoastlines()
	map.drawparallels(np.arange(-90,90,30), labels=[1,0,0,0])
	map.drawmeridians(np.arange(map.lonmin, map.lonmax+30,60), labels=[0,0,0,1])

	# fill continents 'coral' (with zorder=0), color wet areas 'aqua'
	#map.drawmapboundary(fill_color='aqua')
	#map.fillcontinents(color='coral',lake_color='aqua')
	map.bluemarble()

	for key in latslons_dict:
		latslons = latslons_dict[key]
		lons = latslons[:,0]
		lats = latslons[:,1]

		x, y = map(lons, lats)
		n_elem = x.size
		map.plot(x, y, 'ro', markersize=2, latlon=False)
		#map.plot(x[0], y[0], 'ks', markersize=8, latlon=False, label='Start')
		#map.plot(x[-1], y[-1], 'bs', markersize=8, latlon=False, label='End')

	# shade the night areas, with alpha transparency so the
	# map shows through. Use current time in UTC.
	date = datetime.utcnow()
	CS = map.nightshade(date)
	plt.title('Satellite Tracks for %s (UTC)' % date.strftime("%d %b %Y %H:%M:%S"))
	plt.legend()
	plt.show()

	

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



### MAIN EXECUTABLE ###
if __name__ == "__main__":
    #utc_start_time = '2018 02 25 22 4 9'
    #utc_end_time = '2018 02 26 19 38 00'

    root = tk.Tk()
    SetupWindow(root)
    root.mainloop()
