import numpy as np
from CoordinateConversions import *
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import constants as c
import Tkinter as tk
import tkFont
from matplotlib import animation
from kep_to_state import *

fig = plt.figure()
ax = fig.add_subplot(111)

# miller projection
myMap = Basemap(projection = 'mill', llcrnrlat=-90, urcrnrlat=90, \
                      llcrnrlon=-180, urcrnrlon=180, resolution='c')
myMap.drawparallels(np.arange(-90, 90, 30), labels=[1, 0, 0, 0])
myMap.drawmeridians(np.arange(myMap.lonmin, myMap.lonmax+30, 60), labels=[0, 0, 0, 1])
myMap.bluemarble()
x, y = myMap([], [])
point = myMap.plot(x, y, 'r+', markersize=4, latlon=False)[0]
annotation = ax.annotate('', xy=(0, 0))
annotation.set_animated(True)


def animate(ii, selected_tle_dict):    
    latslons_dict = ConvertKepToStateVectors(selected_tle_dict)

    xList = []
    yList = []
    for key in selected_tle_dict:
        latslons = latslons_dict[key]
        lons = latslons[:,0]
        lats = latslons[:,1]
        x, y = myMap(lons, lats)
        xList.append(x)
        yList.append(y)

        ax.annotate(key, xy=(x, y), xytext=(x+4, y+4), color='yellow')
        #annotation = ax.annotate(key, xy=(x, y), xytext=(x+4, y+4), color='yellow')
        #annotation.set_animated(True)
        

    point.set_data(xList, yList)

    # shade the night areas, with alpha transparency so the
    # map shows through. Use current time in UTC.
    #date = datetime.utcnow()
    #CS = myMap.nightshade(date)
    #plt.title('Satellite Tracks for %s (UTC)' % date.strftime("%d %b %Y %H:%M:%S"))
    #plt.legend()
    return point, #annotation
    
    
    

	
#def runPredictionTool(start_time_list, end_time_list, checkbox_dict, tle_dict):
def runPredictionTool(checkbox_dict, tle_dict):

    selected_tle_dict = {}
    for key in checkbox_dict:
        if (checkbox_dict[key].get() == 1):
            selected_tle_dict[key] = tle_dict[key]
            
	# get lat/lons for all selected sats
	#latslons_dict = ConvertKepToStateVectors(selected_tle_dict, utc_start_time, utc_end_time)	
    ani = animation.FuncAnimation(fig, animate, fargs=(selected_tle_dict, ), interval=500, \
                                      blit=False, repeat=False)

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




### MAIN EXECUTABLE ###
if __name__ == "__main__":
    root = tk.Tk()
    SetupWindow(root)
    root.mainloop()
