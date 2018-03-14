import numpy as np

# constants needed for amateur radio satellite tracking

GM = 3.986004418e14                                # mu factor (m^3/s^2)
J2 = 1.0827e-3                                     # oblate spheroid factor
Re = 6.378137e6                                    # Earth radius (m)
deg2rad = np.pi / 180.0
rad2deg = 180.0 / np.pi
twoPi = 2.0*np.pi
day2sec = 1.0 / (24.0 * 3600.0)                    # days to seconds
#num_time_pts = 5000                     
num_time_pts = 1
omega_earth = 0.2506844773746215 * (deg2rad / 60)  # convert (deg/min) -> (rad/s)
lat0 = 45.0*np.pi/180.0
earthEquatorialRadius = Re
earthPolarRadius = 6.356752e6
