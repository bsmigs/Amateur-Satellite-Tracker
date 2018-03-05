import numpy as np
import constants as c

def RAANPrecession(a, e, i):
	eSq = e * e
	precession = np.divide( -1.5 * c.J2 * np.sqrt(c.GM) * (c.Re * c.Re) * np.cos(i), np.power(a, 3.5) * (1.0 - eSq) * (1.0 - eSq) )
	
	#print "RAAN precession = ", np.mean(precession), "rads/s"
	
	return precession

def ArgPerigeePrecession(a, e, i):
	eSq = e * e
	sini = np.sin(i)
	sin_i_sq = sini * sini
	precession = np.divide( 0.75 * c.J2 * np.sqrt(c.GM) * (c.Re * c.Re) * (5.0 * sin_i_sq - 1.0), np.power(a, 3.5) * (1.0 - eSq) * (1.0 - eSq) )
	
	#print "Arg of perigee precession = ", np.mean(precession), "rads/s"
	
	return precession


def ConvertKeplerToECI(a, e, i, Omega, w, nu, time_vec):
    # This converts Keplerian/orbital elements
    # to the ECI frame, (X, Y, Z)

    # a = semi-major axis
    # e = eccentricity
    # i = inclination
    # Omega = right ascension of ascending node
    # w = argument of perigee
    # nu = true anomaly
	
	# update RAAN and argument of perigee
	# based on precession values
	w_precession = ArgPerigeePrecession(a, e, i)
	w = w + (time_vec * (24. * 3600.)) *  w_precession
	
	Omega_precession = RAANPrecession(a, e, i)
	Omega = Omega + (time_vec * (24. * 3600.)) *  Omega_precession

    # pre-calculate trig functions / other vars
	sinnu = np.sin(nu)
	cosnu = np.cos(nu)
	sini = np.sin(i)
	cosi = np.cos(i)
	sin_w_plus_nu = np.sin(w + nu)
	cos_w_plus_nu = np.cos(w + nu)
	sinOmega = np.sin(Omega)
	cosOmega = np.cos(Omega)
	cosw = np.cos(w)
	sinw = np.sin(w)
	esq = e * e

    # get distance to the central body in PQW frame
	r = np.divide( a * (1.0 - esq), 1.0 + e * cosnu )
	x_PQW = r*cosnu
	y_PQW = r*sinnu

    # define rotation matrix (leave out col 3
    # since z-component in rotating frame = 0)
	R11 = cosw*cosOmega - sinw*cosi*sinOmega
	R12 = -(sinw*cosOmega + cosw*cosi*sinOmega)
	R21 = cosw*sinOmega + sinw*cosi*cosOmega
	R22 = -sinw*sinOmega + cosw*cosi*cosOmega
	R31 = sinw*sini
	R32 = cosw*sini
		
    # get Cartesian position coordinates in ECI frame
	X_eci = R11 * x_PQW + R12 * y_PQW
	Y_eci = R21 * x_PQW + R22 * y_PQW
	Z_eci = R31 * x_PQW + R32 * y_PQW

    # get Cartesian velocities in ECI frame
	coeff = np.divide( np.sqrt(c.GM * a), r )
	sinE = np.divide( sinnu * np.sqrt(1 - esq), 1 + e * cosnu)
	cosE = np.divide( e + cosnu, 1 + e * cosnu )
	local_vx = coeff * (-sinE)
	local_vy = coeff * np.sqrt(1.0 - esq) * cosE

    # get velocity vectors
	Xdot_eci = R11 * local_vx + R12 * local_vy
	Ydot_eci = R21 * local_vx + R22 * local_vy
	Zdot_eci = R31 * local_vx + R32 * local_vy

	return X_eci, Y_eci, Z_eci, Xdot_eci, Ydot_eci, Zdot_eci


def ConvertECIToECEF(X_eci, Y_eci, Z_eci, gmst):
    # X_eci, Y_eci, Z_eci = Cartesian coordinates
    # in ECI frame.
    #
    # gmst = Greenwich mean sidereal time
    
    X_ecef = X_eci * np.cos(gmst) + Y_eci * np.sin(gmst)
    Y_ecef = -X_eci * np.sin(gmst) + Y_eci * np.cos(gmst)
    Z_ecef = Z_eci

    return X_ecef, Y_ecef, Z_ecef


def ComputeGeodeticLat(phi, X_ecef, Y_ecef, Z_ecef, a, e):
    # Function to be used in an iterative solver
    # like Newton-Raphson
    #
    # a = semi-major axis
    # e = eccentricity
	p = np.sqrt(X_ecef*X_ecef + Y_ecef*Y_ecef)
	psq = p*p
	kappa = np.divide(p, Z_ecef) * np.tan(phi)
	kappasq = kappa * kappa
	Zsq = Z_ecef * Z_ecef
	esq = e * e
	
	return (kappa - 1.0 - np.divide( esq * a * kappa, np.sqrt(psq + (1.0 - esq) * Zsq * kappasq) ) )

	
def DComputeGeodeticLat(phi, X_ecef, Y_ecef, Z_ecef, a, e):
    # This is derivative of ComputeGeodeticLat function
    # w.r.t. phi. May be useful for Newton-Raphson
    # type solvers
    p = np.sqrt(X_ecef*X_ecef + Y_ecef*Y_ecef)
    psq = p*p
    Zsq = Z_ecef * Z_ecef
    kappa = np.divide(p, Z_ecef) * np.tan(phi)
    kappasq = kappa * kappa
    esq = e * e
    secphi = 1.0 / np.cos(phi)
    alpha = np.divide(p, Z_ecef) * secphi * secphi
	
    numer = psq * esq * a
    denom = np.power( psq + (1.0 - esq) * Zsq * kappasq, 1.5)
	
    deriv_formula = alpha * (1.0 - np.divide( numer, denom ) )

    return deriv_formula


def ComputeGeodeticLat2(X_ecef, Y_ecef, Z_ecef, a, e):
    asq = a*a
    esq = e*e
    b = a * np.sqrt(1.0 - esq)
    bsq = b*b
    p = np.sqrt(X_ecef*X_ecef + Y_ecef*Y_ecef)
    ep = np.sqrt(asq - bsq) / b
    theta = np.arctan2(a * Z_ecef, b * p)
    sintheta = np.sin(theta)
    costheta = np.cos(theta)
    phi = np.arctan2( Z_ecef + ep*ep*b*np.power(sintheta, 3.0), p - esq*a*np.power(costheta, 3.0) )
    
    return phi

	
	
	
def ComputeGeodeticLon(X_ecef, Y_ecef):
	lons = np.arctan2(Y_ecef, X_ecef)
	
	return lons


def ComputeGeodeticAlts(X_ecef, Y_ecef, Z_ecef, lats):
    # this formula obtained from Wiki page
    # on geographic conversion, computing
    # altitude from 
    r = np.sqrt(X_ecef * X_ecef + Y_ecef * Y_ecef + Z_ecef * Z_ecef)
    rSq = r * r
    sinLat = np.sin(lats) # make sure in radians
    cosLat = np.cos(lats) # make sure in radians
    sinLatSq = sinLat * sinLat
    cosLatSq = cosLat * cosLat

    # get constants from the earth
    a = c.earthEquatorialRadius
    b = c.earthPolarRadius
	# flattening factor of earth
	f = 1.0 - (b / a)
	alpha = f*f - 2.0*f
	
	# now compute altitudes for each lat
	alts = np.sqrt( rSq - ( cosLatSq * sinLatSq * alpha ) / (1.0 + alpha * sinLatSq) ) - \
		a * np.sqrt( 1.0 + alpha * sinLatSq )
	
	return alts
	
	
    
