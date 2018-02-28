import numpy as np
import urllib2
from bs4 import BeautifulSoup

'''
Decode 2-line elsets with the following key:
1 AAAAAU 00  0  0 BBBBB.BBBBBBBB  .CCCCCCCC  00000-0  00000-0 0  DDDZ
2 AAAAA EEE.EEEE FFF.FFFF GGGGGGG HHH.HHHH III.IIII JJ.JJJJJJJJKKKKKZ
KEY: A-CATALOGNUM B-EPOCHTIME C-DECAY D-ELSETNUM E-INCLINATION F-RAAN
G-ECCENTRICITY H-ARGPERIGEE I-MNANOM J-MNMOTION K-ORBITNUM Z-CHECKSUM
'''

def ParseTwoLineElementFile():
    # record website we are using to extract TLE
    quote_page = 'https://www.celestrak.com/NORAD/elements/amateur.txt'

    # open the website (use proxy if at MITLL)
    '''
    proxy = urllib2.ProxyHandler({'https': 'https://llproxy:8080'})
    opener = urllib2.build_opener(proxy)
    urllib2.install_opener(opener)
    page = urllib2.urlopen(quote_page)
    '''
    page = urllib2.urlopen(quote_page)
    
    # use BeautifulSoup to parse the data
    soup = BeautifulSoup(page, 'html.parser')
    soup = str(soup.get_text()).splitlines()
	
    counter = 0
    # create numpy array to cache results
    results = np.zeros(9, dtype=float)
    results_dict = {}
    for line in soup:
        split_line = line.split(" ")
        
        if (counter == 0):
            # this is an identifier for a new satellite
            # and len(split_line) == 2
            sat_name = split_line[0].strip('\n')
        elif (counter == 1):
            # remove empty strings from list
            split_line = filter(None, split_line)

            # satellite number
            sat_num = split_line[1]
            designator = split_line[2]

            epoch_info = split_line[3]
            epoch_year = epoch_info[0:2]

            # units of days
            epoch_remainder = epoch_info[2:]

            # ftdmm = first time derivative of mean motion divided by 2
            ftdmm = split_line[4]

            # stdmm = second time derivative of mean motion divided by 6
            stdmm = split_line[5]

            # bstar drag term
            bstar_drag = split_line[6]

            # cache results
            results[0] = epoch_year
            results[1] = epoch_remainder
        
        elif (counter == 2):
            # remove empty strings from list
            split_line = filter(None, split_line)

            # inclination (degrees)
            inclination = split_line[2]

            # right ascension of the ascending node (degrees)
            raan = split_line[3]
 
            # eccentricity (need to add decimal)
            ecc = '.'+split_line[4]
 
            # argument of perigee (degrees)
            arg_perigee = split_line[5]
 
            # mean anomaly (degrees)
            mean_anomaly = split_line[6]

            # mean motion (revolutions per day)
            mean_motion = split_line[7]

            # cache results
            results[2] = inclination
            results[3] = raan
            results[4] = ecc
            results[5] = arg_perigee
            results[6] = mean_anomaly
            results[7] = mean_motion
            results[8] = ftdmm

            # store in dictionary and reset array
            results_dict[sat_name] = results
            results = np.zeros(9, dtype=float)
            
        counter += 1
        counter = np.mod(counter, 3)

    return results_dict


#print ParseTwoLineElementFile()
