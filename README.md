# Amateur-Satellite-Tracker
This repo contains the Python code to track amateur satellites across the Earth.

There are a fair number of Python dependencies you will need which uses Python 2.7. Installed, you must have:

Numpy/Scipy: basic and advanced mathematical functionality
urllib2: Opens the websites needed to HTML scrape the satellite TLE’s
BeautifulSoup: performs the HTML scraping
Basemap: This actually plots the map of the world with the satellites across it
Tkinter: Creates the GUI windows to ingest data from the user
pytz: A package to get local time at your QTH (may already be installed with most distributions)
Pyephem: used for validation of my own answers. Really nice library though to do all things astronomical
