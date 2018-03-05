# Amateur-Satellite-Tracker
This repo contains the Python code to track amateur satellites across the Earth.

## Dependencies
There are a fair number of Python dependencies you will need which uses Python 2.7. Installed, you must have:

* Numpy/Scipy: basic and advanced mathematical functionality
* urllib2: Opens the websites needed to HTML scrape the satellite TLE’s
* BeautifulSoup: performs the HTML scraping
* Basemap: This actually plots the map of the world with the satellites across it
* Tkinter: Creates the GUI windows to ingest data from the user
* pytz: A package to get local time at your QTH (may already be installed with most distributions)
* Pyephem: used for validation of my own answers. Really nice library though to do all things astronomical

## Basic Operation

The steps to using the software are relatively simple. The user first runs the main executable at a terminal prompt:

`python main.py`

From there a GUI window pops up displaying checkboxes next to all the satellites, two entry fields for a start and end UTC time (in 24 hour format), and a Run/Quit button. Once the Run button is pressed, a map will pop up displaying the satellite track across the Earth for the chosen times. Currently, I would only recommend only plotting 1-2 satellites at a time since the map can become quite cluttered with tracks. Also, I would make the stop time only a few hours different from the start time since many of the satellites will have already executed a complete orbit in that timeframe.

## Under the Hood

First the software opens the URL of where the latest and greatest satellite TLEs are stored and scrapes them for all the user choices. Next, if the user chooses a UTC time before the most current TLE, the program defaults to using the epoch UTC from the TLE. As for the end date, I wouldn’t go any further than a few days out if you want the results to be fairly accurate.

Once the data has been gathered, internal routines convert the information in the TLE into the needed orbital or Keplerian elements to find the orbital path. Those Keplerian elements then undergo further processing to two different reference frames until finally the latitude and longitude can be extracted. I’d be happy to go into significant detail on what these algorithms are, but I’m pretty sure you probably don’t care! 🙂 Besides, you could easily check out the source code.

## To-Do

I do plan on updating this a bit, so be sure to pull down new updates from the repo often. But, I’m actually really looking for feedback from one or two people to tell me what functionality they would like to see. For now, however, what I have planned is:

* Display satellite footprints to better see when the satellite is up relative to the QTH
* Have the map run in real time thereby eliminating the entered start/stop time
* Compute satellite Doppler offsets
* Have a function to show which satellites will be up relative to your QTH within a given 48 hour period
* Validate results beyond a few sample cases involving the ISS

So there you have it! Hopefully if I get the above to-do incorporated into the code as well as comms with the rig (to update VFO based on Doppler, for instance), then I should be able to mimic the MacDoppler capability by only paying the price of time. I’d really like this to be an open-source substitute for that software. There’s still much to be done for this to be a really useful tool, but it’s coming along. Also please alert me to any bugs found. I know there are a ton probably. Happy tracking!
