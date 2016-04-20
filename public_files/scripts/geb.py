	#Loads a netCDF file of 30 arc second gridded earth elevation data.

	# using "Python 3.4.2 |Anaconda 2.1.0 (64-bit)| (default, Oct 21 2014, 17:16:37) "
	# and "IPython 2.2.0"


import numpy as np
import scipy
import netCDF4
from netCDF4 import Dataset
import matplotlib
import matplotlib.pyplot as plt
from numpy import reshape
import timeit

tic= timeit.default_timer()	#system time call 1


	# location on disk, vaible link below, but as of Early Spring 2015 no direct ftp, 
	#requires loggin and clicking through several webpages to download NC file.

	# Data: http://www.gebco.net/data_and_products/gridded_bathymetry_data/documents/gebco_2014.pdf

geb= Dataset('/media/brandon/Active/MT_DEEP_SEA/Bathymetry/Gebco_2014/DATA/2D/GEBCO_2014_2D.nc')





ele = geb.variables['elevation'][:]	 # + is up - is down... well duh. 0 is mean tidal state.
lat = geb.variables['lat'][:]
lon = geb.variables['lon'][:]

geb.close() 	# close the file 

avg_elev= ((np.mean(ele))*-1)

	# To show script has run, an interesting factoid.
print()
print("The average elevation of the planet is", np.round(avg_elev,5), "meters below sea-level." )
print()
print("Calculated from Gebco 2014 data")


h= np.reshape(ele,-1)	 # array to vector -1 tells reshape to use all of the array values.

toc= timeit.default_timer()	 #system time call 2

print()
print(" script run time:", np.round(toc-tic,4), "seconds")	 # difference in, system (I think), time between the two calls


		# to plot histogram use the script "plot.py"


#print(numpy.mean(x))
#import mpl_toolkits
#from mpl_toolkits.basemap import Basemap





