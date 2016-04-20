# this script is similar to "geb.py" the purpose being to read in the netcdf files. The second section load the data in and saves a plot of a subeset of the data.

#After they have been read in I have saved them as hd5 files. This is a type of hiearchial data format that is growing in poplualrity with those that analys large datsest  in many fields. In a  nutshecll they are efficecnet and compact files for which you can open and retreive data for analysis. For anything larger  fortran or C would be a better bet, or explore cloud/cluster computing through "Hadoop". Newer Netcdf files are similar/ built upon HD5 file format.

				#see : http://en.wikipedia.org/wiki/Hierarchical_Data_Format
				# &  : https://www.youtube.com/watch?v=wZEFoVUu8h0

import numpy as np
import scipy
import netCDF4
from netCDF4 import Dataset
import h5py
import timeit
import matplotlib
import matplotlib.pyplot as plt

tic= timeit.default_timer()


geb= Dataset('/media/brandon/Active/MT_DEEP_SEA/Bathymetry/Gebco_2014/DATA/2D/GEBCO_2014_2D.nc')
g2 = geb.variables['elevation'][:]
geb.close()

g1= np.roll(g2, 21600, axis=1) 		#DAM Important!!! so that gridded data is the the same layout as strm changed from -180->180 						##to, 0-> 360

#g1.reshape(1,933120000)

del geb 				# deleting variables that are no longer being used is a simple way to be effcient with memory


f= h5py.File("gebco.hdf5", 'w') 	# w for write and read
f.create_dataset("gebco", data=g1)
f.close()
del g1
print('geb done')

strm = Dataset('/media/brandon/Active/MT_DEEP_SEA/Bathymetry/strmbackground/v.11/topo30.grd','r')
s1 = strm.variables['z'][:]
strm.close()
del strm

#s1.reshape(1,933120000)

h= h5py.File("strm.hdf5", 'w')
h.create_dataset("strm", data=s1)
h.close()
del s1


import gc

gc.collect() 
print('strm done')
					# hd5 files created




#part 2 loading and plotting

g= h5py.File("gebco.hdf5", 'r')		# 'r' for read only
s= h5py.File("strm.hdf5", 'r')

st= s['strm']
ge =g['gebco']


							#Plotting:



a= st[:,0][0:]                 # [all rows, column one] [ starting at the first value: all values: step size equals default]
b= ge[:,0][0:]	




del st

del ge


s.close()
g.close()

gc.collect() 				


#plt.scatter(a, b)
#plt.ylabel('GEBCO 2014 - m', fontsize= 13)
#plt.xlabel('STRM 30 PlUS V.11 - m', fontsize= 13)
#plt.grid(True)
#plt.suptitle( 'Transect along Prime Meridian: Comparison between Global Elevation Sources:', fontsize =14)
#plt.title('Elevation in meters. 21600 30-Arc-Second cells', fontsize= 13, y=1.01)




#plt.savefig("Elevation Comparison N-S at 0 longitude ", format='png', dpi= 1200)


toc= timeit.default_timer()

#plt.show()

print()
print("script run time:", np.round(toc-tic,4), "seconds")

































