	#this will work seemlesly with geb.py (not counting importing)
	#plot histogram

plt.hist(h, 20000, cumulative=True)
plt.xlabel('Elevation')
plt.ylabel('Culmalative Frequency')
plt.title('Histogram of Global Elevation in meters')

hmin= h.min()
hmax= h.max()
plt.xlim(hmin, hmax)

plt.show()

	


