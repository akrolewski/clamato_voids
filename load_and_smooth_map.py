import numpy as np
import scipy.ndimage

def load_and_smooth_map(file,shape):
	map = np.fromfile(file)
	map = np.reshape(map,(np.shape(map)[0]/5,5))
	map = np.reshape(map,shape)

	'''First step in void finder.  Find all points with delta_F^{recon} > 0.224.
	Expand them into contiguous regions with <delta_F> = 0.167 (spherical underdensities).
	Find all regions regardless of whether they are overlapping'''

	kernel_size = 2.0
	pix_size = 0.5 #Each pixel is 0.5 h^{-1} Mpc

	map_smoothed = scipy.ndimage.filters.gaussian_filter(map,kernel_size/pix_size,mode='constant',cval=np.mean(map))
	# What smoothing scale to use?
	# How to handle map boundaries?
	# Did Casey smooth?

	xind = np.linspace(0,shape[0]/2,shape[0])
	yind = np.linspace(0,shape[1]/2,shape[1])
	zind = np.linspace(0,shape[2]/2,shape[2])
	allx,ally,allz = np.meshgrid(xind,yind,zind,indexing='ij')

	return map_smoothed, allx, ally, allz
