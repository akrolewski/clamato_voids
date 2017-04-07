import numpy as np
import scipy.ndimage
from load_and_smooth_map import *

mapfile = '/Users/ALEX/Berkeley/IGM_Nyx_CLAMATO/CLAMATO_06_16/map.bin'
map_smoothed, allx, ally, allz = load_and_smooth_map(mapfile,(36,48,680))

'''First step in void finder.  Find all points with delta_F^{recon} > 0.224.
Expand them into contiguous regions with <delta_F> = 0.167 (spherical underdensities).
Find all regions regardless of whether they are overlapping'''

centers = np.array(np.where(map_smoothed>0.224))

min_r = 1.5
max_r = 50.0
nr = 500

pix_size = 0.5

all_radii = np.linspace(min_r,max_r,nr)
all_rho_enclosed = np.zeros_like(all_radii)

vcenter = []
vradius = []

print 'hi'
print np.shape(centers)[1]
for j in range(np.shape(centers)[1]):
	center = centers[:,j]
	dist = np.sqrt((allx-pix_size*center[0])**2.+(ally-pix_size*center[1])**2.+(allz-pix_size*center[2])**2.) 
	totdeltaf_density = np.zeros_like(all_radii)
	totdeltaf = np.cumsum(np.ndarray.flatten(map_smoothed[dist < max_r])[np.argsort(np.ndarray.flatten(dist[dist < max_r]))])
	print center
	dist_sorted = np.ndarray.flatten(dist[dist < max_r])[np.argsort(np.ndarray.flatten(dist[dist < max_r]))]
	totdeltaf_density = totdeltaf/(np.arange(len(totdeltaf))+1)
	vcenter.append(center)
	vradius.append(np.max(dist_sorted[totdeltaf_density > 0.167]))
    
underdensities_file = open('underdensities_2mpc_smoothing_constant_boundary.txt','w')
for i in range(len(vradius)):
	underdensities_file.write('%10.5f %10.5f %10.5f %10.5f\n' % (vradius[i],vcenter[i][0],vcenter[i][1],vcenter[i][2]))
