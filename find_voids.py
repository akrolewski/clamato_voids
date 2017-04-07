import numpy as np
import scipy.ndimage
from load_and_smooth_map import *


mapfile = '/Users/ALEX/Berkeley/IGM_Nyx_CLAMATO/CLAMATO_06_16/map.bin'
map_smoothed, allx, ally, allz = load_and_smooth_map(mapfile,(36,48,680))

'''Second step in void finder.  Remove all overlapping voids.  Start with the largest
void in the map.  Remove all voids that overlap with this void.  Move through the list
until all sub-voids are removed.'''

ud = np.loadtxt('underdensities_2mpc_smoothing_constant_boundary.txt')
vradius = ud[:,0]
vcenter = ud[:,1:]

pix_size = 0.5

# Convert from pixels to h^{-1} Mpc
vcenter = pix_size*np.array(vcenter)
# Radii are already in h^{-1} Mpc
vradius = np.array(vradius)

new_vcenter = np.copy(vcenter)
new_vradius = np.copy(vradius)

non_overlap_vcenter = []
non_overlap_vradius = []

cnt = 0
while np.any(new_vradius):
	maxvoid = new_vcenter[np.argmax(new_vradius)]
	maxradius = new_vradius[np.argmax(new_vradius)]
	inds_to_remove = np.where(np.sqrt(np.sum((new_vcenter-maxvoid)**2.,axis=1)) < maxradius+new_vradius)[0]
	print maxvoid, maxradius, np.shape(new_vcenter)
	non_overlap_vcenter.append(maxvoid)
	non_overlap_vradius.append(maxradius)
	new_vcenter = np.delete(new_vcenter,inds_to_remove,axis=0)
	new_vradius = np.delete(new_vradius,inds_to_remove)
	cnt = cnt + 1

voids_file = open('voids_2mpc_smoothing_constant_boundary.txt','w')
for i in range(len(non_overlap_vcenter)):
	voids_file.write('%10.5f %10.5f %10.5f %10.5f\n' % (non_overlap_vradius[i],non_overlap_vcenter[i][0],non_overlap_vcenter[i][1],non_overlap_vcenter[i][2]))

voids_file.close()