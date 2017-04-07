import numpy as np
import scipy.ndimage
import h5py
import time

def rebin(a, *args):
	# rebins a 3D array to a new shape
    shape = a.shape
    lenShape = len(shape)
    factor = np.asarray(shape)/np.asarray(args)
    evList = ['a.reshape('] + \
             ['args[%d],factor[%d],'%(i,i) for i in range(lenShape)] + \
             [')'] + ['.mean(%d)'%(i+1) for i in range(lenShape)]
    return eval(''.join(evList))


newf = h5py.File("/global/cscratch1/sd/akrolew/Alex/4096-512bin_tot_density_rs.h5")
dm_density = newf["tot_density_rs"][:,:,:]


# Rebin the DM density
rebin_size = 128
dm_density_bin = rebin(dm_density, rebin_size, rebin_size, rebin_size)

# Smoothing kernel: 1 h^-1 Mpc
# smoothing kernel should be ~width of filaments, so that our definition and the "persistence"-based definition employed in Dubois match
kernel_size = 2.0

# Mode=wrap means this is consistent with PBC
#dm_density_smoothed = scipy.ndimage.filters.gaussian_filter(dm_density_bin,kernel_size*(float(rebin_size)/100.0),mode='wrap')
dm_density_smoothed = dm_density_bin

x = np.linspace(0,100,128)
y = np.linspace(0,100,128)
z = np.linspace(0,100,128)
allx,ally,allz = np.meshgrid(x,y,z,indexing='ij')

is_void = np.ones_like(allx)

pix_size = 100./127.

underdensities = np.loadtxt('underdensities_nyx.txt')
vradius = underdensities[:,0]
vcenter = underdensities[:,1:]

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

t0 = time.time()
for i in range(len(non_overlap_vcenter)):
	center = non_overlap_vcenter[i]
	if non_overlap_vradius[i] > 2:
		dist = np.sqrt((allx-pix_size*center[0])**2.+(ally-pix_size*center[1])**2.+(allz-pix_size*center[2])**2.) 
		is_void[dist < non_overlap_vradius[i]] = 0
		print i, time.time()-t0

is_void.tofile('void_indicator_no_overlaps.bin')