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

'''First step in void finder.  Find all points with delta_F^{recon} > 0.224.
Expand them into contiguous regions with <delta_F> = 0.167 (spherical underdensities).
Find all regions regardless of whether they are overlapping'''

centers = np.array(np.where(dm_density_smoothed<0.15))

min_r = 1.5
max_r = 50.0
nr = 500

pix_size = 100./127.

x = np.linspace(0,100,128)
y = np.linspace(0,100,128)
z = np.linspace(0,100,128)
allx,ally,allz = np.meshgrid(x,y,z,indexing='ij')

all_radii = np.linspace(min_r,max_r,nr)
all_rho_enclosed = np.zeros_like(all_radii)

vcenter = []
vradius = []

print 'hi'
print np.shape(centers)[1]
t0 = time.time()
for j in range(np.shape(centers)[1]):
	center = centers[:,j]
	dist = np.sqrt((allx-pix_size*center[0])**2.+(ally-pix_size*center[1])**2.+(allz-pix_size*center[2])**2.) 
	totdeltaf_density = np.zeros_like(all_radii)
	totdeltaf = np.cumsum(np.ndarray.flatten(dm_density_smoothed[dist < max_r])[np.argsort(np.ndarray.flatten(dist[dist < max_r]))])
	print center,j,time.time()-t0
	dist_sorted = np.ndarray.flatten(dist[dist < max_r])[np.argsort(np.ndarray.flatten(dist[dist < max_r]))]
	totdeltaf_density = totdeltaf/(np.arange(len(totdeltaf))+1)
	vcenter.append(center)
	vradius.append(np.max(dist_sorted[totdeltaf_density < 0.3]))
    
underdensities_file = open('underdensities_nyx_2mpc_smoothing.txt','w')
for i in range(len(vradius)):
	underdensities_file.write('%10.5f %10.5f %10.5f %10.5f\n' % (vradius[i],vcenter[i][0],vcenter[i][1],vcenter[i][2]))
