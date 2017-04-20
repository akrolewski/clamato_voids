import numpy as np
import scipy.ndimage
import h5py
import time
from scipy import spatial
from sklearn.neighbors import KDTree
import matplotlib.pyplot as plt
import copy

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

# No smoothing of DM field
dm_density_smoothed = dm_density_bin

# Make some grids
x = np.linspace(0,100,128)
y = np.linspace(0,100,128)
z = np.linspace(0,100,128)
allx,ally,allz = np.meshgrid(x,y,z,indexing='ij')

# Simulation parameters
omegam = 0.3
omegalambda = 1 - omegam
h = 0.685
redshift = 2.4 # Redshift of sims
hz = 100*h*np.sqrt(omegam*(1+redshift)**3.+omegalambda)

boxsize = 100.0 # in Mpc/h
boxsize_mpc = 100.0/h
npix_sim = 512.

pix_size = 100./127.

# Load the underdensities
# These were selected in "max_underdensity_file_nyx_mpi.py" which applies the SO void finder
# to find all voids with minimum overdensity < 0.15 and average overdensity = 0.3 (the Stark'15
# definition of voids in redshift-space)
# Note that these voids may be overlapping
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

# Remove overlapping voids
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

is_void = np.ones_like(allx)
t0 = time.time()
# Define all points within non-overlapping voids larger than 2 Mpc as "void regions"
# We will cross-correlate these "void region" points with the halo catalog
for i in range(len(non_overlap_vcenter)):
	center = non_overlap_vcenter[i]
	if non_overlap_vradius[i]*boxsize_mpc/100. > 2:
		dist = np.sqrt((allx-center[0])**2.+(ally-center[1])**2.+(allz-center[2])**2.) 
		#is_void[np.argmin(dist)] = 0
		is_void[np.where(dist < non_overlap_vradius[i])] = 0
		#is_void[dist < np.min(dist)+0.0005] = 0
		print i, time.time()-t0, np.min(dist), center, np.where(dist < np.min(dist)+0.0005)
		
voids = np.transpose(np.where(is_void == 0))*boxsize_mpc/127.

# Read in the real-space halo catalog
halofile = np.loadtxt("/global/cscratch1/sd/akrolew/Alex/catalog_z24_iso138.txt",usecols=(0,1,2,5))
halopos = halofile[:,0:3] #in Mpc (not Mpc/h)
halomass = halofile[:,3] # in Msun (not Msun/h)
print 'hi'

# Move the halos to redshift space (we are computing the x-corr in redshift space)
f = h5py.File("/global/cscratch1/sd/akrolew/Alex/4096-512bin.h5")
tau_red = f["derived_fields/tau_red"][:,:,:]
vel_z = f["native_fields/velocity_z"][:,:,:]/10**5 #Convert to km/s
# The output velocities are peculiar proper velocities, see 
# eqn. 9 in Almgren et al. (2013) (and Zarija verified that all he did
# to make the hdf5 files was multiply the boxlib output from eqn. 9 by 10^5)
# therefore, we can just divide them by the Hubble constant to get deltachi
deltachi = vel_z*(1+redshift)/hz

halo_inds = np.round(halopos*((npix_sim-1.0)/boxsize_mpc)).astype(int)

deltachi_halos = deltachi[halo_inds[:,0],halo_inds[:,1],halo_inds[:,2]]

halos_redshift_space = copy.copy(halopos)
halos_redshift_space[:,2] = halos_redshift_space[:,2] + deltachi_halos
# Ensure that all halos lie within the box volume
halos_redshift_space[halos_redshift_space > boxsize_mpc] = halos_redshift_space[halos_redshift_space > boxsize_mpc ]-boxsize_mpc
halos_redshift_space[halos_redshift_space < 0] = halos_redshift_space[halos_redshift_space < 0]+boxsize_mpc
halos_redshift_space.tofile('halos_redshift_space.bin')

# Make a KDTree for the halos, and query it to find the cross-correlation of voids around the halos
print 'before tree'
tree = KDTree(halos_redshift_space)
print 'made tree'

#n = 10000
t0 = time.time()
neighbors = tree.query_radius(voids,20,return_distance=True)
print time.time()-t0

# Repeat the same process for a catalog of uniformly distributed randoms
randoms = np.random.uniform(low=0,high=boxsize_mpc,size=(10**6,3))
tree_randoms = KDTree(randoms)
print time.time()-t0
tree_neighbors = tree_randoms.query_radius(voids,20,return_distance=True)
print time.time()-t0
np.save('/global/cscratch1/sd/akrolew/trees_gt2.npy',[neighbors,randoms,tree_neighbors])

# Histogram the galaxy-random pairs and the void-galaxy pairs, and compute the cross-correlation function
# This estimator is xi = VG/VR - 1
gal_dists = np.concatenate(neighbors[1])
random_dists = np.concatenate(tree_neighbors[1])
gal_dist_hist = np.histogram(gal_dists,bins=100)
counts = gal_dist_hist[0]
bedges = gal_dist_hist[1]
bins = 0.5*(bedges[:-1]+bedges[1:])

random_dist_hist = np.histogram(random_dists,bins=100)
# Should actually match up bins here.....
exp = random_dist_hist[0]

halos_redshift_space = np.fromfile('halos_redshift_space.bin').reshape((1209148, 3))
plt.ion()
plt.figure()
plt.plot(bins,(10.**6/np.shape(halos_redshift_space)[0])*counts.astype('float')/exp.astype('float')-1.)
plt.plot(x,np.zeros_like(x),color='r')
plt.xlabel('R (Mpc)')
plt.ylabel('xi')
plt.savefig('void_pixel_xcorr_gt2_mpc.pdf')