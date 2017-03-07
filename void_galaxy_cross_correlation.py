import matplotlib.gridspec as gridspec
gs = gridspec.GridSpec(10,9)
import matplotlib
from matplotlib import font_manager
import numpy as np
from astropy.io import ascii
from astropy.cosmology import FlatLambdaCDM
import matplotlib.pyplot as plt
from load_and_smooth_map import *
from scipy import spatial

'''Big problem here: how to deal with the survey boundaries for zCOSMOS/MOSDEF/3DHST when stacking the voids?'''

# Map parameters and cosmology
cosmo = FlatLambdaCDM(H0=70, Om0=0.31)
ra0 = 149.975
dec0 = 2.15
zmid = 2.35
deg_per_hMpc = 1./ cosmo.h / cosmo.comoving_distance(zmid).value * 180./np.pi
dist_to_center = cosmo.comoving_distance(zmid).value*cosmo.h

# Load the map
mapfile = '/Users/ALEX/Berkeley/IGM_Nyx_CLAMATO/CLAMATO_06_16/map.bin'
map_smoothed, allx, ally, allz = load_and_smooth_map(mapfile)

# Load the galaxies and convert coordinates to h^{-1} Mpc
galaxies = ascii.read('galxcorr_cl2016_v0.dat')
galaxies_x = (galaxies['ra'] - ra0)/deg_per_hMpc
galaxies_y = (galaxies['dec'] - dec0)/deg_per_hMpc
galaxies_z = cosmo.comoving_distance(galaxies['zspec']).value*cosmo.h

# Select certain galaxies
galaxies_x = galaxies_x[galaxies['source'] == '3DHST']
galaxies_y = galaxies_y[galaxies['source'] == '3DHST']
galaxies_z = galaxies_z[galaxies['source'] == '3DHST']

# Load the voids.  Convert from map-centered coordinates to observer-centered coordinates
voids = np.loadtxt('voids.txt')
non_overlap_vradius = voids[:,0]
non_overlap_vcenter = voids[:,1:]
non_overlap_vcenter[:,2] = non_overlap_vcenter[:,2] + dist_to_center-170.

galaxies_pos = np.array(zip(galaxies_x,galaxies_y,galaxies_z))
tree = spatial.cKDTree(galaxies_pos)

#lower_y_voids = real_voids[np.where(non_overlap_vcenter[:,1] < 20)]
lower_y_voids = non_overlap_vcenter[np.where(non_overlap_vcenter[:,0] > 6)]

allmatch = tree.query_ball_point(lower_y_voids,10)
dists = np.zeros_like(allmatch)
for i in range(len(allmatch)):
	dists[i] = np.sqrt(np.sum((galaxies_pos[allmatch[i]]-lower_y_voids[i])**2.,axis=1))
	
dists_flat = np.concatenate(dists)
counts, bedges = np.histogram(dists_flat,bins=np.linspace(0,10,10))
volume = (4.*np.pi/3.)*(bedges[1:]**3-bedges[:-1]**3)
print counts/volume