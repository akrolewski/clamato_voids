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
# Definitely will want to update this for the new maps

# Load the galaxies and convert coordinates to h^{-1} Mpc
galaxies = ascii.read('galxcorr_cl2016_v0.dat')
# I forget if this file includes VUDS or not...
galaxies_x = np.cos(0.5*(dec0+dec1)*np.pi/180.)*(galaxies['ra'] - ra0)/deg_per_hMpc
galaxies_y = (galaxies['dec'] - dec0)/deg_per_hMpc
galaxies_z = cosmo.comoving_distance(galaxies['zspec']).value*cosmo.h

# Load the voids.  Convert from map-centered coordinates to observer-centered coordinates
mapfile = '/Users/ALEX/Berkeley/IGM_Nyx_CLAMATO/CLAMATO_06_16/map.bin'
map_smoothed, allx, ally, allz = load_and_smooth_map(mapfile,(36,48,680))

'''First step in void finder.  Find all points with delta_F^{recon} > 0.224.
Expand them into contiguous regions with <delta_F> = 0.167 (spherical underdensities).
Find all regions regardless of whether they are overlapping'''

void_pixels = np.array(np.where(map_smoothed>0.224))
