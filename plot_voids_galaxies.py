import matplotlib.gridspec as gridspec
gs = gridspec.GridSpec(10,9)
import matplotlib
from matplotlib import font_manager
import numpy as np
from astropy.io import ascii
from astropy.cosmology import FlatLambdaCDM
import matplotlib.pyplot as plt
from load_and_smooth_map import *

'''Plot the voids and galaxies on top of slices of the map.  Each slice is separated by 2 h^{-1} Mpc.
Use integer division to find galaxies in each slice: select all galaxies between a given RA slice and the next one.'''

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
#galaxies_x = galaxies_x[galaxies['source'] == '3DHST']
#galaxies_y = galaxies_y[galaxies['source'] == '3DHST']
#galaxies_z = galaxies_z[galaxies['source'] == '3DHST']

# Load the voids.  Convert from map-centered coordinates to observer-centered coordinates
voids = np.loadtxt('voids.txt')
non_overlap_vradius = voids[:,0]
non_overlap_vcenter = voids[:,1:]
non_overlap_vcenter[:,2] = non_overlap_vcenter[:,2] + dist_to_center-170.

ticks_font = font_manager.FontProperties(family='Helvetica', style='normal',
		size=10, weight='normal', stretch='normal')

matplotlib.rcParams['text.usetex'] = True
matplotlib.rcParams['text.latex.preamble']=[r"\usepackage{amsmath}"]
matplotlib.rcParams['font.family'] = 'serif'
matplotlib.rcParams['font.serif'] = 'cm'

fig = plt.figure()
step = 4
for i in range(9):
	ax = plt.subplot(gs[9-i,:])
	ax.contourf(allz[step*i,:,:]+dist_to_center-170.,ally[step*i,:,:],map_smoothed[step*i,:,:],cmap=plt.cm.jet_r,levels=np.linspace(np.min(map_smoothed),np.max(map_smoothed),10))
	ax.minorticks_on()
	if i != 0:
		plt.setp( ax.get_xticklabels(), visible=False)
	if i == 0:
		ax.set_xlabel(r'z ($h^{-1}$ Mpc)')
	for item in (ax.get_yticklabels()):
		item.set_fontsize(6)
	ax.plot(galaxies_z[np.where((np.array(galaxies_x//(step/2)).astype('int') == i))],galaxies_y[np.where((np.array(galaxies_x//(step/2)).astype('int') == i))],'.',color='k')
	ax.set_xlim(np.array([0,340])+dist_to_center-170.)
	ax.set_ylim([0,24])
	for j in range(len(non_overlap_vcenter)):
		void = non_overlap_vcenter[j]
		radius = non_overlap_vradius[j]
		x = void[0]
		y = void[1]
		z = void[2]
		if int(x//(step/2)) == i:
			print i, x, y, z
			circle1 = plt.Circle((z,y), radius, color='k', fill=False)
			ax.add_artist(circle1)
fig.text(0.04, 0.5, r'y ($h^{-1}$ Mpc)', rotation='vertical')
plt.ion()
plt.show()
plt.savefig('clamato_06_16_voids_galaxies.pdf')

