import matplotlib.gridspec as gridspec
#gs = gridspec.GridSpec(13,12)
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
dec1 = 2.4880706 
deg_per_hMpc = 1./ cosmo.h / cosmo.comoving_distance(zmid).value * 180./np.pi
dist_to_center = cosmo.comoving_distance(zmid).value*cosmo.h

# Load the map
mapfile = 'map_2017tmp_sm2.0.bin'
map_smoothed, allx, ally, allz = load_and_smooth_map(mapfile,(48,48,680))

# Load the galaxies and convert coordinates to h^{-1} Mpc
galaxies = ascii.read('galxcorr_cl2016_v0_with_vuds.dat')
galaxies_x = np.cos(0.5*(dec0+dec1)*np.pi/180.)*(galaxies['ra'] - ra0)/deg_per_hMpc
galaxies_y = (galaxies['dec'] - dec0)/deg_per_hMpc
galaxies_z = cosmo.comoving_distance(galaxies['zspec']).value*cosmo.h

# Select certain galaxies
#galaxies_x = galaxies_x[(galaxies['source'] == 'MOSDEF')]
#galaxies_y = galaxies_y[(galaxies['source'] == 'MOSDEF')]
#galaxies_z = galaxies_z[(galaxies['source'] == 'MOSDEF')]

galaxies_x_mosdef = galaxies_x[(galaxies['source'] == 'MOSDEF')]
galaxies_y_mosdef = galaxies_y[(galaxies['source'] == 'MOSDEF')]
galaxies_z_mosdef = galaxies_z[(galaxies['source'] == 'MOSDEF')]

galaxies_x_vuds = galaxies_x[(galaxies['source'] == 'VUDS')]
galaxies_y_vuds = galaxies_y[(galaxies['source'] == 'VUDS')]
galaxies_z_vuds = galaxies_z[(galaxies['source'] == 'VUDS')]


#galaxies_x = galaxies_x[(galaxies['source'] == 'zDEEP')]
#galaxies_y = galaxies_y[(galaxies['source'] == 'zDEEP')]
#galaxies_z = galaxies_z[(galaxies['source'] == 'zDEEP')]

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

#fig = plt.figure(figsize=(6.09,8))
for j in range(12):
	step = 4
	for i in range(12):
		if i == j:
			plt.clf()
			fig = plt.figure(figsize=(14.17,2.29))
			#ax = plt.subplots(111)
			ax = plt.gca()
			ax.contourf(allz[step*i,:,:]+dist_to_center-170.,ally[step*i,:,:],map_smoothed[step*i,:,:],cmap=plt.cm.jet_r,levels=np.linspace(np.min(map_smoothed),np.max(map_smoothed),10),vmin=-0.25,vmax=0.25)
			ax.minorticks_on()
			ax2 = ax.twiny()
			ax2.minorticks_on()
			ax2.set_xlim([2.15875,2.55805])
			ax2.set_xticks([2.2,2.25, 2.3,2.35,2.4,2.45,2.5])
			ax.set_xlabel(r'Comoving distance ($h^{-1}$ Mpc)',size=18)
			ax.set_ylabel(r'y$_{\textrm{perp}}$ ($h^{-1}$ Mpc)',size=18)
			ax2.set_xlabel(r'redshift',size=18)
			for item in (ax.get_yticklabels()):
				item.set_fontsize(12)
			for item in (ax2.get_xticklabels()):
				item.set_fontsize(12)
			for item in (ax.get_xticklabels()):
				item.set_fontsize(12)
			ax.plot(galaxies_z_mosdef[np.where((np.array(galaxies_x_mosdef//(step/2)).astype('int') == i))],galaxies_y_mosdef[np.where((np.array(galaxies_x_mosdef//(step/2)).astype('int') == i))],'.',color='k',markersize=12)
			ax.plot(galaxies_z_vuds[np.where((np.array(galaxies_x_vuds//(step/2)).astype('int') == i))],galaxies_y_vuds[np.where((np.array(galaxies_x_vuds//(step/2)).astype('int') == i))],'*',color='g',markeredgewidth=0.0,markersize=12)
			ax.set_xlim(np.array([0,340])+dist_to_center-170.)
			ax.set_ylim([0,24])
			ax.set_yticks([0,8,16,24])
			#plt.figaspect(1/14.17)
			'''for j in range(len(non_overlap_vcenter)):
				void = non_overlap_vcenter[j]
				radius = non_overlap_vradius[j]
				x = void[0]
				y = void[1]
				z = void[2]
				#if ((x+2*radius >= 2*i) and (x+2*radius <= 2*(i+1))
				#	or (x-2*radius >= 2*i) and (x-2*radius <= 2*(i+1))):
				if (int(x//(step/2)) == i or int(x//(step/2)) == i-1 or int(x//(step/2)) == i+1
						or int(x//(step/2)) == i+2 or int(x//(step/2)) == i-2
						or int(x//(step/2)) == i+3 or int(x//(step/2)) == i-3):
				#if (2*i > x - radius) and (2*i < x + radius):
					print 2*i, x, y, z, radius
					circle1 = plt.Circle((z,y), radius, color='k', fill=False)
					ax.add_artist(circle1)'''
	#fig.text(0.04, 0.5, r'y ($h^{-1}$ Mpc)', rotation='vertical')

	#Create custom artists
	simArtist = plt.Line2D((0,1),(0,0), color='k', marker='.', linestyle='',markersize=12)
	anyArtist = plt.Line2D((0,1),(0,0), color='g', marker='*', linestyle='',markersize=12)

	#Create legend from custom artist/label lists
	ax.legend([simArtist,anyArtist],
			  ['MOSDEF', 'VUDS'],numpoints=1,bbox_to_anchor=(0.98, 1.02),bbox_transform=plt.gcf().transFigure,frameon=False)

	#fig.subplots_adjust(hspace=0.5)
	plt.tight_layout()

	plt.ion()
	plt.show()
	plt.savefig('clamato_03_17_galaxies_mosdef_vuds_slice%i.pdf' % j)

