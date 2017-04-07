import numpy as np
import scipy.ndimage
import matplotlib.gridspec as gridspec
import matplotlib
from matplotlib import font_manager
import matplotlib.pyplot as plt
from astropy.cosmology import FlatLambdaCDM
from load_and_smooth_map import *

cosmo = FlatLambdaCDM(H0=70, Om0=0.31)

zlya = 2.35
dist_to_center = cosmo.comoving_distance(zlya).value*cosmo.h

mapfile = '/Users/ALEX/Berkeley/IGM_Nyx_CLAMATO/CLAMATO_06_16/map.bin'
map_smoothed, allx, ally, allz = load_and_smooth_map(mapfile,(36,48,680))

'''Plot the voids on top of slices of the map.  Each slice is separated by 2 h^{-1} Mpc'''

voids = np.loadtxt('voids_2mpc_smoothing_constant_boundary.txt')
non_overlap_vradius = voids[:,0]
non_overlap_vcenter = voids[:,1:]
	

gs = gridspec.GridSpec(10,9)


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
	for j in range(len(non_overlap_vcenter)):
		void = non_overlap_vcenter[j]
		radius = non_overlap_vradius[j]
		x = void[0]
		y = void[1]
		z = void[2]+dist_to_center-170.
		if int(x//(step/2)) == i:
			print i, x, y, z
			circle1 = plt.Circle((z,y), radius, color='k', fill=False)
			ax.add_artist(circle1)
			
fig.text(0.04, 0.5, r'y ($h^{-1}$ Mpc)', rotation='vertical')
plt.savefig('clamato_06_16_with_voids_constant_boundary.pdf')