import numpy as np
from astropy.cosmology import FlatLambdaCDM
from astropy.coordinates import SkyCoord
import astropy.units as u

from astropy.io import fits
from astropy.io import ascii
from astropy.table import Table

from fnmatch import fnmatch

cosmos_specz_fil = '/Users/ALEX/Berkeley/COSMOS_Spec/'+\
'OBSERVED_TARGETS_15April2015_withHeader.dat'

cosmos = ascii.read(cosmos_specz_fil,format='fast_commented_header',header_start=-1)

print(cosmos.dtype)
#(set(cosmos['Instr']))

zmin = 2.0
zmax = 2.6
ra0 = 149.975
dec0 = 2.15
ra1 = 150.22856
dec1 = 2.4880706 

cosmo = FlatLambdaCDM(H0=70, Om0=0.31)

deg_per_hMpc = 1./ cosmo.h / cosmo.comoving_distance(2.35).value * 180./np.pi

# Look for galaxies up to 15Mpc/h transverse from the central map region
d_ra  = 15. * deg_per_hMpc
d_dec = 15. * deg_per_hMpc

# make redshift and coordinate cuts
zcut = np.all(np.column_stack([(cosmos['z_spec'] >= zmin) ,
                            (cosmos['z_spec'] < zmax),
                            (cosmos['ra'] >= ra0-d_ra),
                            (cosmos['ra'] <= ra1+d_ra),
                              (cosmos['dec'] >= dec0-d_dec),
                              (cosmos['dec'] <= dec0+d_dec)]),axis=1)              
              
cosmos = cosmos[zcut]
print('%i galaxies after redshift and positional cuts' % len(cosmos))
#print(np.shape(cosmos))
# print(cosmos.meta) # These are the comments

#fig, ax = plt.subplots()
zspec = cosmos['z_spec']
binwidth=0.025
#plt.hist(zspec,bins=np.arange(min(zspec), max(zspec) + binwidth, binwidth))
#plt.show()

#[fnmatch(instr,'MOS*') for instr in cosmos['Instr']]

# get high-confidence (>3) redshifts from zCOSMOSDeep for galaxies only
confcut = np.all(np.column_stack(( (cosmos['Instr'] == 'zDEEP'),
									#(cosmos['z_spec'] > 2.15),
									#(cosmos['z_spec'] <= 2.54),
                                   (cosmos['Q_f'] >= 3.),
                                   (cosmos['Q_f'] < 10.) )), axis=1)

cosmos = cosmos[confcut]
print('%i galaxies within zCOSMOS-Deep after quality cuts' % len(cosmos))

# Match coordinates with MOSDEF sample, and throw away any that match
zDeepCoord = SkyCoord(ra=cosmos['ra']*u.degree,dec=cosmos['dec']*u.degree)

'''idzD, d2d , d3d = zDeepCoord.match_to_catalog_sky(MosdefCoord)

SepUnit= u.arcmin
#plt.hist(d2d.to(SepUnit))
#plt.title('Separation from MOSDEF Objects')
#plt.xlabel("Separation "+"({0.unit:s})".format(1.*SepUnit))
#plt.show()

closematch_zD = d2d < (1*u.arcsec)
notmatch_zD = np.invert(closematch_zD)

print('%i zDeep galaxies are also measured by MOSDEF' % sum(1 for i in closematch_zD if i))

cosmos = cosmos[notmatch_zD]
zDeepCoord = zDeepCoord[notmatch_zD]
print('This leaves %i zDeep galaxies in our sample' % len(cosmos))

fig, ax = plt.subplots()
binwidth=0.025
ax.hist(cosmos['z_spec'],bins=np.arange(min(cosmos['z_spec']),
                                               max(cosmos['z_spec']) + binwidth, binwidth))
ax.set_xlabel(r'$z_{spec}$')
ax.set_ylabel('N')
plt.show()'''

comdist0 = cosmo.comoving_distance(2.15)
binfac = 2.
dcomdist_dz = 2997./cosmo.H(2.35)
xpos = comdist0*(cosmos['ra']-ra0)/np.cos(0.5*np.abs(cosmos['dec']+dec0)*np.pi/180.)*binfac*np.pi/180.
ypos = comdist0*(cosmos['dec']-dec0)*binfac*np.pi/180.
zpos = (cosmos['z_spec']-zmin)*dcomdist_dz*binfac

in_vol = np.where((xpos.value > 1) & (xpos.value < 35) & (ypos.value > 1) & (ypos.value < 47) & (zpos.value > 1) & (zpos.value < 679))
print np.shape(in_vol)