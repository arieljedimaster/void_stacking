import pyfits
import numpy as np
import os

#import matplotlib.pyplot as plt
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
import scipy
from scipy import stats
#from astLib import astWCS
import random as random
#from astLib import astCoords

from void_cat_load import *


from flipper import *

'''
map_coadd_folder ='/Users/shared/ariel/Maps/Coadded/'
map_folder = '/Users/shared/ariel/Maps/'
'''

#map_folder = '/Users/arielamaral/Code/void_code/extracted_BossN/'

map_folder = '/Users/shared/ariel/Maps/Coadded/'

#deep56_S2_c7v5_AR1_night_nomoon_srcsub_srcpoint_noisecut_fixpickupcut_noturn_debuddyperdet_4way_set0123_2pass_noise.fits
#deep56_S2_c7v5_AR1_night_nomoon_srcsub_srcpoint_noisecut_fixpickupcut_noturn_debuddyperdet_4way_set0123_2pass_weights_I.fits
#deep56_S2_c7v5_AR1_night_nomoon_srcsub_srcpoint_noisecut_fixpickupcut_noturn_debuddyperdet_4way_set0123_2pass_wpoly_500_I.fits
#deep56_S2_c7v5_AR2_night_nomoon_srcsub_srcpoint_noisecut_fixpickupcut_noturn_debuddyperdet_4way_set0123_2pass_noise.fits
#deep56_S2_c7v5_AR2_night_nomoon_srcsub_srcpoint_noisecut_fixpickupcut_noturn_debuddyperdet_4way_set0123_2pass_weights_I.fits
#deep56_S2_c7v5_AR2_night_nomoon_srcsub_srcpoint_noisecut_fixpickupcut_noturn_debuddyperdet_4way_set0123_2pass_wpoly_500_I.fits
#deep5_S1_c7v5_night_nomoon_srcsub_srcpoint_noisecut_debuddyperdet_4way_set0123_2pass_noise.fits
#deep5_S1_c7v5_night_nomoon_srcsub_srcpoint_noisecut_debuddyperdet_4way_set0123_2pass_weights_I.fits
#deep5_S1_c7v5_night_nomoon_srcsub_srcpoint_noisecut_debuddyperdet_4way_set0123_2pass_wpoly_500_I.fits
#deep6_S1_c7v5_night_nomoon_srcsub_srcpoint_noisecut_debuddyperdet_4way_set0123_2pass_noise.fits
#deep6_S1_c7v5_night_nomoon_srcsub_srcpoint_noisecut_debuddyperdet_4way_set0123_2pass_weights_I.fits
#deep6_S1_c7v5_night_nomoon_srcsub_srcpoint_noisecut_debuddyperdet_4way_set0123_2pass_wpoly_500_I.fits


#map = 'extracted_boss_tot_ar2_night_tot_sky_div.fits'

output_folder = '/Users/arielamaral/Code/void_code/'

map ='deep6_S1_c7v5_night_nomoon_srcsub_srcpoint_noisecut_debuddyperdet_4way_set0123_2pass_wpoly_500_I.fits'

deep6 = liteMap.liteMapFromFits(map_folder+map)

deep6 = deep6.selectSubMap(28.0, 41.0, -8.0, -1.0)

#void_RA, void_Dec, void_z, void_rlos, void_r_eff, void_r_ang, void_los_size, void_dens, void_f_vol = get_void_info_from_fits()

void_RA, void_Dec, void_r_ang, void_ID, void_e_f = np.loadtxt('void_info_2.txt', unpack=True)



plt.figure()
plt.title('Deep 6 Cropped with voids')
plt.plot(void_RA, void_Dec, 'wo')
deep6.plot(valueRange=[-500,500])
plt.savefig(output_folder+'cropped_deep6_voids.png')



'''

map = 'deep6_S1_c7v5_night_nomoon_srcsub_srcpoint_noisecut_debuddyperdet_4way_set0123_2pass_wpoly_500_I.fits'



flipperfied_map = liteMap.liteMapFromFits(map_coadd_folder+map)

flipperfied_map = flipperfied_map.selectSubMap(28.0, 41.0, -8.0, -1.0)



el = numpy.arange(10000)
Fel = (1. - numpy.exp(-el**2/(2*10.**2)))*numpy.exp(-el**2/(2*300.**2))
filteredMap = flipperfied_map.filterFromList([el,Fel])

'''
'''
#upload new maps
#start with just uploading one of the maps

act_map0010_raw = liteMap.liteMapFromFits(map_folder+'extracted_boss_tot_ar2_night_tot_sky_map0010.fits')

#inverse variance weight map
act_div = liteMap.liteMapFromFits(map_folder+'extracted_boss_tot_ar2_night_tot_sky_div.fits')

#inverse variance weight


act_map0010 = act_map0010_raw.copy()

act_map0010.data = act_map0010_raw.data#*act_div.data
'''




#flipperfied_map0 = liteMap.liteMapFromFits(map_folder+map0)

#diffmap = flipperfied_map.copy()
#print 'co add map'
#print flipperfied_map.info()
#diffmap.data[:] = flipperfied_map.data[:] - flipperfied_map0.data[:]
#print 'original map'
#print flipperfied_map0.info()

#print 'diffmap'
#print diffmap.info()

'''

output_folder = '/Users/arielamaral/Code/void_code/'

plt.figure()
plt.title('Cropped map0010')
act_map0010.plot(valueRange = [-500,500])
plt.savefig(output_folder+'map0010_cropped.png')
'''
'''

plt.figure()
plt.title('Deep 6 Map')
flipperfied_map.plot(valueRange=[-500,500])
plt.savefig(output_folder+'unfiltered_deep6.png')
'''
'''

void_RA, void_Dec, void_z, void_rlos, void_r_eff, void_r_ang, void_los_size, void_dens, void_f_vol = get_void_info_from_fits()

#void_RA, void_Dec, void_r_ang, void_ID, void_e_f = np.loadtxt('void_info_2.txt', unpack=True)

plt.figure()
#plt.title('Locations of Voids on the Sky')
plt.plot(void_RA, void_Dec, 'm*')
plt.xlabel('Right Ascension')
plt.ylabel('Declination')
plt.savefig(output_folder+'void_locations_on_sky_new.png')
'''




