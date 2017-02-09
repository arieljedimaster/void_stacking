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

from flipper import *

#importing functions from other files
from void_reduc import *
from void_filt import *
from void_stack import *

maps_pwd = '/Users/shared/ariel/Maps/Coadded/'

def ACT_upload():

    #uploading ACT cmb maps

	act_map_deep5 = liteMap.liteMapFromFits(maps_pwd+'deep5_S1_c7v5_night_nomoon_srcsub_srcpoint_noisecut_debuddyperdet_4way_set0123_2pass_wpoly_500_I.fits')
	act_map_deep6 = liteMap.liteMapFromFits(maps_pwd+'deep6_S1_c7v5_night_nomoon_srcsub_srcpoint_noisecut_debuddyperdet_4way_set0123_2pass_wpoly_500_I.fits')
	act_map_deep56_AR2 = liteMap.liteMapFromFits(maps_pwd+'deep56_S2_c7v5_AR2_night_nomoon_srcsub_srcpoint_noisecut_fixpickupcut_noturn_debuddyperdet_4way_set0123_2pass_wpoly_500_I.fits')
	act_map_deep56_AR1 = liteMap.liteMapFromFits(maps_pwd+'deep56_S2_c7v5_AR1_night_nomoon_srcsub_srcpoint_noisecut_fixpickupcut_noturn_debuddyperdet_4way_set0123_2pass_wpoly_500_I.fits')

	#old act map for testing purposes
	#act_map_tester = liteMap.liteMapFromFits(maps_pwd+'deep5_S1_c7v5_night_nomoon_srcsub_srcpoint_noisecut_debuddyperdet_4way_set0_2pass_wpoly_500_I.fits')
	#act_map_old = liteMap.liteMapFromFits('ACT_220_equ_season_3_1way_v3_summed.fits')

	return act_map_deep5, act_map_deep6, act_map_deep56_AR2, act_map_deep56_AR1

def ACT_split(act_map_deep5, act_map_deep6, act_map_deep56_AR2, act_map_deep56_AR1):

	#splitting up Deep5 and Deep56 map using submaps

	act_map_deep56_AR2_1 = act_map_deep56_AR2.selectSubMap(0.000005, act_map_deep56_AR2.x0 -0.000005,act_map_deep56_AR2.y0 + 0.000005,act_map_deep56_AR2.y1 -0.000005) #RA = 0 to x0
	act_map_deep56_AR2_2 = act_map_deep56_AR2.selectSubMap(act_map_deep56_AR2.x1+ 0.000005,360- 0.000005,act_map_deep56_AR2.y0 + 0.000005,act_map_deep56_AR2.y1 - 0.000005) #RA = x1 to 360

	act_map_deep56_AR1_1 = act_map_deep56_AR1.selectSubMap(0.000005, act_map_deep56_AR1.x0 -0.000005,act_map_deep56_AR1.y0 + 0.000005,act_map_deep56_AR1.y1 -0.000005) #RA = 0 to x0
	act_map_deep56_AR1_2 = act_map_deep56_AR1.selectSubMap(act_map_deep56_AR1.x1+ 0.000005,360- 0.000005,act_map_deep56_AR1.y0 + 0.000005,act_map_deep56_AR1.y1 - 0.000005) #RA = x1 to 360

	act_map_deep5_1 = act_map_deep5.selectSubMap(0.000005, act_map_deep5.x0 -0.000005,act_map_deep5.y0 + 0.000005,act_map_deep5.y1 -0.000005) #RA = 0 to x0
	act_map_deep5_2 = act_map_deep5.selectSubMap(act_map_deep5.x1+ 0.000005,360- 0.000005,act_map_deep5.y0 + 0.000005,act_map_deep5.y1 - 0.000005) #RA = x1 to 360

	return act_map_deep56_AR2_1, act_map_deep56_AR2_2, act_map_deep56_AR1_1, act_map_deep56_AR1_2, act_map_deep5_1, act_map_deep5_2

def ACT_crop(act_map_deep56_AR2_1, act_map_deep56_AR2_2, act_map_deep56_AR1_1, act_map_deep56_AR1_2, act_map_deep5_1, act_map_deep5_2, act_map_deep6):

    ############################################
    # Re cropping axes to get rid of the edges of the map
    # we think the edges might be contaminating the stacked map

    # we want to crop the maps to get rid of the regions around the edges
    # if I comment this section out I will still use the edges

    act_map_deep5_2 = act_map_deep5_2.selectSubMap(350.0, 360.0 -0.00005, -4, 4)

    act_map_deep5_1 = act_map_deep5_1.selectSubMap(0.0, 4.0, -4.0, 4.0)

    act_map_deep6 = act_map_deep6.selectSubMap(28.0, 41.0, -8.0, -1.0)

    act_map_deep56_AR1_1 = act_map_deep56_AR1_1.selectSubMap(0.0, 45.0, -8.0, 5.0)

    act_map_deep56_AR1_2 = act_map_deep56_AR1_2.selectSubMap(348.0, 360.0 -0.00005, -8.0, 5.0)

    act_map_deep56_AR2_1 = act_map_deep56_AR2_1.selectSubMap(0.0, 45.0, -8.0, 5.0)

    act_map_deep56_AR2_2 = act_map_deep56_AR2_2.selectSubMap(348.0, 360.0 -0.00005, -8.0, 5.0)

    return act_map_deep5_2, act_map_deep5_1, act_map_deep6, act_map_deep56_AR1_1, act_map_deep56_AR1_2, act_map_deep56_AR2_1, act_map_deep56_AR2_2, act_map_deep6

