#import pyfits
import numpy as np
import os
import matplotlib.pyplot as plt
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

pwd_void_cmass1 = '/Users/Ariel/Documents/2016 -Summer/Summer Research Project/Data/void_catalog_2015.03.31/sdss_dr10/redshift/sample_lss.dr10cmass1.dat'
pwd_void_cmass2 = '/Users/Ariel/Documents/2016 -Summer/Summer Research Project/Data/void_catalog_2015.03.31/sdss_dr10/redshift/sample_lss.dr10cmass2.dat'
pwd_void_cmass3 = '/Users/Ariel/Documents/2016 -Summer/Summer Research Project/Data/void_catalog_2015.03.31/sdss_dr10/redshift/sample_lss.dr10cmass3.dat'
pwd_void_lowz1 = '/Users/Ariel/Documents/2016 -Summer/Summer Research Project/Data/void_catalog_2015.03.31/sdss_dr10/redshift/sample_lss.dr10lowz1.dat'
pwd_void_lowz2 = '/Users/Ariel/Documents/2016 -Summer/Summer Research Project/Data/void_catalog_2015.03.31/sdss_dr10/redshift/sample_lss.dr10lowz2.dat'
pwd_void_lowz3 = '/Users/Ariel/Documents/2016 -Summer/Summer Research Project/Data/void_catalog_2015.03.31/sdss_dr10/redshift/sample_lss.dr10lowz3.dat'
pwd_void_lowz4 = '/Users/Ariel/Documents/2016 -Summer/Summer Research Project/Data/void_catalog_2015.03.31/sdss_dr10/redshift/sample_lss.dr10lowz4.dat'

#opening void positions

###cmass1

void_pos_cmass1 = np.loadtxt(pwd_void_cmass1+'/sky_positions_all_lss.dr10cmass1.dat.out', unpack=True)

void_RA_cmass1 = void_pos_cmass1[0]
void_dec_cmass1 = void_pos_cmass1[1]
void_r_eff_cmass1 = void_pos_cmass1[2]
void_ID_cmass1 = void_pos_cmass1[3]

void_centres_cmass1 = np.loadtxt(pwd_void_cmass1+'/centers_all_lss.dr10cmass1.dat.out', unpack=True)
void_cent_cmass1 = void_centres_cmass1[0]
void_vol_cmass1 = void_centres_cmass1[1]
void_r_eff_2_cmass1 = void_centres_cmass1[2]
void_ID_cmass1 = void_centres_cmass1[3]
void_dens_cmass1 = void_centres_cmass1[4] #void density contrast

###cmass2

void_pos_cmass2 = np.loadtxt(pwd_void_cmass2+'/sky_positions_all_lss.dr10cmass2.dat.out', unpack=True)

void_RA_cmass2 = void_pos_cmass2[0]
void_dec_cmass2 = void_pos_cmass2[1]
void_r_eff_cmass2 = void_pos_cmass2[2]

void_centres_cmass2 = np.loadtxt(pwd_void_cmass2+'/centers_all_lss.dr10cmass2.dat.out', unpack=True)
void_cent_cmass2 = void_centres_cmass2[0]
void_vol_cmass2 = void_centres_cmass2[1]
void_r_eff_2_cmass2 = void_centres_cmass2[2]
void_ID_cmass2 = void_centres_cmass2[3]
void_dens_cmass2 = void_centres_cmass2[4] #void density contrast

###cmass3

void_pos_cmass3 = np.loadtxt(pwd_void_cmass3+'/sky_positions_all_lss.dr10cmass3.dat.out', unpack=True)

void_RA_cmass3 = void_pos_cmass3[0]
void_dec_cmass3 = void_pos_cmass3[1]
void_r_eff_cmass3 = void_pos_cmass3[2]

void_centres_cmass3 = np.loadtxt(pwd_void_cmass3+'/centers_all_lss.dr10cmass3.dat.out', unpack=True)
void_cent_cmass3 = void_centres_cmass3[0]
void_vol_cmass3 = void_centres_cmass3[1]
void_r_eff_2_cmass3 = void_centres_cmass3[2]
void_ID_cmass3 = void_centres_cmass3[3]
void_dens_cmass3 = void_centres_cmass3[4] #void density contrast

###lowz1

void_pos_lowz1 = np.loadtxt(pwd_void_lowz1+'/sky_positions_all_lss.dr10lowz1.dat.out', unpack=True)

void_RA_lowz1 = void_pos_lowz1[0]
void_dec_lowz1 = void_pos_lowz1[1]
void_r_eff_lowz1 = void_pos_lowz1[2]

void_centres_lowz1 = np.loadtxt(pwd_void_lowz1+'/centers_all_lss.dr10lowz1.dat.out', unpack=True)
void_cent_lowz1 = void_centres_lowz1[0]
void_vol_lowz1 = void_centres_lowz1[1]
void_r_eff_2_lowz1 = void_centres_lowz1[2]
void_ID_lowz1 = void_centres_lowz1[3]
void_dens_lowz1 = void_centres_lowz1[4] #void density contrast

###lowz2

void_pos_lowz2 = np.loadtxt(pwd_void_lowz2+'/sky_positions_all_lss.dr10lowz2.dat.out', unpack=True)

void_RA_lowz2 = void_pos_lowz2[0]
void_dec_lowz2 = void_pos_lowz2[1]
void_r_eff_lowz2 = void_pos_lowz2[2]

void_centres_lowz2 = np.loadtxt(pwd_void_lowz2+'/centers_all_lss.dr10lowz2.dat.out', unpack=True)
void_cent_lowz2 = void_centres_lowz2[0]
void_vol_lowz2 = void_centres_lowz2[1]
void_r_eff_2_lowz2 = void_centres_lowz2[2]
void_ID_lowz2 = void_centres_lowz2[3]
void_dens_lowz2 = void_centres_lowz2[4] #void density contrast

###lowz3

void_pos_lowz3 = np.loadtxt(pwd_void_lowz3+'/sky_positions_all_lss.dr10lowz3.dat.out', unpack=True)

void_RA_lowz3 = void_pos_lowz3[0]
void_dec_lowz3 = void_pos_lowz3[1]
void_r_eff_lowz3 = void_pos_lowz3[2]

void_centres_lowz3 = np.loadtxt(pwd_void_lowz3+'/centers_all_lss.dr10lowz3.dat.out', unpack=True)
void_cent_lowz3 = void_centres_lowz3[0]
void_vol_lowz3 = void_centres_lowz3[1]
void_r_eff_2_lowz3 = void_centres_lowz3[2]
void_ID_lowz3 = void_centres_lowz3[3]
void_dens_lowz3 = void_centres_lowz3[4] #void density contrast

###lowz4

void_pos_lowz4 = np.loadtxt(pwd_void_lowz4+'/sky_positions_all_lss.dr10lowz4.dat.out', unpack=True)

void_RA_lowz4 = void_pos_lowz4[0]
void_dec_lowz4 = void_pos_lowz4[1]
void_r_eff_lowz4 = void_pos_lowz4[2]

void_centres_lowz4 = np.loadtxt(pwd_void_lowz4+'/centers_all_lss.dr10lowz4.dat.out', unpack=True)
void_cent_lowz4 = void_centres_lowz4[0]
void_vol_lowz4 = void_centres_lowz4[1]
void_r_eff_2_lowz4 = void_centres_lowz4[2]
void_ID_lowz4 = void_centres_lowz4[3]
void_dens_lowz4 = void_centres_lowz4[4] #void density contrast


#gather together all void RA/ DEC from all of the void files


void_RA = np.concatenate((void_RA_cmass1, void_RA_cmass2, void_RA_cmass3, void_RA_lowz1, void_RA_lowz2, void_RA_lowz3, void_RA_lowz4))
void_Dec = np.concatenate((void_dec_cmass1, void_dec_cmass2, void_dec_cmass3, void_dec_lowz1, void_dec_lowz2, void_dec_lowz3, void_dec_lowz4))
void_r_eff = np.concatenate((void_r_eff_cmass1, void_r_eff_cmass2, void_r_eff_cmass3, void_r_eff_lowz1, void_r_eff_lowz2, void_r_eff_lowz3, void_r_eff_lowz4))

void_cent = np.concatenate((void_cent_cmass1, void_cent_cmass2, void_cent_cmass3, void_cent_lowz1, void_cent_lowz2, void_cent_lowz3, void_cent_lowz4)) 
void_vol = np.concatenate((void_vol_cmass1, void_vol_cmass2, void_vol_cmass3, void_vol_lowz1, void_vol_lowz2, void_vol_lowz3, void_vol_lowz4))
void_r_eff_2 = np.concatenate((void_r_eff_2_cmass1, void_r_eff_2_cmass2, void_r_eff_2_cmass3, void_r_eff_2_lowz1, void_r_eff_2_lowz2, void_r_eff_2_lowz3, void_r_eff_2_lowz4))
void_ID = np.concatenate((void_ID_cmass1, void_ID_cmass2, void_ID_cmass3, void_ID_lowz1, void_ID_lowz2, void_ID_lowz3, void_ID_lowz4))
void_dens = np.concatenate((void_dens_cmass1, void_dens_cmass2, void_dens_cmass3, void_dens_lowz1, void_dens_lowz2, void_dens_lowz3, void_dens_lowz4))

#creating the matrix to save data in a txt file:

all_voids = [] #contains all info about all voids

for i in np.arange(0, len(void_RA)):
    all_voids.append([void_RA[i], void_Dec[i], void_r_eff[i], void_cent[i], void_vol[i], void_r_eff_2[i], void_ID[i], void_dens[i]])

np.savetxt('void_info.txt', all_voids)
