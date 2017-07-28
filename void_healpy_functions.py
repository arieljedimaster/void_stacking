# ALL FUNCTIONS FOR HEALPY STACKING

import pyfits
import numpy as np
import os

#import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import scipy
from scipy import stats
#from astLib import astWCS
import random as random
#from astLib import astCoords

# to time how long each step takes
import time


#importing healpy

import healpy as hp

#importing fits
from astropy.io import fits




def DeclRaToIndex(NSIDE,decl,RA):
    return hp.pixelfunc.ang2pix(NSIDE,np.radians(-decl+90.),np.radians(RA))

def DeclRaToThetaPhi(decl,RA):
    return np.radians(-decl+90.),np.radians(RA)

def IndexToDeclRa(index, NSIDE):
    theta,phi=hp.pixelfunc.pix2ang(NSIDE,index)
    return np.degrees(pi/2.-theta),np.degrees(phi)

def CartviewRA(RA):
    RA[RA>=180.] =  RA[RA>=180.] - 360.
    return RA


def get_annulus(RA,Dec,R,hp_map,nside):
    theta,phi = DeclRaToThetaPhi(Dec,RA)
    vec = hp.ang2vec(theta,phi)
    R_rad_inner = (np.pi/180.)*R
    R_rad_outer = R_rad_inner*np.sqrt(2)
    pix_disc_out = hp.query_disc(nside,vec,R_rad_outer,inclusive = True)
    pix_disc_inner = hp.query_disc(nside,vec,R_rad_inner, inclusive = True)
    pix_annulus = np.setdiff1d(pix_disc_out,pix_disc_inner)
    vals_annulus = hp_map[pix_annulus]
    return pix_annulus, vals_annulus

def get_annulus_galactic(phi,theta,R,hp_map,nside):
    vec = hp.ang2vec(theta,phi)
    R_rad_inner = (np.pi/180.)*R
    R_rad_outer = R_rad_inner*np.sqrt(2)
    pix_disc_out = hp.query_disc(nside,vec,R_rad_outer,inclusive = True)
    pix_disc_inner = hp.query_disc(nside,vec,R_rad_inner, inclusive = True)
    pix_annulus = np.setdiff1d(pix_disc_out,pix_disc_inner)
    vals_annulus = hp_map[pix_annulus]
    return pix_annulus, vals_annulus

def get_void_region(RA,Dec,R,hp_map,nside):
    theta,phi = DeclRaToThetaPhi(Dec,RA)
    vec = hp.ang2vec(theta,phi)
    R_rad = (180./np.pi)*R
    pix_disc = hp.query_disc(nside,vec,R_rad, inclusive = True)
    vals_disc = hp_map[pix_disc]
    return pix_disc, vals_disc
def get_void_region_galactic(phi,theta,R,hp_map,nside):
    vec = hp.ang2vec(theta,phi)
    R_rad = (180./np.pi)*R
    pix_disc = hp.query_disc(nside,vec,R_rad, inclusive = True)
    vals_disc = hp_map[pix_disc]
    return pix_disc, vals_disc

def chunks(seq, num):
  avg = len(seq) / float(num)
  out = []
  last = 0.0

  while last < len(seq):
    out.append(seq[int(last):int(last + avg)])
    last += avg

  return out

def data_radial_binning(condition, void_r, void_RA, void_Dec, radial_bins):
    '''note: the input void_RA, void_Dec, and void_r are LISTS'''
    #if the first argument in the function radialbin() is 'effective radius', it means we'll bin by effective radius
    #if the argument is 'angular radius', we'll bin by angular radius
    
    if condition == 'angular radius':
        #meaning we're binning by the angular radius

        #set up a dictionary which related void_index to other void properties
        #dictionary of void_RA, void_Dec, void_r
        void_dict = dict((z[0],list(z[1:])) for z in zip(void_r,void_RA,void_Dec))
        void_r_sort = sorted(void_r, key=float)
        void_r_split = chunks(void_r_sort, radial_bins)
        #now len(void_r_split) should be n
        #empty pre beinned void_RA and Dec lists, based off of the input bin number
        void_RA_split = [[] for _ in range(radial_bins)]
        void_Dec_split = [[] for _ in range(radial_bins)]
        void_r_avg = []

    
        for i in np.arange(0,len(void_r_split)):
            #compute average void_r for each bin
            #this will be used for the R_filt = 0.7*void_r_avg
            void_r_avg.append(np.mean(void_r_split[i]))
            for k in np.arange(0, len(void_r_split[i])):
                #splitting up void_RA based on bin
                void_RA_split[i].append(void_dict[void_r_split[i][k]][0])
                void_Dec_split[i].append(void_dict[void_r_split[i][k]][1])

        return void_r_split, void_RA_split, void_Dec_split, void_r_avg





