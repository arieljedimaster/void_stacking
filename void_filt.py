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

'''
All functions which relate to filtering and reducing voids ready for stacking analysis
uses flipper to do things.
'''

def mask(ra, dec, m, r):
    '''
    masks out data within a disk, this does the same as the
    flipper mask function
    '''
    x,y = m.skyToPix(ra,dec)
    dat = m.data
    array_size = dat.shape
    mask = np.zeros(array_size)
    num_pix = 0 #number of non-zero pixels
    for i in np.arange(0,len(mask[0])):
        for j in np.arange(0, len(mask)):
            dist = np.sqrt((i-x)**2 + (j-y)**2)
            if dist < r:
                mask[j][i] = 1
                num_pix +=1
    return mask*dat, num_pix

def generate_smap(void_RA, void_Dec, bound_x, bound_y, m):
    '''
    Returns a list of submaps stamps around void locations on
    the desired CMB map
    '''
    #smap = m.selectSubMap(void_RA - bound_x, void_RA + bound_x, void_Dec - bound_y, void_Dec + bound_y, safe = False)
    smap = m.selectSubMap(void_RA - bound_x, void_RA + bound_x, void_Dec - bound_y, void_Dec + bound_y)
    return smap

def generate_filtsmap(smap):
    '''
    Returns a list of filtered submap stamps around the void locations
    from smap_list()
    '''
    el = numpy.arange(10000)
    Fel = (1. - numpy.exp(-el**2/(2*10.**2)))*numpy.exp(-el**2/(2*300.**2))
    #Fel = 1/(1 + np.exp(-(el-1000))) # I just changed this, previously this was 1000
    filteredMap = smap.filterFromList([el,Fel])
    return filteredMap

def delta_T_list(void_RA, void_Dec, filt_smaps, r):
    '''
    r in pixels
    '''
    delta_T = []
    for i in np.arange(0, len(void_RA)):
        maskk1, num_pix1 = mask(void_RA[i], void_Dec[i], filt_smaps[i], r)
        avg1 = np.sum(maskk1)/num_pix1
        maskk2, num_pix2 = mask(void_RA[i], void_Dec[i], filt_smaps[i], np.sqrt(2)*r)
        annulus = maskk2 - maskk1
        annulus_pix = num_pix2 - num_pix1
        avg_annulus = np.sum(annulus)/annulus_pix
        delta_T.append(avg1 - avg_annulus)
    return delta_T

def delt_T_filt(void_RA, void_Dec, r, m):
    '''
    returns a list of a list of delta_T filters meant to be applied to
    mode filtered maps
    '''
    #change r to r_pix
    #input r in degrees
    r_arcmin = 60*r #degrees to arcmins
    pix_scale = (m.pixScaleX)*((180/np.pi) *(60))
    r_pix = r_arcmin / pix_scale 
    mapout = m.mask(void_RA,void_Dec,r_pix,mask_lo=15, mask_hi=25)
    mapin = m.mask(void_RA,void_Dec,np.sqrt(2)*r_pix,mask_lo=15, mask_hi=25)
    outer_circle = m.copy()
    outer_circle.data = m.data - mapin.data
    inner_circle = m.copy()
    inner_circle.data = m.data - mapout.data
    diffmap = inner_circle.copy()
    diffmap.data[:]= outer_circle.data[:]-inner_circle.data[:]
    return diffmap, inner_circle

def num_pix_within_annulus():
    '''
    returns the number of pixels within the masked region
    '''
    return

def num_pix_within_circle():
    '''
    returns the number of pixels within the masked region
    '''
    return

def mod_filt(filtmap, delTmap):
    '''
    takes in the mode filtered maps, and the delta T filter and subtracts them
    to get a modified T and filtered submap ready for stacking
    '''
    modmap = filtmap.copy()
    #fix the annulus mean bit!
    num_pix_index = np.where(delTmap.data != 0.)
    num_pix = len(num_pix_index[0])
    annulus_mean = np.sum(delTmap.data)/num_pix
    modmap.data = filtmap.data - annulus_mean
    return modmap


def del_T_measurement(inner_circle_list, diff_list):
    '''
    Similar to delta_T_list, but I'm using the mask function in flipper
    instead of the one I made above!
    '''
    del_T_list = []
    for i in np.arange(0, len(diff_list)):
        inner_circ = inner_circle_list[i]
        #num_pix_circ_i = np.where(inner_circ.data != 0.) #indeces where the data is non zero
        #num_pix_circ = len(num_pix_circ_i[0])
        num_pix_circ = len(np.transpose(np.nonzero(inner_circ.data)))
        inner_circ_mean = np.sum(inner_circ.data)/num_pix_circ
        annulus = diff_list[i]
        #num_pix_annulus_i = np.where(annulus.data != 0.) #indeces where this occurs
        #num_pix_annulus = len(num_pix_annulus_i[0])
        num_pix_annulus = len(np.transpose(np.nonzero(annulus.data)))
        annulus_mean = np.sum(annulus.data)/num_pix_annulus
        del_T = inner_circ_mean - annulus_mean
        del_T_list.append(del_T)
    return del_T_list
    

