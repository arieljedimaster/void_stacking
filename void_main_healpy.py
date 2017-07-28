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

data_pwd = '/Users/arielamaral/Documents/void_stacking_healpy/Data/'
code_pwd = '/Users/arielamaral/Documents/void_stacking_healpy/Code/'
plots_pwd = '/Users/arielamaral/Documents/void_stacking_healpy/Plots/'


#upload Planck Maps 

planck_nilc = fits.open(data_pwd+'COM_CompMap_YSZ_R2.00.fits/nilc_ymaps.fits')

planck_nilc_header = planck_nilc[1].header

planck_nilc_data = planck_nilc[1].data

planck_nilc_full = planck_nilc_data['FULL    ']

planck_nilc_first = planck_nilc_data['FIRST   ']

planck_nilc_last = planck_nilc_data['LAST    ']

planck_nilc_mask = planck_nilc_data['MASK    ']


plt.figure()
hp.mollview(planck_nilc_full, min = -0.000001, max = +0.000001)
plt.title('Planck NILC "FULL" map')
plt.savefig(plots_pwd + 'planck_nilc_full_plot.png')
#plt.show()

plt.figure()
hp.mollview(planck_nilc_first, min = -0.000001, max = +0.000001)
plt.title('Planck NILC "FIRST" map')
plt.savefig(plots_pwd + 'planck_nilc_first_plot.png')
#plt.show()

plt.figure()
hp.mollview(planck_nilc_last, min = -0.000001, max = +0.000001)
plt.title('Planck NILC "LAST" map')
plt.savefig(plots_pwd + 'planck_nilc_last_plot.png')
#plt.show()

plt.figure()
hp.mollview(planck_nilc_mask)
plt.title('Planck NILC "MASK" map')
plt.savefig(plots_pwd + 'planck_nilc_mask_plot.png')
#plt.show()


#mask the full healpy map using mask

planck_nilc_full_masked = hp.pixelfunc.ma(planck_nilc_full)
planck_nilc_full_masked.mask = 1 - planck_nilc_mask

plt.figure()
hp.mollview(planck_nilc_full_masked, min = np.ma.min(planck_nilc_full_masked), max = np.ma.max(planck_nilc_full_masked))
plt.title('Planck NILC masked "FULL" map')
plt.savefig(plots_pwd + 'planck_nilc_full_masked_plot.png')
#plt.show()


#healpix info

nside= 2048 # found from planck_nilc[1].header
npix = 12*nside**2
res = int(np.log(nside)/np.log(2))

#pixel indeces for planck maps
planck_nilc_pix = np.arange(0, npix+1)

map_healpy_vals = planck_nilc_full_masked #just to make it easier, this is the map we're using for now


#load in void info
#Nadathur (2016) BOSS catalogue from: http://www.icg.port.ac.uk/stable/nadathur/voids/

void_RA, void_Dec, void_r_ang, void_ID, void_e_f = np.loadtxt(data_pwd+'void_info_2.txt', unpack=True)

#loading additional void info (from that same catalogue - can't remember why I put it into two different arrays)

void_ID, void_centre_x, void_centre_y, void_centre_z, void_min_dens, void_WtdAvgDens, void_r_eff, void_dens025, void_edge_flag, void_dens_ratio = np.loadtxt(data_pwd+'void_info_BOSS.txt', unpack = True)

'''
#lets try first with only like 100 sources
void_RA = void_RA[:500]
void_Dec = void_Dec[:500]
void_r_ang = void_r_ang[:500]
'''




print "Total length of the catalogue: ", len(void_RA)


#want to change void_RA and void_Dec into healpy pixel locations
#what pixel each void falls in    
void_pix_index = DeclRaToIndex(nside,void_Dec,void_RA)


num_in_stack = 0 #starting the count

#empty lists of smaps
smaps = [] #[[]for _ in range(void_RA)]

len_smaps = np.zeros(len(void_RA))

stacked_smap = 'hello'

total_voids_num = 0

#plotting fake data and fake map on top of eachother:

plt.figure()
hp.cartview(map_healpy_vals, min = -0.000001, max = +0.000001)
plt.plot(CartviewRA(void_RA),void_Dec, 'g.', label = 'Voids')
plt.title('Planck NILC masked full array with BOSS void locations')
plt.legend()
plt.savefig(plots_pwd+'planck_nilc_full_masked_plot_voids.png')






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

def get_void_region(RA,Dec,R,hp_map,nside):
    theta,phi = DeclRaToThetaPhi(Dec,RA)
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


radial_bins = 10 #??? can change this if I want

dT_radial = np.zeros(radial_bins)
dT_radial_err = np.zeros(radial_bins)

#sorting voids by radial bins

void_r_ang_split, void_RA_split, void_Dec_split, void_r_avg_split= data_radial_binning('angular radius', void_r_ang, void_RA, void_Dec, radial_bins)


#looping through each individual void which falls in the map boundary

for j in np.arange(0, radial_bins):

    print "On void radial bin # " + str(j+1)
    void_RA_bin = void_RA_split[j] #voids in that radial bin
    void_Dec_bin = void_Dec_split[j]
    void_r_ang_bin = void_r_ang_split[j]

    #setting size of smap and also size of annulus:

    R = 0.7*void_r_avg_split[j] #0.7 is from nadathur et al I think
    R_smap = void_r_avg_split[j]*3.0

    dT_array = np.zeros(len(void_RA_bin))

    #deleting all voids which fall void_RA + R_smap or void_Dec + R_smap of the boundary
    del_indeces_RA = np.argwhere((void_RA_bin >= (360. - R_smap)) ^ (void_RA_bin <= (0. + R_smap)))[:,0]
    del_indeces_Dec = np.argwhere(((void_Dec_bin + R_smap) >= 90.) ^ ((void_Dec_bin - R_smap) <= -90.))[:,0]
    
    #concatenate these arrays
    del_indeces = np.concatenate((del_indeces_RA,del_indeces_Dec))
    
    void_RA_bin = np.delete(void_RA_bin, del_indeces)
    void_Dec_bin = np.delete(void_Dec_bin, del_indeces)
    void_r_ang_bin = np.delete(void_r_ang_bin, del_indeces)



    for k in np.arange(0,len(void_RA_bin)):

        if (k % 100 == 0) & (k!=0):
            print "Stacking Void # " + str(k) + ", we are " + str(float(k)/float(len(void_RA_bin))*100.) + " percent done for this bin."

        RA_min = void_RA_bin[k] - R_smap
        RA_max = void_RA_bin[k] + R_smap
        Dec_min = void_Dec_bin[k] - R_smap
        Dec_max = void_Dec_bin[k] + R_smap

        void_Dec_rad, void_RA_rad = DeclRaToThetaPhi(void_Dec_bin[k],void_RA_bin[k])

        Dec_min_rad, RA_min_rad = DeclRaToThetaPhi(Dec_min,RA_min)
        Dec_max_rad, RA_max_rad = DeclRaToThetaPhi(Dec_max,RA_max)

        #get vertice vectors for the postage stamp from the RA/Dec

        vect1 = hp.ang2vec(Dec_min_rad,RA_min_rad)
        vect2 = hp.ang2vec(Dec_min_rad,RA_max_rad)
        vect3 = hp.ang2vec(Dec_max_rad,RA_max_rad)
        vect4 = hp.ang2vec(Dec_max_rad,RA_min_rad)

        #put into format for other function

        vertices = np.array([vect1,vect2,vect3,vect4])

        #get polygon of pixels arround void pix

        pix_matrix = hp.query_polygon(nside, vertices, inclusive=True)

        vals_matrix = map_healpy_vals[pix_matrix]

        #get rid of any voids that fall within that masked region
        #a quick way to fix it for now.
        if np.any(np.isnan(vals_matrix)) == False: #there are no nan values

            #not going the append this for now - will get error

            #smaps.append(vals_matrix)

            annulus_pix, annulus_vals = get_annulus(void_RA_bin[k],void_Dec_bin[k],R,map_healpy_vals,nside)


            void_region_pix, void_region_vals = get_void_region(void_RA_bin[k],void_Dec_bin[k],R,map_healpy_vals,nside)


            dT_reduced = np.average(void_region_vals) - np.average(annulus_vals)

            dT_array[k] = dT_reduced

            num_in_stack +=1
            total_voids_num += 1
        else:
            continue

        #len_smaps[k] = len(vals_matrix)
    print " "
    print "The average deltaT for this radial bin is: "+ str(np.average(dT_array))+ " +/- " + str(np.std(dT_array))
    print " "

    dT_radial[j] = np.average(dT_array)
    dT_radial_err[j] = np.std(dT_array)

#print len_smaps

#print "Maximum # pixels in stamp: ", np.amax(len_smaps)
#print "Minimum # pixels in stamp: ", np.amin(len_smaps)



plt.figure()
plt.title(r"Void $\Delta$T profile for " + str(radial_bins) + " radial bins")
plt.plot(void_r_avg_split, dT_radial, 'bo', markersize = 5,alpha = 0.7)
plt.plot(np.arange(np.amin(void_r_avg_split) - 1, np.amax(void_r_avg_split) +1), np.arange(np.amin(void_r_avg_split) - 1, np.amax(void_r_avg_split) +1) * 0., 'm--',linewidth = 4, alpha = 0.3, label = "Zero Temperature Line")
plt.errorbar(void_r_avg_split, dT_radial, yerr=dT_radial_err, fmt='none', ecolor='red', elinewidth=2, markeredgewidth=2)
plt.xlim((np.amin(void_r_avg_split) - 0.5, np.amax(void_r_avg_split) +0.5))
plt.xlabel("Radial Bin [degrees]")
plt.ylabel(r'$\Delta$T [$\mu$K]')
plt.legend()
plt.savefig(plots_pwd+'temperature_profile.png')
plt.show()

