import healpy as hp
import numpy as np
import pylab as pyl

from astropy.io import fits

import matplotlib.pyplot as plt

import void_healpy_functions as vhf

#this is adapted from tessas code

data_pwd = '/Users/arielamaral/Documents/void_stacking_healpy/Data/'
code_pwd = '/Users/arielamaral/Documents/void_stacking_healpy/Code/'
plots_pwd = '/Users/arielamaral/Documents/void_stacking_healpy/Plots/'


nside= 2048 # found from planck_nilc[1].header
npix = 12*nside**2
res = int(np.log(nside)/np.log(2))

#pixel indeces
pix_indeces = np.arange(npix)

#theta and phi projections fror planck map
planck_theta,planck_phi = hp.pix2ang(nside,pix_indeces)


#reading in my planck map
planck_nilc = fits.open(data_pwd+'COM_CompMap_YSZ_R2.00.fits/nilc_ymaps.fits')
planck_nilc_header = planck_nilc[1].header
planck_nilc_data = planck_nilc[1].data
planck_nilc_full = planck_nilc_data['FULL    ']
planck1 = planck_nilc_full


#load in void info
#Nadathur (2016) BOSS catalogue from: http://www.icg.port.ac.uk/stable/nadathur/voids/

void_RA, void_Dec, void_r_ang, void_ID, void_e_f = np.loadtxt(data_pwd+'void_info_2.txt', unpack=True)

#loading additional void info (from that same catalogue - can't remember why I put it into two different arrays)

void_ID, void_centre_x, void_centre_y, void_centre_z, void_min_dens, void_WtdAvgDens, void_r_eff, void_dens025, void_edge_flag, void_dens_ratio = np.loadtxt(data_pwd+'void_info_BOSS.txt', unpack = True)

'''
#only using the first 50 for now

void_RA = void_RA[:500]
void_Dec = void_Dec[:500]
void_r_ang = void_r_ang[:500]
'''

#nsource
nsources = len(void_RA)

void_theta, void_phi = vhf.DeclRaToThetaPhi(void_Dec,void_RA)


#split by radius:
radial_bins = 3 #??? can change this if I want

#sorting voids by radial bins
void_r_ang_split, void_phi_split, void_theta_split, void_r_avg_split= vhf.data_radial_binning('angular radius', void_r_ang, void_phi, void_theta, radial_bins)


#setting that I want 50 pixels resolution
#now we dont even have to resize!
nx=50
ny=50


#all the cutouts of some shape nx ny, because I set this.
#ims=np.zeros((nsources,nx,ny))

stack_per_bin = np.zeros((radial_bins,nx,ny))

num_voids_in_bin = np.zeros(radial_bins)

#setting up deltaT arrays

dT_radial = np.zeros(radial_bins)
dT_radial_err = np.zeros(radial_bins)


#looping through each individual void which falls in the map boundary

for j in np.arange(0, radial_bins):
	print "On void radial bin # " + str(j+1)
	void_phi_bin = void_phi_split[j] #voids in that radial bin
	void_theta_bin = void_theta_split[j]
	void_r_ang_bin = void_r_ang_split[j]

	#getting the radius of the postage stamp in radians
	radius_void_deg = 0.7*void_r_avg_split[j]
	radius_deg = radius_void_deg*3.0 #radius of postage stamp
	radius_rad = radius_deg*(np.pi/180.)

	#diameter (or total side length) of the postage stamp.
	#how many radians per pixel.
	drx=radius_rad*2/nx
	dry=radius_rad*2/ny

	#locations in radians of each of the pixels in the stamp in x and y
	xd=np.arange(0,radius_rad*2.,drx)-radius_rad
	yd=np.arange(0,radius_rad*2.,dry)-radius_rad

	#putting into a 2D array
	xxd,yyd=np.meshgrid(xd,yd)

	num_voids_in_bin_j = 0

	for i in np.arange(0,len(void_phi_bin)):
		if (i % 100 == 0) & (i!=0): 
		    print ">>> Stacking Void # " + str(i) + ", we are " + str(float(i)/float(len(void_phi_bin))*100.) + " percent done for this bin."
		
		#coordinates in healpy radians
		src_phi=void_phi_bin[i]
		src_theta=void_theta_bin[i]

	    #dimension of stamp in theta and phi
		phi_dist=xxd+src_phi
		theta_dist=yyd+src_theta

		#gets all pixel values within the distances specified
		vals_stamp=hp.get_interp_val(planck_nilc_full,theta_dist.flatten(),phi_dist.flatten())
		pix_annulus, vals_annulus = vhf.get_annulus_galactic(src_phi,src_theta,radius_void_deg,planck_nilc_full,nside)
		#taking the average of the annulus
		annulus_avg = np.average(vals_annulus)
		mod_vals_stamp = vals_stamp - annulus_avg

		#calculate individual dT from this
		radius_void_rad = (np.pi/180.)*radius_void_deg

		#need to find out how many pixels is in 
		dr_void=radius_void_rad*2/nx #nx is from above, should be the same


		#reshaping into 2D dimensions
		stack_per_bin[j] += mod_vals_stamp.reshape(nx,ny)
		num_voids_in_bin_j +=1

	num_voids_in_bin[j] = num_voids_in_bin_j
	stack_per_bin[j] = stack_per_bin[j]/float(num_voids_in_bin_j)

#total stack
void_stack_total = np.average(stack_per_bin, axis = 0)

plt.figure()
plt.title("Stacked Planck NILC Map for " + str(nsources) + " Voids")
plt.imshow(void_stack_total)
plt.colorbar()
plt.savefig(plots_pwd + "stacked_visual_plot_allvoids.png")
plt.show()

#plotting the stacked voids per radial bin:

for r in np.arange(0,radial_bins):

	plt.figure()
	plt.title("Stacked Void map for Radial Bin " + str(r + 1) + " with " + str(num_voids_in_bin[r]) + " voids")
	plt.imshow(stack_per_bin[r])
	plt.colorbar()
	plt.savefig(plots_pwd + "stacked_healpy_bin_"+str(r+1)+".png")
	plt.show()


'''
plt.figure()
plt.title(r"Void $\Delta$T profile for " + str(radial_bins) + " radial bins")
plt.plot(void_r_avg_split, dT_radial, 'bo', markersize = 5,alpha = 0.7)
plt.plot(np.arange(np.amin(void_r_avg_split) - 1, np.amax(void_r_avg_split) +1), np.arange(np.amin(void_r_avg_split) - 1, np.amax(void_r_avg_split) +1) * 0., 'm--',linewidth = 4, alpha = 0.3, label = "Zero Temperature Line")
plt.errorbar(void_r_avg_split, dT_radial, yerr=dT_radial_err, fmt='none', ecolor='red', elinewidth=2, markeredgewidth=2)
plt.xlim((np.amin(void_r_avg_split) - 0.5, np.amax(void_r_avg_split) +0.5))
plt.xlabel("Radial Bin [degrees]")
plt.ylabel(r'$\Delta$T [$\mu$K]')
plt.legend()
plt.savefig(plots_pwd+'temperature_profile_from_image.png')
plt.show()
'''

