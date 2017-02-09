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
from matplotlib.patches import Circle

from flipper import *

#importing functions from other files
from void_reduc import *
from void_filt import *
from void_stack import *


maps_pwd = '/Users/shared/ariel/Maps/Coadded/'

code_pwd = '/Users/arielamaral/Code/void_code/'

map_dump = '/Users/arielamaral/Code/void_code/old_map_old_cat/'

deep5_map_dump = code_pwd+'deep5_plots/' #all plots relating to deep 5 go here

deep6_map_dump = code_pwd+'deep6_plots/' #all plots relating to deep 6 go here

deep56_map_dump = code_pwd+'deep56_plots/' #all plots relating to deep 56 go here

############################################################
# All the plots here are PER RADIAL BIN, so you need to input the bin number
# These plots go inside the loop while looping through the various bins in main.py

def degrees_to_pix(r_deg):
    r_arcmin = 60*r_deg #degrees to arcmins
    pix_scale_x = 0.00014532290052643554 # for act_map_deep5_2.pixScaleX
    pix_scale = (pix_scale_x)*((180/np.pi) *(60)) #I chose a random specific map to use
    r_pix = r_arcmin / pix_scale 
    return r_pix

def void_stack_bin_plot(stacked_voids_bin,R,num_in_stack,j,dT_void_stack_in_bin):
    #saving a figure of the stacked voids
    
    #plt.figure()
    #plt.title('Stacked Void CMB Locations in Bin ' + str(j) + ' with '+ str(num_in_stack)+ ' voids')
        
    void_cent_x = stacked_voids_bin.shape[0]/2.0
    void_cent_y = stacked_voids_bin.shape[1]/2.0
    #plt.plot(void_cent_x,void_cent_y, 'wo')
    
    print "R in pixel scale", degrees_to_pix(R) 
    print "Outer annulus in pixel scales", np.sqrt(2)*degrees_to_pix(R)
        
    #stacked_voids_bin.plot()
    #yes = plt.imshow(stacked_voids_bin, origin = 'lower')
        
    #printing apparent centre of map
    #print 'The centre of the map is (x,y) = ('+str(void_cent_x)+','+str(void_cent_y)+')'
    print "The matrix size of this bin is: ", stacked_voids_bin.shape
    print 'The centre of the map is (x,y) = (' + str(stacked_voids_bin.shape[0]/2.0) + ',' + str(stacked_voids_bin.shape[1]/2.0) + ')'
        
        
    #circle1 = plt.Circle((void_cent_x,void_cent_y), radius=degrees_to_pix(R), color='k', linewidth = 2, fill=False)
    #circle2 = plt.Circle((void_cent_x,void_cent_y), radius=np.sqrt(2)*degrees_to_pix(R), color='k', linewidth = 2, fill=False)
    x_shape = stacked_voids_bin.shape[0]
    y_shape = stacked_voids_bin.shape[1]
    
    fig,ax = plt.subplots(1)
    ax.set_aspect('equal')
    ax.imshow(stacked_voids_bin, origin = 'lower')
    #ax.plot(void_cent_x,void_cent_y, 'wo')
    circle1 = Circle((void_cent_x,void_cent_y), radius=degrees_to_pix(R), color='k', linewidth = 2, fill=False)
    circle2 = Circle((void_cent_x,void_cent_y), radius=np.sqrt(2)*degrees_to_pix(R), color='k', linewidth = 2, fill=False)
    
    ax.add_patch(circle1)
    ax.add_patch(circle2)
    ax.text(stacked_voids_bin.shape[0]*(1./10.), stacked_voids_bin.shape[1]*(1./10.), 'deltaT = ' +str(dT_void_stack_in_bin), style='italic',
        bbox={'facecolor':'white', 'alpha':0.5, 'pad':10})
    plt.title('Stacked Void CMB Locations in Bin ' + str(j) + ' with '+ str(num_in_stack)+ ' voids')
    plt.xlim(0,x_shape)
    plt.ylim(0,y_shape)
        
    plt.savefig(map_dump+'void_stack_all_maps_bin_' + str(j)+'.png')

        
    #total number of voids in the bin stack
    print 'The total number of random locations stacked in bin '+ str(j) + ' is ' + str(num_in_stack)
    return

def rand_stack_bin_plot(stacked_rand_bin,R,num_in_stack,j,dT_rand_stack_in_bin):
    #saving a figure of the stacked voids

    print "R in pixel scale", degrees_to_pix(R) 
    print "Outer annulus in pixel scales", np.sqrt(2)*degrees_to_pix(R)
        
    rand_cent_x = stacked_rand_bin.shape[0]/2.
    rand_cent_y = stacked_rand_bin.shape[1]/2.
    
    x_shape = stacked_rand_bin.shape[0]
    y_shape = stacked_rand_bin.shape[1]
    
    print "The matrix size of this bin is: ", stacked_rand_bin.shape
    print 'The centre of the map in pixels is (x,y) = (' + str(stacked_rand_bin.shape[0]/2.) + ',' + str(stacked_rand_bin.shape[1]/2.) + ')'
    
    fig,ax = plt.subplots(1)
    ax.set_aspect('equal')
    ax.imshow(stacked_rand_bin, origin = 'lower')
    #ax.plot(rand_cent_x,rand_cent_y, 'wo')
    circle1 = Circle((rand_cent_x,rand_cent_y), radius=degrees_to_pix(R), color='k', linewidth = 2, fill=False)
    circle2 = Circle((rand_cent_x,rand_cent_y), radius=np.sqrt(2)*degrees_to_pix(R), color='k', linewidth = 2, fill=False)
    
    ax.add_patch(circle1)
    ax.add_patch(circle2)
    
    ax.text(stacked_rand_bin.shape[0]*(1./10.), stacked_rand_bin.shape[1]*(1./10.), 'deltaT = ' +str(dT_rand_stack_in_bin), style='italic',
        bbox={'facecolor':'white', 'alpha':0.5, 'pad':10})
    
    plt.xlim(0,x_shape)
    plt.ylim(0,y_shape)
    plt.title('Stacked Random CMB Locations for Bin '+ str(j)+ ' with '+ str(num_in_stack)+ ' locations')
    
        
    plt.savefig(map_dump+'rand_stack_all_maps_bin_' + str(j)+'.png')
    return


def delta_T_hist_bin_plot(del_T_list,delT_mean,del_T_list_rand,delT_mean_rand,j):
    plt.figure()
    yrange = np.zeros(300)
    plt.hist(del_T_list, bins=20, histtype='stepfilled', color ='b',alpha=0.5, label = 'Void Filtered Temp')
    plt.hist(del_T_list_rand, histtype='stepfilled', color = 'r',alpha=0.5, label='Random Filtered Temp')
    plt.plot(yrange + delT_mean, np.arange(0, len(yrange)),'b--',linewidth = 2)
    plt.plot(yrange + delT_mean_rand, np.arange(0, len(yrange)),'r--',linewidth = 2)
    plt.title("Histogram of Filtered Void Tempertures for Bin "+ str(j))
    plt.xlabel('Filtered Void Temperature, Delta_T')
    plt.ylabel("Frequency")
    plt.legend()
    plt.savefig(map_dump+'filtered_temp_hist_bin_'+str(j)+'_filt_10_300.png')

    return

def delta_T_vs_r_eff_bin_plot(del_T_list, void_r_eff_in_bin, j):            
    #plot of delta_T versus void radius
    #recreating Figure 2 from Cai et al. (2016) for each bin
    plt.figure()
    plt.title("Void Effective Radius vs. Temperature")
    plt.ylabel("Void Effective Radius [Mpc/h]")
    plt.xlabel("Temperature [micro-K]")
    plt.plot(del_T_list, void_r_eff_in_bin, 'ro')
    plt.savefig(map_dump+'delT_vs_r_eff_bin_'+str(j)+'.png')

    return

#######################################################################################
# These plots are for ALL of the voids, they do not depend or radial bins, or maps
# They go outside all of the loops in main.py

def rescaled_stacked_map_per_bin_plot(stacked_map_bin,rescaled_map,z,R_filt,R_filt_ref):
    '''
    Plots both the rescaled map and the original map to see the pixel number difference
    all_stacked_maps[z] = stack_map_bin
    '''
    
    x_original = stacked_map_bin.shape[0]
    y_original = stacked_map_bin.shape[1]
    x_rescaled = rescaled_map.shape[0]
    y_rescaled = rescaled_map.shape[1]
    
    rescale_param = R_filt_ref/R_filt
    
    void_cent_x = stacked_map_bin.shape[0]/2.
    void_cent_y = stacked_map_bin.shape[1]/2.
    
    rescaled_void_cent_x = rescaled_map.shape[0]/2.
    rescaled_void_cent_y = rescaled_map.shape[1]/2.
    
    fig = plt.figure()
    ax = fig.add_subplot(211)
    
    #original stacked map
    
    ax.imshow(stacked_map_bin, origin = 'lower')
    #circ= Circle((1,3), 1,color='k', linewidth = 2, fill=False)
    circle1 = Circle((void_cent_x,void_cent_y), radius=degrees_to_pix(R_filt), color='k', linewidth = 2, fill=False)
    circle2 = Circle((void_cent_x,void_cent_y), radius=np.sqrt(2)*degrees_to_pix(R_filt), color='k', linewidth = 2, fill=False)
    ax.add_patch(circle1)
    ax.add_patch(circle2)
    #ax.plot(void_cent_x,void_cent_y,'wo')
    plt.xlim(0,x_original)
    plt.ylim(0,y_original)
    plt.title("Original Void Stacked Map, Bin "+str(z))
    
    #rescaled stacked map
    
    ax1 = fig.add_subplot(212)
    ax1.imshow(rescaled_map, origin = 'lower')
    circ= Circle((2,3), 1,color='k', linewidth = 2, fill=False)
    circle1 = Circle((rescaled_void_cent_x,rescaled_void_cent_y), radius=rescale_param*degrees_to_pix(R_filt), color='k', linewidth = 2, fill=False)
    circle2 = Circle((rescaled_void_cent_x,rescaled_void_cent_y), radius=rescale_param*np.sqrt(2)*degrees_to_pix(R_filt), color='k', linewidth = 2, fill=False)
    ax1.add_patch(circle1)
    ax1.add_patch(circle2)
    #ax1.plot(rescaled_void_cent_x,rescaled_void_cent_y,'wo')
    plt.xlim(0,x_rescaled)
    plt.ylim(0,y_rescaled)
    plt.title("Rescaled Void Stacked Map, Bin "+str(z))
    fig.savefig(map_dump + "rescale_orig_bin_"+str(z)+"_void.png")
    return
    
def rescaled_stacked_map_per_bin_plot_rand(stacked_rand_bin,rescaled_map,z,R_filt,R_filt_ref):
    '''
    Plots both the rescaled map and the original map to see the pixel number difference
    all_stacked_maps[z] = stack_map_bin
    '''
    
    rescale_param = R_filt_ref/R_filt
    
    x_original = stacked_rand_bin.shape[0]
    y_original = stacked_rand_bin.shape[1]
    x_rescaled = rescaled_map.shape[0]
    y_rescaled = rescaled_map.shape[1]
    
    rand_cent_x = stacked_rand_bin.shape[0]/2.
    rand_cent_y = stacked_rand_bin.shape[1]/2.
    
    rescaled_rand_cent_x = rescaled_map.shape[0]/2.
    rescaled_rand_cent_y = rescaled_map.shape[1]/2.
    
    fig = plt.figure()
    ax = fig.add_subplot(211)
    
    #original stacked map
    
    ax.imshow(stacked_rand_bin, origin = 'lower')
    #circ= Circle((1,3), 1,color='k', linewidth = 2, fill=False)
    circle1 = Circle((rand_cent_x,rand_cent_y), radius=degrees_to_pix(R_filt), color='k', linewidth = 2, fill=False)
    circle2 = Circle((rand_cent_x,rand_cent_y), radius=np.sqrt(2)*degrees_to_pix(R_filt), color='k', linewidth = 2, fill=False)
    ax.add_patch(circle1)
    ax.add_patch(circle2)
    #ax.plot(rand_cent_x,rand_cent_y,'wo')
    plt.xlim(0,x_original)
    plt.ylim(0,y_original)
    plt.title("Original Random Stacked Map, Bin "+str(z))
    
    #rescaled stacked map
    
    ax1 = fig.add_subplot(212)
    ax1.imshow(rescaled_map, origin = 'lower')
    circ= Circle((2,3), 1,color='k', linewidth = 2, fill=False)
    circle1 = Circle((rescaled_rand_cent_x,rescaled_rand_cent_y), radius=rescale_param*degrees_to_pix(R_filt), color='k', linewidth = 2, fill=False)
    circle2 = Circle((rescaled_rand_cent_x,rescaled_rand_cent_y), radius=rescale_param*np.sqrt(2)*degrees_to_pix(R_filt), color='k', linewidth = 2, fill=False)
    ax1.add_patch(circle1)
    ax1.add_patch(circle2)
    #ax1.plot(rescaled_rand_cent_x,rescaled_rand_cent_y,'wo')
    plt.title("Rescaled Random Stacked Map, Bin "+str(z))
    plt.xlim(0,x_rescaled)
    plt.ylim(0,y_rescaled)
    fig.savefig(map_dump + "rescale_orig_bin_"+str(z)+"_rand.png")
    
    return

def total_rescaled_stacked_map_plot(total_rescaled_stack,R,rescaled_dT_void):
    void_cent_x = total_rescaled_stack.shape[0]/2.0
    void_cent_y = total_rescaled_stack.shape[1]/2.0
    
    fig,ax = plt.subplots(1)
    ax.set_aspect('equal')
    ax.imshow(total_rescaled_stack, origin = 'lower')
    circle1 = Circle((void_cent_x,void_cent_y), radius=degrees_to_pix(R), color='k', linewidth = 2, fill=False)
    circle2 = Circle((void_cent_x,void_cent_y), radius=np.sqrt(2)*degrees_to_pix(R), color='k', linewidth = 2, fill=False)
    
    ax.add_patch(circle1)
    ax.add_patch(circle2)
    
    ax.text(total_rescaled_stack.shape[0]*(1./10.), total_rescaled_stack.shape[1]*(1./10.), 'deltaT = ' +str(rescaled_dT_void), style='italic',
        bbox={'facecolor':'white', 'alpha':0.5, 'pad':10})
    
    plt.title("Rescaled Stacked Void Locations")
    plt.imshow(total_rescaled_stack, origin = 'lower')
    plt.savefig(map_dump+"rescaled_stack_voids.png")
    return
    
def total_rescaled_stacked_map_plot_rand(total_rescaled_stack, R,rescaled_dT_rand):

    void_cent_x = total_rescaled_stack.shape[0]/2.0
    void_cent_y = total_rescaled_stack.shape[1]/2.0
    fig,ax = plt.subplots(1)
    ax.set_aspect('equal')
    ax.imshow(total_rescaled_stack, origin = 'lower')
    circle1 = Circle((void_cent_x,void_cent_y), radius=degrees_to_pix(R), color='k', linewidth = 2, fill=False)
    circle2 = Circle((void_cent_x,void_cent_y), radius=np.sqrt(2)*degrees_to_pix(R), color='k', linewidth = 2, fill=False)
    
    ax.add_patch(circle1)
    ax.add_patch(circle2)
    
    ax.text(total_rescaled_stack.shape[0]*(1./10.), total_rescaled_stack.shape[1]*(1./10.), 'deltaT = ' +str(rescaled_dT_rand), style='italic',
        bbox={'facecolor':'white', 'alpha':0.5, 'pad':10})

    plt.title("Rescaled Stacked Random Locations")

    plt.savefig(map_dump+"rescaled_stack_random.png")

    return

def void_temp_profile_plot(void_r_eff_avg, avg_delT_bin, err_delT_bin_upper, err_delT_bin_lower,n):
    plt.figure()
    
    plt.plot(void_r_eff_avg, avg_delT_bin, 'ko-', linewidth = 2)
    plt.plot(void_r_eff_avg, err_delT_bin_upper, 'r--')
    plt.plot(void_r_eff_avg, err_delT_bin_lower, 'r--')

    plt.title("Void Temperature profile for "+str(n)+" bins")
    plt.xlabel("Void Effective Radius [Mpc/h]")
    plt.ylabel("Temperature [micro-K]")
    #plt.gca().invert_xaxis()
    plt.savefig(map_dump + "temp_void_profile_for_"+str(n)+"bins.png")
    return
    
def delT_significance_plot_bin(sig_list, j):
    plt.figure()
    #yrange = np.zeros(300)
    plt.hist(sig_list, bins=20, histtype='stepfilled', color ='b',alpha=0.5, label = 'Void Filtered Temp')
    #plt.hist(del_T_list_rand, histtype='stepfilled', color = 'r',alpha=0.5, label='Random Filtered Temp')
    #plt.plot(yrange + delT_mean, np.arange(0, len(yrange)),'b--',linewidth = 2)
    #plt.plot(yrange + delT_mean_rand, np.arange(0, len(yrange)),'r--',linewidth = 2)
    plt.title("Histogram of Delta(T)/sigma for bin "+ str(j))
    plt.xlabel('Significance of DeltaT [deltaT/sigma]')
    plt.ylabel("Frequency")
    plt.legend()
    plt.savefig(map_dump+'temp_significance_hist_bin_'+str(j)+'.png')
    return
