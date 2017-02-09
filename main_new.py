import pyfits
import numpy as np
import os

#import matplotlib.pyplot as plt
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
import scipy
from scipy import stats
from astLib import astWCS
import random as random
from astLib import astCoords

from flipper import *

#importing functions from other files
from void_reduc import *
from void_filt import *
from void_stack import *
from void_plots import *
from void_ACT_upload_crop import *
from void_cat_load import *

import sys


maps_pwd = '/Users/arielamaral/Code/void_code/extracted_BossN/' #fill in where the maps are on scinet

code_pwd = '/Users/arielamaral/Code/void_code/'

#map dump should be on the scatch pwd
#map_dump = '/Users/arielamaral/Code/void_code/new_figures'

map_dump = '/Users/arielamaral/Code/void_code/new_map_new_cat_figures/'



#upload new maps
#start with just uploading one of the maps

act_map0010_raw = liteMap.liteMapFromFits(maps_pwd+'extracted_boss_tot_ar2_night_tot_sky_map0010.fits')

#inverse variance weight map
act_div = liteMap.liteMapFromFits(maps_pwd+'extracted_boss_tot_ar2_night_tot_sky_div.fits')

#inverse variance weight


act_map0010 = act_map0010_raw.copy()

act_map0010.data = act_map0010_raw.data*act_div.data


#Split new maps
#crop new maps
act_map0010 = act_map0010.selectSubMap(125, 240, -5, 20)


maps = [act_map0010] #from scinet
'''

void_RA, void_Dec, void_r_eff, void_cent, void_vol, void_r_eff_2, void_ID, void_dens = np.loadtxt('void_info.txt', unpack=True)

'''

map_names = ['map0010'] #fill this in

map_dict = dict(zip(maps,map_names))

#pixel scale
delta_xi = 0.494350/60 #degrees
delta_yi = 0.494986/60 #degrees

'''
#load in void info

void_RA, void_Dec, void_r_ang, void_ID, void_e_f = np.loadtxt('void_info_2.txt', unpack=True)



#loading additional void info

void_ID, void_centre_x, void_centre_y, void_centre_z, void_min_dens, void_WtdAvgDens, void_r_eff, void_dens025, void_edge_flag, void_dens_ratio = np.loadtxt('void_info_BOSS.txt', unpack = True)
'''

#load in void info from fits files

void_RA, void_Dec, void_z, void_rlos, void_r_eff, void_r_ang, void_los_size, void_dens, void_f_vol = get_void_info_from_fits()


print "The void effective radius range of the catalogue is: ", np.amin(void_r_eff), np.amax(void_r_eff)

void_r_ang = arcmin_to_degrees(void_r_ang)

#change all void_RA to positive ones
void_RA = neg_void_RA(void_RA)

print "Total length of the catalogue: ", len(void_RA)


#input these void_RA void_Dec into radial_bin()
n = 3
# n = 1 means voids are stacked in one bin, no void radius dependence

#if the first argument in the function radialbin() is 'effective radius', it means we'll bin by effective radius
#if the argument is 'angular radius', we'll bin by angular radius

void_r_ang, void_RA, void_Dec, void_r_avg, void_r_eff, void_r_eff_avg = radial_bin('angular radius', void_r_ang, void_RA, void_Dec, void_r_eff, n)
# then loop through each bin in the for loop below

print "void_r_eff_avg before loop: ", void_r_eff_avg


#for random locations use the same (!) void radii
rand_r_ang = void_r_ang
rand_r_avg = void_r_avg

#keeping track of the delta T lists of all the voids in each bin
delt_T_for_bin = [[] for _ in range(n)]

#all stacked maps
#this is so that I have them all in an array ready to rescale and restack into a single graph

all_stacked_maps = []
all_stacked_maps_rand = []

#a list of this is also necessary for rescaling stacked voids later on in the code

R_filt_list = []


for j in np.arange(0,len(void_RA)):
	print "length of bin " + str(j) + " is: ", len(void_RA[j])

#sys.exit()

#looping through each bin
for j in np.arange(0,len(void_RA)):
    print "                    "
    print "############# Stacking for Bin ", j
    print "                    "
    
    num_in_stack = 0
    
    #double check void_r is in degrees
    # it looks like it is!
    
    R = 0.7 * void_r_avg[j]
    
    R_filt_list.append(R)
    
    #R = 0.5 # choosing constant radius right now
    
    print "The filter radius is ", R
    
    bound_x = 2.0 * R 
    bound_y = 2.0 * R
    
    #empty lists which contains the void stamp maps at different stages
    smaps = []
    filt_smaps = []
    annulus_list = []
    inner_circ_list = []
    mod_smaps = []
    
    #empty lists which contains the random stamp maps at different stages
    smaps_rand = []
    filt_smaps_rand = []
    annulus_list_rand = []
    inner_circ_list_rand = []
    mod_smaps_rand = []
    
    void_r_eff_in_bin =[]    
    

    #so we can stack all of the maps together
    #will loop through all of the maps individually
    for map in maps:
        print "..... Map ", map_dict[map]

        # filter the whole map
        map = generate_filtsmap(map)
        
        #get rid of void locations who are located less than the sub map bound
        #distance from the corner
        print("Total Number of Voids ", len(void_RA[j]))
        void_RA_new, void_Dec_new, void_r_ang_new, void_r_eff_new  = del_voids_close_to_submap_bounds(void_RA[j], void_Dec[j], void_r_ang[j], void_r_eff[j], bound_x, bound_y, map)
        print("Voids within map bounds ",len(void_RA_new))
        
        for a in np.arange(0,len(void_r_eff_new)):
            void_r_eff_in_bin.append(void_r_eff_new[a])
        
        #setting random locations with the same length as void_RA_new and void_dec_new
        #must take into account the sizes of patches
        
        rand_RA, rand_Dec = gen_rand_bounds(map, len(void_RA_new), bound_x, bound_y)
        
        #counting the number of voids going into this bin stack
        num_in_stack += len(void_RA_new)
        
        
        for k in np.arange(0,len(void_RA_new)):
            
            #stacking for void locations
            
            #generating a submap
            smap = map.selectSubMap(void_RA_new[k]-bound_x, void_RA_new[k] + bound_x, void_Dec_new[k] - bound_y, void_Dec_new[k] + bound_y)
            smaps.append(smap) #appending to giant list of all smaps
            annulus, inner_circ = delt_T_filt(void_RA_new[k], void_Dec_new[k], R, smap)
            annulus_list.append(annulus)#appending to giant list of all annuli
            inner_circ_list.append(inner_circ) #appending to giant list of all inner circles
            mod_smap = mod_filt(smap, annulus)
            mod_smaps.append(mod_smap) #appending to giant list of all modified smaps
            
            #stacking for random locations
            smap_rand = map.selectSubMap(rand_RA[k]-bound_x, rand_RA[k] + bound_x, rand_Dec[k] - bound_y, rand_Dec[k] + bound_y)
            smaps_rand.append(smap_rand) #appending to giant list of all smaps
            annulus_rand, inner_circ_rand = delt_T_filt(rand_RA[k], rand_Dec[k], R, smap_rand)
            annulus_list_rand.append(annulus_rand)#appending to giant list of all annuli
            inner_circ_list_rand.append(inner_circ_rand) #appending to giant list of all inner circles
            mod_smap_rand = mod_filt(smap_rand, annulus_rand)
            mod_smaps_rand.append(mod_smap_rand) #appending to giant list of all modified smaps
        
        
               
    #stacking all voids within this specific bin    
    stacked_voids_bin = stacked_voids(mod_smaps)
    stacked_rand_bin = stacked_voids(mod_smaps_rand)
    
    #appending stacked map to list of all stacked maps so I can rescale and stack them all on top of eachother
    
    all_stacked_maps.append(stacked_voids_bin)
    
    all_stacked_maps_rand.append(stacked_rand_bin)
    
    #total number of voids in the bin stack
    print 'The total number of voids stacked in bin '+ str(j) + ' is ' + str(num_in_stack)
    
    if mod_smaps == []:
        print "No plot, no voids to stack"
        
    else:
        
        #### LOTS OF PLOTTING SHIT NOW
        
        
        #saving a figure of the stacked voids
        dT_void_stack_in_bin = array_mask(stacked_voids_bin, R)
        void_stack_bin_plot(stacked_voids_bin,R,num_in_stack,j,dT_void_stack_in_bin)
        
        #####################################################################
        #plotting the stacked random locations
        dT_rand_stack_in_bin = array_mask(stacked_rand_bin, R)
        rand_stack_bin_plot(stacked_rand_bin,R,num_in_stack,j,dT_rand_stack_in_bin)
        

        
        #calculating average filtered temperatures for random and void locations
        del_T_list = del_T_measurement(inner_circ_list, annulus_list)
        delT_mean = np.mean(del_T_list)
        del_T_list_rand = del_T_measurement(inner_circ_list_rand, annulus_list_rand)
        delT_mean_rand = np.mean(del_T_list_rand)
        
        #appending the Delta T list to a list which will contain all of this info for all bins
        delt_T_for_bin[j].append(del_T_list)
        
        

        #generating a histogram of the temperature distribution in void and random locations
        delta_T_hist_bin_plot(del_T_list,delT_mean,del_T_list_rand,delT_mean_rand,j)
        #do a KS test here?
        
        #delta_T vs. R_eff plot
        delta_T_vs_r_eff_bin_plot(del_T_list, void_r_eff_in_bin, j)
        
        #temp_std = SN_std(mod_smaps, 'whole map', R)
        
        #sig_list = SN_significance(del_T_list, temp_std)
        
        #plot the significance histogram
        
        #delT_significance_plot_bin(sig_list, j)
        
        
        
        
        


#Now it's time to look at the list with all of the stacked maps
#and rescale them so that we can stack them all in one plot!

print "                            "

print "----------------------------"
print "Rescaling all binned stacked maps......"
print "----------------------------"

rescaled_stacked_maps = []
rescaled_stacked_maps_rand = []

#corresponds to the largest void bin
#because we want to re scale all smaller voids larger, that's just the way the function works

 
R_filt_ref = R_filt_list[0]

for z in np.arange(0, len(all_stacked_maps)):
    rescaled_map = stack_rescale(all_stacked_maps[z], R_filt_list[z], R_filt_ref)
    rescaled_stacked_maps.append(rescaled_map)
    #plotting each of the individual rescaled stacked maps:
    rescaled_stacked_map_per_bin_plot(all_stacked_maps[z],rescaled_map,z,R_filt_list[z], R_filt_ref)
    
    #lets do the same for random locations as well
    rescaled_map_rand = stack_rescale(all_stacked_maps_rand[z], R_filt_list[z], R_filt_ref)
    rescaled_stacked_maps_rand.append(rescaled_map_rand)
    #plotting the individual random rescaled maps
    rescaled_stacked_map_per_bin_plot_rand(all_stacked_maps_rand[z],rescaled_map_rand,z,R_filt_list[z], R_filt_ref)
    

#now we will put together all of the maps into one rescaled stack

total_rescaled_stack =  stack_all_rescaled_maps(rescaled_stacked_maps) 

total_rescaled_stack_rand = stack_all_rescaled_maps(rescaled_stacked_maps_rand)

#calculating the deltaT in the inner circle for the two rescaled stacked maps (void and random locations)

rescaled_dT_void = array_mask(total_rescaled_stack, R_filt_ref)

rescaled_dT_rand = array_mask(total_rescaled_stack_rand, R_filt_ref)

#now we can plot what this total rescaled stacked map looks like

total_rescaled_stacked_map_plot(total_rescaled_stack,R_filt_ref,rescaled_dT_void)

total_rescaled_stacked_map_plot_rand(total_rescaled_stack_rand, R_filt_ref,rescaled_dT_rand)

#calculating the void temperature profile from the various radial bins

avg_delT_bin, err_delT_bin_upper, err_delT_bin_lower = void_temperature_profile(delt_T_for_bin)


#plotting the void temperature profile
void_temp_profile_plot(void_r_eff_avg, avg_delT_bin, err_delT_bin_upper, err_delT_bin_lower,n)

