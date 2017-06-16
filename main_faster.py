import pyfits
import numpy as np
import os

#import matplotlib.pyplot as plt
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import scipy
from scipy import stats
#from astLib import astWCS
import random as random
#from astLib import astCoords

from flipper import *

# to time how long each step takes
import time


#importing functions from other files
from void_reduc import *
from void_filt import *
from void_stack import *
from void_plots import *
from void_ACT_upload_crop import *
from void_cat_load import *

maps_pwd = '/Users/shared/ariel/Maps/Coadded/'

code_pwd = '/Users/arielamaral/Code/void_stacking/'

### !!! CHANGE THIS GIVEN WHERE YOU WANT TO PUT ALL THE OUTPUT PLOTS:
#right now this is to the folder "old_map_old_cat_new_code", because I'm using older ACT maps and catalogues, but running it with my faster and newer code

map_dump = '/Users/arielamaral/Code/void_stacking/old_map_old_cat_new_code/'


#UPLOADING ACT CMB MAPS


act_map_deep_5, act_map_deep_6, act_map_deep56_AR2, act_map_deep56_AR1 = ACT_upload()


#splitting up Deep5 and Deep56 map using submaps, because flipper has trouble when the RA coordinates wrap around

act_map_deep56_AR2_1, act_map_deep56_AR2_2, act_map_deep56_AR1_1, act_map_deep56_AR1_2, act_map_deep5_1, act_map_deep5_2 = ACT_split(act_map_deep_5, act_map_deep_6, act_map_deep56_AR2, act_map_deep56_AR1)


#CROPPING ACT MAPS TO GET RID OF NOISY EDGES
#the dimensions on where to crop them was determined by eye

act_map_deep5_2 = act_map_deep5_2.selectSubMap(350.0, 360.0 -0.00005, -4, 4)

act_map_deep5_1 = act_map_deep5_1.selectSubMap(0.0, 4.0, -4.0, 4.0)

act_map_deep_6 = act_map_deep_6.selectSubMap(28.0, 41.0, -8.0, -1.0)

act_map_deep56_AR1_1 = act_map_deep56_AR1_1.selectSubMap(0.0, 45.0, -8.0, 5.0)

act_map_deep56_AR1_2 = act_map_deep56_AR1_2.selectSubMap(348.0, 360.0 -0.00005, -8.0, 5.0)

act_map_deep56_AR2_1 = act_map_deep56_AR2_1.selectSubMap(0.0, 45.0, -8.0, 5.0)

act_map_deep56_AR2_2 = act_map_deep56_AR2_2.selectSubMap(348.0, 360.0 -0.00005, -8.0, 5.0)

act_map_deep6 = act_map_deep_6

#we'll use AR2 instead of AR1 for the deep56 maps

maps = [act_map_deep5_2, act_map_deep6, act_map_deep5_1, act_map_deep56_AR1_1, act_map_deep56_AR1_2]


map_names = ["act_map_deep5_2", "act_map_deep6", "act_map_deep5_1", "act_map_deep56_AR1_1", "act_map_deep56_AR1_2"]


map_dict = dict(zip(maps,map_names))


'''

maps_pwd = '/Users/arielamaral/Code/extracted_BossN/extracted_BossN/' 


#map_dump = '/Users/arielamaral/Code/void_code/new_figures'

map_dump = '/Users/arielamaral/Code/void_stacking/new_map_new_cat_new_code/'


#UPLOADING NEWER MAP

act_map0010_raw = liteMap.liteMapFromFits(maps_pwd+'extracted_boss_tot_ar2_night_tot_sky_map0010.fits')

#inverse variance weight map
act_div = liteMap.liteMapFromFits(maps_pwd+'extracted_boss_tot_ar2_night_tot_sky_div.fits')

#inverse variance weight

act_map0010 = act_map0010_raw.copy()

act_map0010.data = act_map0010_raw.data*act_div.data


#Cropping newer ACT Map
act_map0010 = act_map0010.selectSubMap(125, 240, -5, 20)

maps = [act_map0010] #from scinet

map_names = ['map0010'] #fill this in

map_dict = dict(zip(maps,map_names))

'''

#void_RA, void_Dec, void_r_eff, void_cent, void_vol, void_r_eff_2, void_ID, void_dens = np.loadtxt('void_info.txt', unpack=True)

#pixel scale
delta_xi = 0.494350/60 #degrees
delta_yi = 0.494986/60 #degrees


#load in void info
#Nadathur (2016) BOSS catalogue from: http://www.icg.port.ac.uk/stable/nadathur/voids/

void_RA, void_Dec, void_r_ang, void_ID, void_e_f = np.loadtxt('void_info_2.txt', unpack=True)

#loading additional void info (from that same catalogue - can't remember why I put it into two different arrays)

void_ID, void_centre_x, void_centre_y, void_centre_z, void_min_dens, void_WtdAvgDens, void_r_eff, void_dens025, void_edge_flag, void_dens_ratio = np.loadtxt('void_info_BOSS.txt', unpack = True)

#larger void catalogue:
#catalogue from Joseph Clampitt, which was given to me by an email from Renee:

'''
#load in void info from fits files
void_RA, void_Dec, void_z, void_rlos, void_r_eff, void_r_ang, void_los_size, void_dens, void_f_vol = get_void_info_from_fits()
void_r_ang = arcmin_to_degrees(void_r_ang)
'''

#change all void_RA to positive ones
void_RA = neg_void_RA(void_RA)

print "Total length of the catalogue: ", len(void_RA)

print "The void effective radius range of the catalogue is: ", np.amin(void_r_eff), np.amax(void_r_eff)


#number of radial bins

n = 3

#if the first argument in the function radialbin() is 'effective radius', it means we'll bin by effective radius
#if the argument is 'angular radius', we'll bin by angular radius
#this splits up the voids into n lists, index containing a list of the void details for that radial bin!

void_r_ang, void_RA, void_Dec, void_r_avg, void_r_eff, void_r_eff_avg = radial_bin('angular radius', void_r_ang, void_RA, void_Dec, void_r_eff, n)


print "void_r_eff_avg before loop: ", void_r_eff_avg


#for random locations use the same void radii
rand_r_ang = void_r_ang
rand_r_avg = void_r_avg

#keeping track of the delta T lists of all the voids in each bin
delt_T_for_bin = [[] for _ in range(n)]

#all stacked maps
#this is so that I have them all in an array ready to rescale and restack into a single graph

all_stacked_maps = []

all_stacked_maps_mod = []

void_bin_temp_rescale = []

#a list of this is also necessary for rescaling stacked voids later on in the code

R_filt_list = []

total_voids_num = 0

ref_map = maps[0]

#looping through each radial bin
for j in np.arange(0,len(void_RA)):
    print "                    "
    print "############# Stacking for Bin ", j
    print "                    "
    
    num_in_stack = 0 #starting the count
    
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
    
    t1=time.time()
    
    stacked_mod_smap = 'hello'
    stacked_smap = 'hello'
    
    #now looping through each of the CMB maps so we can stack
    
    for orig_map in maps:
        print "..... Map ", map_dict[orig_map]

        #filter the whole map
        #check to make sure that the filter is what you want!
        map = generate_filtsmap(orig_map)
        
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
        
        #looping through each individual void which falls in the map boundary
        
        if len(void_RA_new) !=0:
        	if num_in_stack == 0:
        		print "first map"
        		smap = map.selectSubMap(void_RA_new[0]-bound_x, void_RA_new[0] + bound_x, void_Dec_new[0] - bound_y, void_Dec_new[0] + bound_y)
        		smap_data = smap.data[:]
        		stacked_mod_smap = np.zeros(smap_data.shape) #will need to set shape I think
        		stacked_smap = np.zeros(smap_data.shape) #will neep to set shape I think
        		void_bin_temp = 0.
        	
        	for k in np.arange(0,len(void_RA_new)):
				num_in_stack +=1
				total_voids_num += 1
				smap = map.selectSubMap(void_RA_new[k]-bound_x, void_RA_new[k] + bound_x, void_Dec_new[k] - bound_y, void_Dec_new[k] + bound_y)
				annulus, inner_circ = delt_T_filt(void_RA_new[k], void_Dec_new[k], R, smap)
				mod_smap = mod_filt(smap, annulus)
				smap_data = smap.data[:]
				mod_smap_data = mod_smap.data[:]
				arrays = [smap_data,stacked_smap]
				dim_list = min_dim_of_list_of_arrays(arrays)
				#cropping the stacked map
				stacked_smap = crop_array(stacked_smap, dim_list)
				stacked_mod_smap = crop_array(stacked_mod_smap, dim_list)
				#cropping the smaps
				smap_data = crop_array(smap_data, dim_list)
				mod_smap_data = crop_array(mod_smap_data, dim_list)
				stacked_mod_smap += mod_smap_data
				stacked_smap += smap_data
				inner_circ_data = inner_circ.data[:]
				annulus_data = annulus.data[:]
				void_bin_temp += (np.mean(inner_circ_data[inner_circ_data != 0.]) - np.mean(annulus_data[annulus_data != 0.]))
        else:
            continue
        
    if num_in_stack != 0:
        #dividing by counter
        print 'stacking all maps for bin: ', j
        stacked_mod_smap *= 1./float(num_in_stack)
        stacked_smap *= 1./float(num_in_stack)
        void_bin_temp *= 1./float(num_in_stack)
        void_bin_temp_rescale.append(void_bin_temp)
        
        all_stacked_maps.append(stacked_smap)
        all_stacked_maps_mod.append(stacked_mod_smap)
        
        
        #plotting the stacked maps per bin!
        
        void_cent_x = stacked_smap.shape[0]/2.0
    	void_cent_y = stacked_smap.shape[1]/2.0
   		 
    	fig,ax = plt.subplots(1)
    	ax.set_aspect('equal')
    	ax.imshow(stacked_smap, origin = 'lower')
    	circle1 = Circle((void_cent_x,void_cent_y), radius=degrees_to_pix(R, ref_map), color='k', linewidth = 2, fill=False)
    	circle2 = Circle((void_cent_x,void_cent_y), radius=np.sqrt(2)*degrees_to_pix(R, ref_map), color='k', linewidth = 2, fill=False)
    	ax.text(stacked_smap.shape[0]*(1./10.), stacked_smap.shape[1]*(1./10.), r'$\Delta$T = ' +str(void_bin_temp) + r' $\mu$K', style='italic',bbox={'facecolor':'white', 'alpha':0.5, 'pad':10})
    	
    	ax.add_patch(circle1)
    	ax.add_patch(circle2)
    	cax = ax.imshow(stacked_smap, origin = 'lower')
    	cbar = fig.colorbar(cax)
    	fig.canvas.draw()
    	image_size_x = stacked_smap.shape[1]
    	image_size_y = stacked_smap.shape[0]
    	num_ticks = 7
        xlabels, ylabels = pix_to_deg_ticks(stacked_smap,ref_map)
        ticks_num_x = np.linspace(0,len(xlabels)-1, num_ticks).astype(int)
        ticks_num_y = np.linspace(0,len(ylabels)-1, num_ticks).astype(int)
        xlabels = xlabels[ticks_num_x]
        ylabels = ylabels[ticks_num_y]
        xlabels = np.around(xlabels, decimals=1)
        ylabels = np.around(ylabels, decimals=1)
        plt.xticks(ticks_num_x, xlabels, rotation='horizontal')
        plt.yticks(ticks_num_y, ylabels, rotation='vertical')
    	plt.title("Stacked Void Locations for Bin " + str(j)+ " for " + str(num_in_stack) + " voids")
    	plt.xlabel("[Degrees]")
    	plt.ylabel("[Degrees]")
    	plt.savefig(map_dump+'stacked_bin_'+str(j))
           
           
    else:
        continue
    
    
    #total number of voids in the bin stack
    print 'The total number of voids stacked in bin '+ str(j) + ' is ' + str(num_in_stack)
    t2=time.time()
    print "Stacking for this bin took " + str((t2-t1)/60.) + " minutes."


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

#I'm doing the rescale only for the *NON* modified (annulus subtracted) bins
# If you want to change that just change "all_stacked_maps" to "all_stacked_maps_mod" below

for z in np.arange(0, len(all_stacked_maps)):
	#rescaling for individual radial bin 'z'
	#we're rescaling all maps to the smallest sized map
    rescaled_map = stack_rescale(all_stacked_maps[z], R_filt_list[z], R_filt_ref)
    rescaled_stacked_maps.append(rescaled_map)
    #plotting each of the individual rescaled stacked maps:
    rescaled_stacked_map_per_bin_plot(all_stacked_maps[z],rescaled_map,z,R_filt_list[z], R_filt_ref,map_dump, ref_map)

#now we will put together all of the maps into one rescaled stack
total_rescaled_stack =  stack_all_rescaled_maps(rescaled_stacked_maps) 

#getting delta_T for the new rescaled map

#computing the "outer circle" of the rescaled void data array
outer_circ = array_mask(total_rescaled_stack, np.sqrt(2)*R_filt_ref, ref_map)

#computing the "inner circle" of the rescaled void data array
inner_circ = array_mask(total_rescaled_stack, R_filt_ref, ref_map)

#computing the annulus by subtracting the two (to get everything between the inner and the outer circle
annulus =  outer_circ - inner_circ

#annulus subtracting the inner circle
void_bin_temp_rescale = inner_circ - annulus


print "Rescaled delta T: ", void_bin_temp_rescale


total_rescaled_stacked_map_plot2(total_rescaled_stack, R_filt_ref, void_bin_temp_rescale,map_dump, total_voids_num, ref_map)


print "DONE CHECK PLOTS FOLDER"

