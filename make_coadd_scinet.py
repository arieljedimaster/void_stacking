import numpy as np
from flipper import *


folder = '/Users/Shared/ariel/Maps/bossN/' #scinet folder in renees directory
folder_dump = '/Users/Shared/ariel/Maps/bossN_coadd/' #scinet folder in my directory
#scinet maps
'''

maps = ['boss_s14_ar2_night_1way_0_sky_'] #only the one map

ways = ['map0001','map0002','map0005','map0010','map0020','map0050','map0100','map0200']
'''

###############

maps = ['boss_tot_ar2_night_tot_sky_']

ways = ['map0010']


for map in maps:
    #cmb = liteMap.liteMapFromFits(folder+map+ways[0]+'.fits')
    cmb = liteMap.liteMapFromFits(folder+'boss_tot_ar2_night_tot_sky_map0010.fits')
 
    cmb.data[:]=0.
    weights = cmb.copy(); weights.data[:]=0.;
    var = cmb.copy(); var.data[:]=0.;
   
    
    print('Coadding ', map)

    for set in ways: #iterating through the different sets
        print('Set ', set)
        
        #corr is the actual data
        corr   = liteMap.liteMapFromFits(folder+map+set+'.fits')
        print 'input map set ', set
        #gotta figure out which maps to input here 
        noise  = liteMap.liteMapFromFits(folder+map+'div.fits') #fill in the stars
        weight = 1./(noise)#liteMap.liteMapFromFits(folder+map+'**************.fits') #fill in the stars

        ii = np.where(noise.data[:]!=0.); #where there is data
        cmb.data[ii] +=(corr.data[ii])/((noise.data[ii])**2.)
        weights.data[ii] += weight.data[ii]
        var.data[ii] += 1./((noise.data[ii])**2.) #1/sigma^2 for all maps
         
		  
        jj = np.where(var.data[:]!=0) #where there is data
        var.data[jj] = 1./(var.data[jj]) 
        cmb.data[jj] *= var.data[jj]
        var.data[jj]=np.sqrt(var.data[jj])

        set_new='coadd'

        print(folder_dump+map+'map_'+set_new+'.fits')
        print(folder_dump+map+'weight_'+set_new+'.fits')

        cmb.writeFits(folder_dump+map+'map_'+set_new+'.fits', overWrite=True)
        weights.writeFits(folder_dump+map+'weight_'+set_new+'.fits', overWrite=True)
        var.writeFits(folder_dump+map+'var_'+set_new+'.fits', overWrite=True)
        
