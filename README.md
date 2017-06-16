# void_stacking

The most recent version of the "main" python file is: main_faster.py 
Use main_faster.py to run code, this one will not take up too much memory on the computer as it does not save the individual CMB postage stamps of each void. It also contains update plots (better looking plot code), and a fix to a few mistakes made while calculating deltaT in the other code.
One thing this code does not do it calculating the stacked void maps for random locations, but can be easily added in using the same way in main.py (literally just copy and paste!), the code which generates the random void locations and sizes is there, but it doesnt stack or anything.

Before running:
1) make sure map_dump is going to the directory in which you want to output your files
2) If you want to use different CMB maps, you must upload, and crop them in the mannar done in the code
3) To upload new void maps, do so in the same way done in the code (follow the same method)
4) If you want to change the /ell modes the code filters out of the CMB maps go to the generate_filtsmap(smap)function in void_filt.py and change the filter function! (instructions are in that function)

All of these are in lines 1 - 138 in main_faster.py (these lines literally just correspond to loading in the maps and catalogues, and cropping them.


The code contains the following python files:

1) main.py: old main file, will run out of memory if using large maps/catalogues, dont use!
2) main_faster.py: USE THIS MAIN FILE
3) main_new.py: older main file which used newer catalogues and maps (but still uses the old stacking method) dont use!
4) make_coadd_scinet.py: code to coadd maps, haven't used it in a while
5) map_plot.py: random plotting code, dont use
6) void_ACT_upload_crop.py: functions which uploads ACT maps (deep5,6,56), crops them
7) void_cat_load.py: function which loads void catalogue info from a fits file (only for the Clampett void catalogue)
8) void_filt.py: all functions which involve filtering/masking/annulus generating 
9) void_locations.py: uploads and puts together void info from an old catalogue, outdated dont use!
10) void_plots.py: all functions which involve plotting, or converting from pixels to degrees for plotting purposes
11) void_reduc.py: all functions which involve cutting out voids which do not fit in map boundaries, sorting voids into radial bins generating random void locations etc.
12) void_stack.py: all functions which involve physically stacking voids, resizing stacked void stamps to get a final resized void map, void temperature profile, deltaT statistics, signal to noise stuff (which wasnt used)

To run the code:

ipython: "run main_faster.py -"
or
"python2.7 main_faster.py"
That's it! you just run the main_faster python file and it'll do the rest!
The only thing you'll ever have to change is the first few lines of code in this file which is outlined above


