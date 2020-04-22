
"""
This Fiji Jython script will convert a two channel .oir xxx file into 2 (two) channel .tif files

Run this inside Fiji by selecting menu:
	Plugins - Macro - Run...

Assumes:
	The file you select is a proper .oir and has two channels

Caveats:
	1) This assumes all paths including folder names and file names do NOT have spaces !!!
	2) if you have a finder window open, the new files will not immediately appear
	3) This will overwite _ch1/_ch2 .tif file if they already exist
	
Recipe:
=======
run("Bio-Formats Importer", "open=/Users/cudmore/box/data/sami/Cell_1/1_5ADVMLEG1L1.oir color_mode=Default rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT");
run("8-bit");
run("Split Channels");
selectWindow("C1-1_5ADVMLEG1L1.oir");
run("Save", "save=/Users/cudmore/box/data/sami/Cell_1/c12_ch1.tif");
selectWindow("C2-1_5ADVMLEG1L1.oir");
run("Save", "save=/Users/cudmore/box/data/sami/Cell_1/c12_ch2.tif");

run("Z Project...", "projection=[Max Intensity]");
selectWindow("8_5ADVMLEG1L1_ch2_dvSkel_1.tif");
run("Z Project...", "projection=[Max Intensity]");
run("Merge Channels...", "c1=MAX_8_5ADVMLEG1L1_ch2_dvSkel_1.tif c2=MAX_8_5ADVMLEG1L1_ch2.tif create keep");
selectWindow("8_5ADVMLEG1L1_ch2_dvSkel_1.tif");
run("Merge Channels...", "c1=8_5ADVMLEG1L1_ch2_dvSkel_1.tif c2=8_5ADVMLEG1L1_ch2.tif create keep");
"""

# has to be the first line
from __future__ import print_function

import os

from ij import IJ, WindowManager
from ij.io import OpenDialog

#from loci.plugins import BF

#path = '/Users/cudmore/box/data/sami/Cell_1/1_5ADVMLEG1L1.oir'
path = '/Users/cudmore/box/data/sami/200108/WT_Female/Cell_8/8_5ADVMLEG1L1_ch2.tif'

# ask user for file. I do not know how to handle when users hits cancel??? Script will just fail
if path is None or not os.path.isfile(path):
	notUsedBy_oxs = 'Open a _ch2.tif file'
	path = OpenDialog(notUsedBy_oxs).getPath()
	if path is None:
		exit()
	
print('    user selected path:', path)
fileName = os.path.basename(path)

# load
imp = IJ.openImage(path)
imp.show()

# save channel 1
windowName = 'C1-' + fileName
IJ.selectWindow(windowName)
windowName_notSure = ij.WindowManager.getImage(windowName)

saveFilePath = os.path.splitext(path)[0] + '_ch1.tif' # this is getting garbled when we have spaces (' ') in path !!!!!!!!!
print('    saving', saveFilePath)
#windowName_notSure.close()

print('    done')

