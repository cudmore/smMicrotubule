"""
docs:
    https://javadoc.scijava.org/Fiji/ij3d/Image3DUniverse.html
    
    https://javadoc.scijava.org/Fiji/ij3d/Content.html
    
"""
   
# has to be the first line
from __future__ import print_function

import os

from ij import IJ, WindowManager
from ij.io import OpenDialog

#from loci.plugins import BF

from ij3d import Image3DUniverse  

from org.scijava.vecmath import Color3f

stackPath = '/Users/cudmore/box/data/sami/191230/BIN1_smKO_Male/Cell_12/12_5ADVMLEG1L1_ch2.tif'
maskPath = '/Users/cudmore/box/data/sami/191230/BIN1_smKO_Male/Cell_12/12_5ADVMLEG1L1_ch2_dvMask_1.tif'
skelPath = '/Users/cudmore/box/data/sami/191230/BIN1_smKO_Male/Cell_12/12_5ADVMLEG1L1_ch2_dvSkel_1.tif'

# call("ij3d.ImageJ3DViewer.add", "12_5ADVMLEG1L1_ch2.tif", "Green", "12_5ADVMLEG1L1_ch2.tif-2", "0", "true", "true", "true", "2", "0");

imp = IJ.openImage(stackPath)
cal = imp.getCalibration() 
print(cal)
imp.show()
IJ.run("Enhance Contrast", "saturated=0.35");
IJ.run("Apply LUT", "stack");
imp.changes = 0 # will not ask to save on close


impMask = IJ.openImage(maskPath)
impMask.setCalibration(cal.copy())
impMask.show()
IJ.run("Enhance Contrast", "saturated=0.35");
IJ.run("Apply LUT", "stack");
impMask.changes = 0 # will not ask to save on close

impSkel = IJ.openImage(skelPath)
impSkel.setCalibration(cal.copy())
impSkel.show()
IJ.run("Enhance Contrast", "saturated=0.35");
IJ.run("Apply LUT", "stack");
impSkel.changes = 0 # will not ask to save on close


# Create a universe and show it
univ = Image3DUniverse(512, 512)  


c1 = univ.addVoltex(imp)
myColor = Color3f((0.0,1.0,0.0))
c1.setColor(myColor)
c1.setThreshold(5)
#print('c1:', type(c1), c1)



c2 = univ.addVoltex(impMask)
myColor = Color3f((0.0,0.0,0.0))
c2.setColor(myColor)



c3 = univ.addVoltex(impSkel)
myColor = Color3f((1.0,0.0,0.0))
c2.setColor(myColor)




univ.show();