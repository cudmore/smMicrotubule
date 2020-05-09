# Author: Robert Cudmore
# Date: 20200424

import os

from ij import IJ
from ij.plugin import ChannelSplitter

from loci.plugins import BF

original_OIR = '/Users/cudmore/box/share/deconvolded-oir-errors/8_5ADVMLEG1L1.oir'
deconvolved_OIR = '/Users/cudmore/box/share/deconvolded-oir-errors/8_5ADVMLEG1L1.oir'

# choose either oringal or deconvolved
#filePath = original_OIR
filePath = deconvolved_OIR

# set up a base name to save channels 
# We will append one of ('_ch1.tif', '_ch2.tif', '_ch3.tif') to saveFileNoExtension
folderPath, fileName = os.path.split(filePath)
fileNameNoExtension = os.path.splitext(fileName)[0]
saveFileNoExtension = os.path.join(folderPath, fileNameNoExtension)

# open the .oir
imps = BF.openImagePlus(filePath) # list with one imp
if len(imps)==1:
	imp = imps[0]

# split the channel
channelArray = ChannelSplitter.split(imp) # returns an array of imp

for channelIdx, channelImp in enumerate(channelArray):
	saveFilePath = saveFileNoExtension + '_ch' + str(channelIdx+1) + '.tif'
	print('saving saveFilePath:', saveFilePath)
	IJ.save(channelImp, saveFilePath)
