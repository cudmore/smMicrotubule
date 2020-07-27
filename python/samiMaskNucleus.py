"""

	mask the nucleus
	
	todo: check scipy morphology for structure parameter ...
			this is causes problems above/below first/last slice???
"""

from collections import OrderedDict
import numpy as np

import napari

import bimpy

import samiUtil

def samiMaskNucleus(path, paramDict, doSave=False, doNapari=False):
	"""
	"""
	print('  samiMaskNucleus() path:', path)
	
	stackResults = OrderedDict()

	#
	# load
	stackData, tiffHeader = bimpy.util.bTiffFile.imread(path)
	stackResults['stackData'] = stackData
	
	#
	# smooth
	kernelSize = (3,5,5)
	smoothedData = bimpy.util.morphology.gaussianFilter(stackData, kernelSize=kernelSize)
	stackResults['smoothedData'] = smoothedData
	
	#
	# threshold
	# remember, my threshold otsu does slice by slice, does not seem to be a 3d otsu?
	thresholdData = bimpy.util.morphology.threshold_otsu(smoothedData, verbose=False)
	stackResults['thresholdData'] = thresholdData
	
	#
	# fill holes
	# important to do this first, as erosion (next) will actually expand inner holes
	filledStack = bimpy.util.morphology.binary_fill_holes(thresholdData)
	stackResults['filledStack'] = filledStack
	
	#
	# erode
	# we erode here because our initial gaussian filter smoothed/moved the outer edges of the cell
	# versus using a median filter which would preserve the edges but leave lots of holes (in the mask)
	# using erodedIterations of 4 because our gaussian was (5,5) in x/y
	# change erodedIterations if we change gaussian filter
	#erodedIterations = 4
	erodedData = bimpy.util.morphology.fillHoles(filledStack, doDilation=False)

	erodeIterations = paramDict['nucleus']['erodeIterations'] # 1
	if erodeIterations > 0:
		# usually do not do this
		erodedData = bimpy.util.morphology.binary_erosion(erodedData, erodeIterations)
	
	#
	# remove top/bottom 3x slices (fixes z-spread?)
	removeNumSlicesFromMask = paramDict['nucleus']['removeNumSlicesFromMask']
	erodedData = samiUtil.removeZSpreadFromMask(erodedData, removeNumSlices=removeNumSlicesFromMask, myName='nucleus')
		
	stackResults['erodedData'] = erodedData

	if doNapari:
		myNapari(stackResults)
	
	return stackResults
	
def myNapari(results):
	print('myNapari()')
	myScale = (3, 1, 1)
	with napari.gui_qt():
		viewer = napari.Viewer(ndisplay=2)
		for k,v in results.items():
			minContrast = 0
			maxContrast = 1
			theMax = np.max(v)
			if theMax == 1:
				maxContrast = 1 # binary mask
			elif theMax>250:
				maxContrast = 60 # 8-bit image
			else:
				maxContrast = theMax + 1 # labeled stack
			colormap = 'gray'
			if k == 'stackData':
				colormap = 'red'
			elif k == 'smoothedData':
				colormap = 'blue'
			elif k == 'thresholdData':
				colormap = 'inferno'
			# colormap
			viewer.add_image(v, scale=myScale, contrast_limits=(minContrast,maxContrast), opacity=0.5, colormap=colormap, visible=True, blending='additive', name=k)

if __name__ == '__main__':

	#path = '/Users/cudmore/Sites/smMicrotubule/data/200108/WT_Female/Cell_1/1_5ADVMLEG1L1_ch1.tif'
	#path = '/Users/cudmore/Sites/smMicrotubule/data/200108/WT_Female/Cell_7/7_5ADVMLEG1L1_ch1.tif'
	
	# check what is going on with this cell
	# 3 channel, channel 1 is corrupt?
	#path = '/Users/cudmore/Sites/smMicrotubule/data/200108/WT_Female/Cell_8/8_5ADVMLEG1L1_ch1.tif'

	# don't escape the spaces
	path = '/Users/cudmore/Sites/smMicrotubule/data/200414/Bin1smKO Male/Cell 1/Male KO cell 1_3ADVMLEG1L1_ch1.tif'
		
	samiMaskNucleus(path, doNapari=True)