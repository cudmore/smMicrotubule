"""
#
# try and fill holes from original mask

Working on this for 4+ hours, I have no idea what is going on?

# maybe this
# scipy.ndimage.measurements.watershed_ift
"""

import os
import tifffile

import numpy as np

#from scipy.ndimage.morphology import binary_fill_holes
import scipy.ndimage.morphology

#from skimage.morphology import local_maxima

import skimage.morphology
#from skimage import measure

from scipy import ndimage

def myImageStat(name, data, find=None):
	print(name, ',', np.min(data), ',', np.max(data), ',', np.mean(data), ',', 'find=', find, ',' , np.count_nonzero(data==find), ',', data.shape, ',', data.dtype)
	
def myFillHoles(filePath):
	
	folderPath, fileName = os.path.split(filePath)
	fileNameNoExtension, extension = os.path.splitext(fileName)
	
	maskFile = fileNameNoExtension + '_dvMask.tif'
	maskPath = os.path.join(folderPath, maskFile)
	
	print(',', 'min', ',', 'max', ',', 'mean', ',', 'find', ',', 'count', ',', 'shape', ',' 'type')
	
	# remember, _dvMask.tif is always 0/255
	maskData = tifffile.imread(maskPath)
	n1 = np.count_nonzero(maskData==255)
	myImageStat('maskData', maskData, find=255)
	
	# only works on 2d
	#contours = measure.find_contours(maskData, 0.8)

	#
	# find_objects
	tmpOut = np.zeros(maskData.shape, dtype=np.uint8)
	slicesLoc = ndimage.find_objects(maskData) # returns list of slice where slice(start,stop,step)
	for oneSlice in slicesLoc:
		if oneSlice is None:
			continue
		slicesImg = maskData[oneSlice]
		tmpOut[oneSlice] = 1
	myImageStat('tmpOut', tmpOut, find=1)
	
	#tmpOut = tmpOut.astype(int8)
	#myImageStat('tmpOut2', tmpOut, find=1)
	
	findObjectPath = os.path.join(folderPath, fileNameNoExtension + '_dvFindObjects.tif')
	tifffile.imsave(findObjectPath, tmpOut)

	#
	# watershed
	distance = ndimage.distance_transform_edt(maskData)
	myImageStat('distance', distance, find=True)
	# required 2d
	#from skimage.feature import peak_local_max
	#local_maxi = peak_local_max(distance, indices=False, footprint=np.ones((3, 3)), labels=maskData)
	local_maxi = skimage.morphology.local_maxima(distance)
	myImageStat('local_maxi', local_maxi)
	
	footPrintTuple = (3,3,3)
	markers = ndimage.label(local_maxi, structure=np.ones(footPrintTuple))[0] # markers is a stack
	labels = skimage.morphology.watershed(-distance, markers, mask=maskData)
	labels = labels.astype(np.uint8)
	myImageStat('labels', labels)

	labelsPath = os.path.join(folderPath, fileNameNoExtension + '_dvLabels.tif')
	tifffile.imsave(labelsPath, labels)


	markersPath = os.path.join(folderPath, fileNameNoExtension + '_dvMarkers.tif')
	tifffile.imsave(markersPath, markers)


	#
	# fill holes
	'''
	maskFilled = scipy.ndimage.morphology.binary_fill_holes(maskData)
	nFilled = np.count_nonzero(maskFilled==255)
	myImageStat('maskFilled', maskFilled, find=True)
	'''
	
	#
	# dilate
	'''
	maskDilated = scipy.ndimage.morphology.binary_dilation(maskData)
	nDilated = np.count_nonzero(maskDilated==255)
	myImageStat('maskDilated', maskDilated, find=True)
	
	dilatedPath = os.path.join(folderPath, fileNameNoExtension + '_dvDilated.tif')
	tifffile.imsave(dilatedPath, maskDilated)
	'''
	
if __name__ == '__main__':

	filePath = '/Users/cudmore/box/data/sami/191230/BIN1_smKO_Male/Cell_12/12_5ADVMLEG1L1_ch2.tif'	
	myFillHoles(filePath)