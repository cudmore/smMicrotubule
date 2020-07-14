"""
20200530
test aics segmentation on vasculature

on 'pip install napari'

	ERROR: aicsimageio 0.6.4 has requirement tifffile==0.15.1, but you'll have tifffile 2020.5.25 which is incompatible.

ERROR: aicsimageio 0.6.4 has requirement tifffile==0.15.1, but you'll have tifffile 2020.5.25 which is incompatible.
ERROR: aicsimageprocessing 0.7.3 has requirement aicsimageio>=3.1.2, but you'll have aicsimageio 0.6.4 which is incompatible.

"""

import sys, time

import numpy as np
import scipy

import tifffile

import napari

#import bimpy

from aicssegmentation.core.vessel import filament_3d_wrapper
from aicssegmentation.core.pre_processing_utils import intensity_normalization, edge_preserving_smoothing_3d, image_smoothing_gaussian_3d

from my_suggest_normalization_param import my_suggest_normalization_param # clone of aics segmentation

def _printStackParams(name, myStack):
	print('  ', name, myStack.shape, myStack.dtype, 'dtype.char:', myStack.dtype.char,
		'min:', np.min(myStack),
		'max:', np.max(myStack),
		'mean:', np.mean(myStack),
		'std:', np.std(myStack),
		)

def slidingZ(imageStack, upDownSlices=1, verbose=False):
	print('  .slidingZ() upDownSlices:', upDownSlices)

	if verbose: _printStackParams('input', imageStack)

	slidingz = imageStack.copy()
	m = imageStack.shape[0]
	for i in range(m):
		firstSlice = i - upDownSlices
		lastSlice = i + upDownSlices
		if firstSlice < 0:
			firstSlice = 0
		if lastSlice > m-1:
			lastSlice = m-1
		
		slidingz[i,:,:] = np.max(imageStack[firstSlice:lastSlice+1,:,:], axis=0)
		
	if verbose: _printStackParams('output', slidingz)

	return slidingz

def medianFilter(imageStack, kernelSize=(2,2,2), verbose=False):
	"""
	"""
	print('  .myMedianFilter() kernelSize:', kernelSize, '... please wait')
	
	if verbose: _printStackParams('input', imageStack)
		
	startTime = time.time()
	
	result = scipy.ndimage.median_filter(imageStack, size=kernelSize)
	
	if verbose: _printStackParams('output', result)

	stopTime = time.time()
	print('    myMedianFilter took', round(stopTime-startTime,2), 'seconds')
	
	return result

def myRun(path):

	"""
		scale_x is set based on the estimated thickness of your target filaments.
		For example, if visually the thickness of the filaments is usually 3~4 pixels,
		then you may want to set scale_x as 1 or something near 1 (like 1.25).
		Multiple scales can be used, if you have filaments of very different thickness.
	
		cutoff_x is a threshold applied on the actual filter reponse to get the binary result.
		Smaller cutoff_x may yielf more filaments, especially detecting more dim ones and thicker segmentation,
		while larger cutoff_x could be less permisive and yield less filaments and slimmer segmentation.
	"""
	f3_param=[[1, 0.01]]
	f3_param=[[5, 0.001], [3, 0.001]]
	
	stackData = tifffile.imread(path)
	numSlices = stackData.shape[0]
	_printStackParams('stackData', stackData)
	
	
	stackData = slidingZ(stackData, upDownSlices=1)
	
	stackData = medianFilter(stackData)

	# give us a guess for our intensity_scaling_param parameters
	low_ratio, high_ratio = my_suggest_normalization_param(stackData)

	# try per slice
	normData = stackData.astype(np.float64)
	normData[:] = 0
	for i in range(numSlices):
		oneSlice = stackData[i,:,:]
		low_ratio, high_ratio = my_suggest_normalization_param(oneSlice)
		
		
		print(i, low_ratio, high_ratio)
		
		#low_ratio = 0.2
		low_ratio -= 0.2
		high_ratio -= 1
		
		intensity_scaling_param = [low_ratio, high_ratio]
		sliceNormData = intensity_normalization(oneSlice, scaling_param=intensity_scaling_param)
		normData[i,:,:] = sliceNormData
		
	#sys.exit()
	
	'''
	#intensity_scaling_param = [0.0, 22.5]
	intensity_scaling_param = [low_ratio, high_ratio]
	print('    === intensity_normalization() intensity_scaling_param:', intensity_scaling_param)
	
	# intensity normalization
	print('    === calling intensity_normalization()')
	normData = intensity_normalization(stackData, scaling_param=intensity_scaling_param)
	_printStackParams('normData', normData)
	'''
	
	# smoothing with edge preserving smoothing 
	print('    === calling edge_preserving_smoothing_3d()')
	smoothData = edge_preserving_smoothing_3d(normData)
	_printStackParams('smoothData', smoothData)

	print('    === calling filament_3d_wrapper() f3_param:', f3_param)
	filamentData = filament_3d_wrapper(smoothData, f3_param)
	#filamentData = filamentData > 0
	_printStackParams('filamentData', filamentData)

	#filamentData2 = slidingZ(filamentData, upDownSlices=1)
	
	#
	# napari
	print('opening in napari')
	scale = (1, 0.6, 0.6)
	with napari.gui_qt():
		viewer = napari.Viewer(title='xxx')
		
		minContrast = 0
		maxContrast = 255
		myImageLayer = viewer.add_image(stackData, scale=scale, contrast_limits=(minContrast, maxContrast), colormap='green', visible=True, name='stackData')
		
		minContrast = 0
		maxContrast = 1
		myImageLayer = viewer.add_image(normData, scale=scale, contrast_limits=(minContrast, maxContrast), opacity=0.6, colormap='gray', visible=True, name='normData')
		
		minContrast = 0
		maxContrast = 1
		myImageLayer = viewer.add_image(filamentData, scale=scale, contrast_limits=(minContrast, maxContrast), opacity=0.6, colormap='blue', visible=True, name='filamentData')
		
		'''
		minContrast = 0
		maxContrast = 1
		myImageLayer = viewer.add_image(filamentData2, scale=scale, contrast_limits=(minContrast, maxContrast), opacity=0.6, colormap='red', visible=True, name='filamentData2')
		'''

if __name__ == '__main__':

	path = '/Users/cudmore/box/data/nathan/20200518/analysis2/20200518__A01_G001_0003_ch2_raw.tif'
	myRun(path)