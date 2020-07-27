"""

20200716, rewrite of microtubule mask to get (cytosol mask, membrane mask)

Problem: now is fullmask sometimes gets portion of nucleus
and then on subtracting the nucleus, they are left over

Solution: add another step (1) label mask 1,2,3,... and remove small labels?

"""

import os, sys, math
from collections import OrderedDict

import numpy as np
import pandas as pd

import multiprocessing as mp

import napari

import bimpy

from gAnalysisPath import gSami_Params
from samiMaskNucleus import samiMaskNucleus

import samiUtil

def getDefaultParamDict():
	paramDict = OrderedDict()
	paramDict['gaussianKernel'] = (3,4,4)
	paramDict['minThreshold'] = 1 # filtered values > than this go into mask
	paramDict['erodeIterations'] = 1
	paramDict['removeNumSlicesFromMask'] = 0 # 0 will remove one oon top and one of bottom
	paramDict['cytosolErode'] = 4
	
	paramDict['nucleus'] = OrderedDict()
	paramDict['nucleus']['erodeIterations'] = 2
	paramDict['nucleus']['removeNumSlicesFromMask'] = 4
	
	return paramDict
	
def samiMaskMicrotubules(path, paramDict, doSave=False, doNapari=False):
	"""
	"""
	print('=== samiMaskMicrotubules() path:', path)

	# todo: add path to this
	# for now, assuming all keys are stack/mask data ???
	stackResults = OrderedDict()
	
	outPath, filename = os.path.split(path) # so we save to same analysis folder, like '/Users/cudmore/Desktop/samiVolume3'
	fileNameNoExtension, tmpExt = filename.split('.')
	baseSavePath = os.path.join(outPath, fileNameNoExtension)
	print('  baseSavePath:', baseSavePath)

	#
	# load and analyze the nucleus channel (we will eventually mask it)
	nucleusPath = path.replace('_ch2.tif', '_ch1.tif')
	nucleusResults = samiMaskNucleus(nucleusPath, paramDict)

	stackResults['nucleusMask'] = nucleusResults['erodedData']
	stackResults['nucleusData'] = nucleusResults['stackData']
	
	#
	# load _ch2 from path
	currentStack, tiffHeader = bimpy.util.bTiffFile.imread(path)
	stackResults['stackData'] = currentStack
	
	print('  xVoxel:', tiffHeader['xVoxel'], 'yVoxel:', tiffHeader['yVoxel'], 'zVoxel:', tiffHeader['zVoxel'])
	
	#
	# smooth
	kernelSize = paramDict['gaussianKernel'] #(3,4, 4)
	currentStack = bimpy.util.morphology.gaussianFilter(currentStack, kernelSize=kernelSize)
	stackResults['smoothedData'] = currentStack
	
	minSmooth = np.min(currentStack)
	maxSmooth = np.max(currentStack)
	rangeSmooth = maxSmooth - minSmooth + 1
	print('  minSmooth:', minSmooth, 'maxSmooth:', maxSmooth, 'rangeSmooth:', rangeSmooth, currentStack.dtype)
	
	# todo: report min/max of filtered
	
	#
	# threshold
	# remember, my threshold otsu does slice by slice, does not seem to be a 3d otsu?
	'''
	# otsu misses a bunch of dim microtubules
	currentStack = bimpy.util.morphology.threshold_otsu(currentStack)
	stackResults['thresholdOtsuData'] = currentStack
	'''
	
	# simple min, hard to empirically determine value, that is what otsu is for?
	minThreshold = paramDict['minThreshold'] #3 # mask is where image is strictly greater than threshold_min
	print('  minThreshold:', minThreshold)
	currentStack = bimpy.util.morphology.threshold_min(currentStack, minThreshold)
	#stackResults['thresholdData'] = currentStack
	
	#
	# fill holes
	# important to do this first, as erosion (next) will actually expand inner holes
	currentStack = bimpy.util.morphology.binary_fill_holes(currentStack)
	#stackResults['filledStack'] = currentStack
	
	#
	# erode
	# we erode here because our initial gaussian filter smoothed/moved the outer edges of the cell
	# versus using a median filter which would preserve the edges but leave lots of holes (in the mask)
	# using erodedIterations of 4 because our gaussian was (5,5) in x/y
	# change erodedIterations if we change gaussian filter
	erodeIterations = paramDict['erodeIterations'] #1
	#erodedData = bimpy.util.morphology.fillHoles(filledStack, doDilation=False)
	currentStack = bimpy.util.morphology.binary_erosion(currentStack, erodeIterations)
	
	#
	# remove top/bottom slices from mask to get rid of z-spread artifact
	removeNumSlicesFromMask = paramDict['removeNumSlicesFromMask']
	currentStack = 	samiUtil.removeZSpreadFromMask(currentStack,
							removeNumSlices=removeNumSlicesFromMask, myName='microtubules')

	stackResults['fullMask'] = currentStack

	#
	# remove nucleus mask from fullMask
	fullMaskNoNucleus = stackResults['fullMask'] + stackResults['nucleusMask']
	print('  fullMaskNoNucleus:', np.min(fullMaskNoNucleus), np.max(fullMaskNoNucleus))
	fullMaskNoNucleus[fullMaskNoNucleus==2] = 0
	stackResults['fullMaskNoNucleus'] = fullMaskNoNucleus

	# label fullMaskNoNucleus and print size of each label
	# we are getting isolated portions of the nucleus erroneously in the mask
	labeledMask = bimpy.util.morphology.labelMask(fullMaskNoNucleus)
	labelMin = np.min(labeledMask)
	labelMax = np.max(labeledMask)
	print('  labeledMask min:', labelMin, 'max:', labelMax)
	for i in range(labelMax+1):
		numInLabel = np.count_nonzero(labeledMask==i)
		print('   label:', i, 'has', numInLabel)
	
	# todo: report number of labels
	
	#
	# make ring and cytosol masks
	cytosolErode = paramDict['cytosolErode']
	cytosolMask = bimpy.util.morphology.binary_erosion(currentStack, cytosolErode)
	#stackResults['cytosolMask'] = cytosolMask
	
	ringMask = cytosolMask + currentStack
	ringMask[ringMask==2] = 0
	stackResults['ringMask'] = ringMask
	
	#
	# remove nucleus mask from cytosol mask
	cytosolMask = cytosolMask + stackResults['nucleusMask']
	cytosolMask[cytosolMask==2] = 0
	#stackResults['cytosolMask'] = cytosolMask
	stackResults['erodedMask'] = cytosolMask
	
	if doSave:
		# _fullMask 
		fullMaskPath = baseSavePath + '_fullMask.tif'
		print('  saving fullMaskPath:', fullMaskPath, stackResults['fullMask'].dtype)
		bimpy.util.bTiffFile.imsave(fullMaskPath, stackResults['fullMask'],
									tifHeader=tiffHeader, overwriteExisting=True)
		
		# _fullMaskNoNucleus 
		fullMaskNoNucleusPath = baseSavePath + '_fullMaskNoNucleus.tif'
		print('  saving fullMaskNoNucleusPath:', fullMaskNoNucleusPath, stackResults['fullMaskNoNucleus'].dtype)
		bimpy.util.bTiffFile.imsave(fullMaskNoNucleusPath, stackResults['fullMaskNoNucleus'],
									tifHeader=tiffHeader, overwriteExisting=True)
		
		# _erodedMask (e.g. cytosol mask)
		erodedMaskPath = baseSavePath + '_erodedMask.tif'
		print('  saving erodedMaskPath:', erodedMaskPath, stackResults['erodedMask'].dtype)
		bimpy.util.bTiffFile.imsave(erodedMaskPath, stackResults['erodedMask'],
									tifHeader=tiffHeader, overwriteExisting=True)
		
		# _ringMask
		ringMaskPath = baseSavePath + '_ringMask.tif'
		print('  saving ringMaskPath:', ringMaskPath, stackResults['ringMask'].dtype)
		bimpy.util.bTiffFile.imsave(ringMaskPath, stackResults['ringMask'],
									tifHeader=tiffHeader, overwriteExisting=True)
		
	if doNapari:
		myNapari(paramDict['index'], path, stackResults)
		
	return True
	
def myNapari(idx, path, results):
	print('  myNapari()')
	myScale = (3, 1, 1)
	with napari.gui_qt():
		title = str(idx) + ' ' + path
		viewer = napari.Viewer(ndisplay=2, title=title)
		for k,v in results.items():
			include = True
			minContrast = 0
			maxContrast = 1
			theMax = np.max(v)
			if theMax == 1:
				maxContrast = 1 # binary mask
			elif theMax>250:
				maxContrast = 60 # 8-bit image
			else:
				maxContrast = theMax + 1 # labeled stack
			visible = True
			colormap = 'gray'
			if k == 'stackData':
				colormap = 'green'
			elif k == 'smoothedData':
				colormap = 'blue'
				visible = False
				minContrast = 0
				maxContrast = 5
			elif k in ['thresholdData', 'thresholdOtsuData']:
				colormap = 'inferno'
				visible = False
			elif k == 'fullMask':
				include = False
				colormap = 'gray'
				visible = False
			
			elif k == 'erodedMask': # e.g. cytosolMask
				include = False
				colormap = 'yellow'
			elif k == 'ringMask':
				include = False
				visible = False
				colormap = 'cyan'

			elif k == 'nucleusData':
				colormap = 'red'
			elif k == 'nucleusMask':
				colormap = 'gray_trans'
				visible = False

			if maxContrast == 1:
				opacity = 0.5
			else:
				opacity = 1.0

			if include:
				viewer.add_image(v, scale=myScale, contrast_limits=(minContrast,maxContrast),
					opacity=opacity, colormap=colormap, visible=visible, blending='additive', name=k)

if __name__ == '__main__':

	if 0:
		path = '/Users/cudmore/Sites/smMicrotubule/data/200108/WT_Female/Cell_1/1_5ADVMLEG1L1_ch2.tif'
		#path = '/Users/cudmore/Sites/smMicrotubule/data/200108/WT_Female/Cell_7/7_5ADVMLEG1L1_ch2.tif'
		
		# don't escape the spaces
		#path = '/Users/cudmore/Sites/smMicrotubule/data/200414/Bin1smKO Male/Cell 1/Male KO cell 1_3ADVMLEG1L1_ch2.tif'
		#path = '/Users/cudmore/Sites/smMicrotubule/data/191230/BIN1_smKO_Male/Cell_12/12_5ADVMLEG1L1_ch2.tif'
		
		samiMaskMicrotubules(path, doNapari=True)
			
	if 1:
		# for output analysis
		dataPath = gSami_Params['gAnalysisPath']
		#dataPath = '/Users/cudmore/Sites/smMicrotubule/data'
		
		batchFileList = gSami_Params['gBatchFileList']
				
		paramDict = getDefaultParamDict()
		
		# to specify just one file from row of ../analysis/maskParams.csv
		thisFileNum = None

		####
		####
		# file 30 needs to use less erosion on nucleus ?????
		thisFileNum = 135 # if specified will NOT run in parallel
		####
		####
		
		paramDict['index'] = thisFileNum
		
		#
		doSave = False
		doNapari = True
		
		maskParamsPath = '../analysis/maskParams.csv'
		dfParams = pd.read_csv(maskParamsPath, header=0)
		#print(dfParams)
		numRows = dfParams.shape[0]

		#
		# parallel
		cpuCount = mp.cpu_count()
		cpuCount -= 2
		pool = mp.Pool(processes=cpuCount)
		pResults = []
		
		for i in range(numRows):
			# do just one file
			if thisFileNum is not None and i != thisFileNum:
				continue

			dfRow = dfParams.iloc[i]
			include = dfRow['include']
			minThreshold = dfRow['minThreshold']
			path = dfRow['path']
			
			# this MUST be specified
			#print('minThreshold:', str(minThreshold), minThreshold.dtype, str(minThreshold)=='nan')
			if str(minThreshold) == 'nan':
				# don't use
				print('!!! ignoring file', i, 'minThreshold is not specified, path:', path)
				#pass
				continue
			else:
				paramDict['minThreshold'] = minThreshold
			
			# include or not
			doInclude = True
			if isinstance(include, bool):
				doInclude = include
			if not doInclude:
				print('!!! ignoring file', i, 'doInclude:', doInclude, 'path:', path)
				continue
				
			#
			# actualPath can be in Desktop/samiVolume3/'
			actualPath = os.path.join(dataPath,path) # dataPath is like '/Users/cudmore/Sites/smMicrotubule/data'
			# make sure it exists
			if not os.path.isfile(actualPath):
				print('ERROR: did not find file', actualPath)
				continue
				
			print(maskParamsPath, 'row:', i, 'minThreshold:', minThreshold, actualPath)
		
			#samiMaskMicrotubules(actualPath, paramDict, doSave=True, doNapari=doNapari)
			
			if thisFileNum is None:
				# parallel
				args = [actualPath, paramDict, doSave, doNapari]
				oneAnswer = pool.apply_async(samiMaskMicrotubules, args=args)
				pResults.append(oneAnswer)
			else:
				# do one file
				samiMaskMicrotubules(actualPath, paramDict, doSave, doNapari)
				
		#sys.exit()
		
		# run in parallel
		for pResult in pResults:
			oneResult = pResult.get()
			#if tiffFileInfo is not None:
			#	fileDictList.append(tiffFileInfo)

