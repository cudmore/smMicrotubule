"""
Robert Cudmore
20200330

Goal here is to run aics segmentation with different parameters and ask if we get different results.

Parameters that go into making a 3D mask:
	intensity_scaling_param : used by intensity_normalization()
	edge_preserving_smoothing_3d() has no parameters
	f3_param : used by filament_3d_wrapper()
	
In particular:
	Different number of segments?
	Statistically different segment lengths?
	Statistically different segment tortuosity
	
What this is missing is a comparison of 3D mask versus raw image data, for example:
(Not sure how to fix this ... keep working ...)
	Are we tracing things not in the image?
	Are we not tracing things that are really in the image?
	

First thing to try is
	f3_param = [[1, 0.01]]
	f3_param = [[1.25, 0.01]]

To run this code:

	on my home desktop, virtual env is in:
		/Users/cudmore/Sites/my-aics/aics-segmentation/aicssegmentation

	to start virtual env:
		cd /Users/cudmore/Sites/my-aics/aics-segmentation
		conda activate segmentation

	change back into folder where this file lives:
		cd /Users/cudmore/box/data/sami/
	
	run it:
		#python samiAnalysis.py
		python samiAnalysis.py batch=analysis/wt-male.txt 
	
	remember, I was using original recipe from jupyter notebook in:
		notebooks/playground_filament3d.ipynb
"""

import os, sys, time
from collections import OrderedDict

import numpy as np
import pandas as pd
import tifffile

#from scipy.stats import norm # used by my_suggest_normalization_param

# statistical tests between analysis parameters
# also try mannwhitneyu
#from scipy.stats import wilcoxon # used test sig difference between (length, euclidean distance)
from scipy.stats import ttest_ind # used test sig difference between (length, euclidean distance)


# package for 3d visualization
from itkwidgets import view							  
from aicssegmentation.core.visual import seg_fluo_side_by_side,  single_fluorescent_view, segmentation_quick_view
import matplotlib.pyplot as plt
plt.rcParams["figure.figsize"] = [16, 12]

# package for io 
from aicsimageio import AICSImage, omeTifWriter							

# function for core algorithm
from aicssegmentation.core.vessel import filament_3d_wrapper
from aicssegmentation.core.pre_processing_utils import intensity_normalization, edge_preserving_smoothing_3d, image_smoothing_gaussian_3d
from skimage.morphology import remove_small_objects	

from skimage import feature
from skimage import morphology # to convert 3d mask to 1-pixel skeleton
import skan

from IPython.display import display # to show pandas dataframe as a table

# my code
from my_suggest_normalization_param import my_suggest_normalization_param
from readVoxelSize import readVoxelSize

def myAnalyzeSkeleton(out=None, maskPath=None, imagePath=None):
	"""
	out: numpy array with 1-pixel skeleton
	maskPath : full path to _dvMask.tif file (can include appended _0.tif
	"""
	
	# load x/y/z voxel size (assumes .tif was saved with Fiji
	# we use this to scale length
	xVoxel, yVoxel, zVoxel = readVoxelSize(imagePath)

	# load the mask
	if out is not None:
		maskData = out
	else:
		#maskPath = os.path.splitext(path)[0] + '_dvMask_' + str(saveNumber) + '.tif'
		maskData = tifffile.imread(maskPath)
	
	# was used by shape_index
	#imageData = tifffile.imread(imagePath)
	
	print('=== myAnalyzeSkeleton() maskData.shape:', maskData.shape)
	
	# make a 1-pixel skeleton from volume mask (similar to Fiji Skeletonize)
	mySkeleton = morphology.skeletonize_3d(maskData)

	'''
	# shape_index() does not work for 3D images !!!
	scale = 1
	threshold_radius = 1 # from AICS
	smooth_radius =  0.01 # from AICS
	pixel_threshold_radius = int(np.ceil(threshold_radius / scale))
	pixel_smoothing_radius = smooth_radius * pixel_threshold_radius
	quality = feature.shape_index(imageData, sigma=pixel_smoothing_radius, mode='reflect')
	#skeleton = morphology.skeletonize(thresholded) * quality
	mySkeleton = morphology.skeletonize_3d(maskData) * quality
	'''
	
	# analyze the skeleton (similar to Fiji Analyze Skeleton)
	mySkanSkel = skan.Skeleton(mySkeleton)

	# look at the results
	branch_data = skan.summarize(mySkanSkel) # branch_data is a pandas dataframe
	nBranches = branch_data.shape[0]
	
	'''
	print('    number of branches:', branch_data.shape[0])
	display(branch_data.head())
	'''
	
	#
	# convert everything to nump arrays
	branchDistance = branch_data['branch-distance'].to_numpy()
	euclideanDistance = branch_data['euclidean-distance'].to_numpy()
	branchType = branch_data['branch-type'].to_numpy()
	#tortuosity = branchDistance / euclideanDistance # this gives divide by 0 warning
	tmpOut = np.full_like(branchDistance, fill_value=np.nan)
	tortuosity = np.divide(branchDistance, euclideanDistance, out=tmpOut, where=euclideanDistance!=0)



	"""
	
	Sunday 20200405
	HERE I AM RUNNING CODE TWICE and APPENDING TO summary2.xlsx AFTER RUNNING samiMetaAnalysis.py
	
	https://jni.github.io/skan/_modules/skan/pipe.html#process_images
	in the Skan Pipe code, they multiply the binary skeleton as follows
	maybe I can implement this with scale=1 and threshold_radius taken from AICS Segmentaiton?
	
		scale = 1
		threshold_radius = 1 # from AICS
		smooth_radius =  0.01 # from AICS
		pixel_threshold_radius = int(np.ceil(threshold_radius / scale))
		pixel_smoothing_radius = smooth_radius * pixel_threshold_radius
		quality = skimage.feature.shape_index(image, sigma=pixel_smoothing_radius,
							  mode='reflect')
		skeleton = morphology.skeletonize(thresholded) * quality
	"""
	'''
	# 20200407, no longer needed as i am now saving 'branch-type', do this in meta analysis
	print('\n\n\t\tREMEMBER, I AM ONLY INCLUDING junction-to-junction !!!!!!!!!!!!!! \n\n')
	#
	# do again just for junction-to-junction
	# 'mean-pixel-value' here is 'mean shape index' in full tutorial/recipe
	# if I use ridges, I end up with almost no branche?
	#ridges = ((branch_data['mean-pixel-value'] < 0.625) & (branch_data['mean-pixel-value'] > 0.125))
	j2j = branch_data['branch-type'] == 2 # returns True/False pandas.core.series.Series
	#datar = branch_data.loc[ridges & j2j].copy()
	datar = branch_data.loc[j2j].copy()

	branchDistance = datar['branch-distance'].to_numpy()
	euclideanDistance = datar['euclidean-distance'].to_numpy()
	#tortuosity = branchDistance / euclideanDistance # this gives divide by 0 warning
	tmpOut = np.full_like(branchDistance, fill_value=np.nan)
	tortuosity = np.divide(branchDistance, euclideanDistance, out=tmpOut, where=euclideanDistance != 0)
	'''
	
	#
	# organize a return dicitonary
	retDict = OrderedDict()

	retDict['data'] = OrderedDict()
	#retDict['data']['nBranches'] = nBranches
	retDict['data']['branchLength'] = branchDistance
	retDict['data']['euclideanDistance'] = euclideanDistance
	retDict['data']['branchType'] = branchType
	#retDict['data']['tortuosity'] = tortuosity

	# todo: search for 0 values in (branchDistance, euclideanDistance)
	
	# stats
	'''
	print('***** THIS IS NOT SCALED ***')
	print('    branchDistance mean:', np.mean(branchDistance), 'SD:', np.std(branchDistance), 'n:', branchDistance.size)
	#
	decimalPlaces = 2
	retDict['stats'] = OrderedDict()
	retDict['stats']['branchLength_mean'] = round(np.mean(branchDistance),decimalPlaces)
	retDict['stats']['branchLength_std'] = round(np.std(branchDistance),decimalPlaces)
	retDict['stats']['branchLength_n'] = branchDistance.shape[0]
	tmpCount = branchDistance[branchDistance<=2]
	retDict['stats']['branchLength_n_2'] = tmpCount.shape[0]
	#
	retDict['stats']['euclideanDistance_mean'] = round(np.mean(euclideanDistance),decimalPlaces)
	retDict['stats']['euclideanDistance_std'] = round(np.std(euclideanDistance),decimalPlaces)
	retDict['stats']['euclideanDistance_n'] = euclideanDistance.shape[0]
	#
	retDict['stats']['tortuosity_mean'] = round(np.nanmean(tortuosity),decimalPlaces)
	retDict['stats']['tortuosity_std'] = round(np.nanstd(tortuosity),decimalPlaces)
	retDict['stats']['tortuosity_n'] = tortuosity.shape[0]
	'''
	
	return retDict, mySkeleton # returning mySkeleton so we can save it
	
def run(path, f3_param=[[1, 0.01]], minArea=20, saveNumber=0):
	"""
	use aicssegmentation to pre-process raw data and then make/save a 3D mask
	"""
	print('=== path:', path)
	
	# load x/y/z voxel size (assumes .tif was saved with Fiji
	xVoxel, yVoxel, zVoxel = readVoxelSize(path)
	print('    xVoxel:', xVoxel, 'yVoxel:', yVoxel, 'zVoxel:', zVoxel)
	
	# load the data
	reader = AICSImage(path) 
	IMG = reader.data.astype(np.float32)
	print('    IMG.shape:', IMG.shape)

	structure_channel = 0
	struct_img0 = IMG[0,structure_channel,:,:,:].copy()

	# give us a guess for our intensity_scaling_param parameters
	#from aicssegmentation.core.pre_processing_utils import suggest_normalization_param
	#suggest_normalization_param(struct_img0)
	low_ratio, high_ratio = my_suggest_normalization_param(struct_img0)

	#intensity_scaling_param = [0.0, 22.5]
	intensity_scaling_param = [low_ratio, high_ratio]
	print('*** intensity_normalization() intensity_scaling_param:', intensity_scaling_param)
	
	# intensity normalization
	print('=== calling intensity_normalization()')
	struct_img = intensity_normalization(struct_img0, scaling_param=intensity_scaling_param)

	# smoothing with edge preserving smoothing 
	print('=== calling edge_preserving_smoothing_3d()')
	structure_img_smooth = edge_preserving_smoothing_3d(struct_img)

	#
	"""
	see: notebooks/playground_filament3d.ipynb

	scale_x is set based on the estimated thickness of your target filaments.
		For example, if visually the thickness of the filaments is usually 3~4 pixels,
		then you may want to set scale_x as 1 or something near 1 (like 1.25).
		Multiple scales can be used, if you have filaments of very different thickness.
	cutoff_x is a threshold applied on the actual filter reponse to get the binary result.
		Smaller cutoff_x may yielf more filaments, especially detecting more dim ones and thicker segmentation,
		while larger cutoff_x could be less permisive and yield less filaments and slimmer segmentation.
	"""
	#f3_param = [[1, 0.01]] # [scale_1, cutoff_1]
	print('=== calling filament_3d_wrapper() f3_param:', f3_param)
	bw = filament_3d_wrapper(structure_img_smooth, f3_param)
		
	#
	#minArea = 20 # from recipe
	print('=== calling remove_small_objects() minArea:', minArea)
	seg = remove_small_objects(bw>0, min_size=minArea, connectivity=1, in_place=False)

	#
	# save original file again (with saveNumber
	saveNumberStr = ''
	if saveNumber>1:
		saveNumberStr = '_' + str(saveNumber)
		
	#
	# save mask
	seg = seg >0
	out=seg.astype(np.uint8)
	out[out>0]=255
	
	# save _dvMask
	maskPath = os.path.splitext(path)[0] + '_dvMask' + saveNumberStr + '.tif'
	print('=== saving 3D mask [WILL FAIL IF FILE EXISTS] as maskPath:', maskPath)
	try:
		writer = omeTifWriter.OmeTifWriter(maskPath)
		writer.save(out)
	except(OSError) as e:
		print('    error: file already exists, di dnot resave, maskPath:', maskPath)
		
	#
	# analyze skeleton, take a 3d mask and analyze as a 1-pixel skeleton
	retDict0, mySkeleton = myAnalyzeSkeleton(out=out, imagePath=path)
	retDict = OrderedDict()
	retDict['tifPath'] = path
	retDict['maskPath'] = maskPath
	retDict['tifFile'] = os.path.basename(path)
	retDict['xVoxel'] = xVoxel
	retDict['yVoxel'] = yVoxel
	retDict['zVoxel'] = zVoxel
	#
	retDict['params'] = OrderedDict()
	retDict['params']['saveNumber'] = saveNumber
	retDict['params']['intensity_scaling_param'] = intensity_scaling_param # calculated in my_suggest_normalization_param
	retDict['params']['f3_param'] = f3_param[0] # cludge, not sure where to put this. f3_param is a list of list but screws up my .csv output !!!
	retDict['params']['minArea'] = minArea

	retDict.update( retDict0 )

	# save 1-pixel skeleton: mySkeleton
	# save _dvSkel
	skelPath = os.path.splitext(path)[0] + '_dvSkel' + saveNumberStr + '.tif'
	print('=== saving 3D skel [WILL FAIL IF FILE EXISTS] as maskPath:', skelPath)
	try:
		writer = omeTifWriter.OmeTifWriter(skelPath)
		writer.save(mySkeleton)
	except(OSError) as e:
		print('    error: file already exists, di dnot resave, skelPath:', skelPath)
			
	return retDict
	
def parseBatchFile(path):
	with open(path) as f:
		content = f.readlines()
	# remove whitespace characters like `\n` at the end of each line
	content = [x.strip() for x in content] 
	return content
	
if __name__ == '__main__':
	
	#path = '/Users/cudmore/box/data/sami/Cell_1/1_5ADVMLEG1L1_ch2.tif'
	
	path = None
	pathList = None
	savePath = ''
	doBatchFile = None
	sheetName = None
	for i, arg in enumerate(sys.argv):
		if i == 0:
			# my .py filename (not full path)
			continue
		#print('i:', i, 'arg:', arg)
		if arg.startswith('file='):
			doBatchFile = False
			tmp, filePath = arg.split('=')
			path = filePath
			savePath = os.path.splitext(path)[0]
			fileNameNoExtension = os.path.splitext(os.path.basename(path))[0]
			sheetName = fileNameNoExtension
		elif arg.startswith('batch='):
			doBatchFile = True
			# open specified file and make a list of path
			#os.path.splitext(path)[0]
			tmp, batchFilePath = arg.split('=')
			pathList = parseBatchFile(batchFilePath)
			savePath = os.path.splitext(batchFilePath)[0]
			batchFileNoExtension = os.path.splitext(os.path.basename(batchFilePath))[0]
			sheetName = batchFileNoExtension
			
		#elif i == 2:
		#	saveNumber = arg
	
	if pathList is None:
		pathList = [path]
	
	resultsList = []
	#pathList = [path]
	for path in pathList:	
	
		#fileName = os.path.basename(path)
		#fileNameNoExtension = os.path.splitext(os.path.basename(path))[0]
		
		print('\n***')
		print('*** path:', path)
		#print('    fileName:', fileName)
		#print('    fileNameNoExtension:', fileNameNoExtension)
	
		#1
		saveNumber = 1
		f3_param=[[1, 0.01]]
		f3_param=[[0.5, 0.005]]
		minArea=20
		retDict = run(path, f3_param, minArea=minArea, saveNumber=saveNumber)
		resultsList.append(retDict)
		
		'''
		#2
		saveNumber = 2
		f3_param=[[1, 0.01]]
		#f3_param=[1, 0.01]
		#minArea=40
		minArea=20
		retDict = run(path, f3_param, minArea=minArea, saveNumber=saveNumber)
		resultsList.append(retDict)
		'''
		
		'''
		#3
		saveNumber = 3
		f3_param=[[1.25, 0.01]]
		minArea=20
		retDict = run(path, f3_param, minArea=minArea, saveNumber=saveNumber)
		resultsList.append(retDict)
	
		#4
		saveNumber = 4
		f3_param=[[1.25, 0.01]]
		minArea=40
		retDict = run(path, f3_param, minArea=minArea, saveNumber=saveNumber)
		resultsList.append(retDict)
		'''
		
	#
	# print results, this works
	'''
	print('RESULTS ...')
	for result in resultsList:
		print(result['maskPath'])
		print('    ', result['params'])
		print('    ', result['stats'])
		print('')
	'''
	
	#
	# stats
	#
	# 20200407, THESE ARE NOT SCALED
	#
	'''
	print('STATS ...')
	decimalPlaces = 2
	for idx, result in enumerate(resultsList):
		if idx == len(resultsList)-1:
			continue # could be break
		a = resultsList[idx]['data']['branchLength']
		b = resultsList[idx+1]['data']['branchLength']
		#fStat: If alternative is “two-sided”, the sum of the ranks of the differences above or below zero, whichever is smaller. Otherwise the sum of the ranks of the differences above zero.
		#fStat, pValue = wilcoxon(a, b, alternative='two-sided')
		# exception: ValueError: zero_method 'wilcox' and 'pratt' do not work if the x - y is zero for all elements.
		fStat, pValue = ttest_ind(a, b)
		result['bl f'] = round(fStat,decimalPlaces)
		result['bl p'] = round(pValue,decimalPlaces)
		#print('    branchLength f:', fStat, 'p:', pValue)

		a = resultsList[idx]['data']['euclideanDistance']
		b = resultsList[idx+1]['data']['euclideanDistance']
		fStat, pValue = ttest_ind(a, b)
		result['ed f'] = round(fStat,decimalPlaces)
		result['ed p'] = round(pValue,decimalPlaces)
		#print('    euclideanDistance f:', fStat, 'p:', pValue)
	
		a = resultsList[idx]['data']['tortuosity']
		b = resultsList[idx+1]['data']['tortuosity']
		a = a[~np.isnan(a)] # remove nan values, these occur when original euclidean distance was 0
		b = b[~np.isnan(b)]
		fStat, pValue = ttest_ind(a, b)
		result['t f'] = round(fStat,decimalPlaces)
		result['t p'] = round(pValue,decimalPlaces)
		#print('    tortuosity f:', fStat, 'p:', pValue)
	
		#print('result:', result)
	'''
	
	#
	# make flat output for viewing in excel
	# remember, commas (',') in list screw up auto open in excel !!!!
	# todo: simplify this, WAY to complicated !!!!!!!!!!!
	skipTheseKeys = ['data'] # skip keys that are actually numpy arrays

	headerStr = ''
	valueStr = ''
	thisDelim = ','
	#dfList0 = []
	for idx, result in enumerate(resultsList):
		for k1,v1 in result.items():
			# skip ['data']
			if k1 in skipTheseKeys:
				continue

			'''
			print('k1:', k1)
			print('v1:', str(v1))
			'''
			
			#df0 = pd.DataFrame({k1:[v1]})
			#dfList0.append(df0)
			
			if isinstance(result[k1], dict):
				# append all k2,v2 in dict
				for k2,v2 in result[k1].items():
					'''
					print('    k2:', k2)
					print('    v2:', str(v2))
					'''
					if idx == 0:
						headerStr += k2 + thisDelim
					if isinstance(v2, list) or isinstance(v2, tuple):
						valueStr += '['
						for tmpItem in v2:
							valueStr += str(tmpItem) + ' ' 
						valueStr += ']' + thisDelim
					else:
						valueStr += str(v2) + thisDelim
			else:
				# just append k1, v1
				if idx == 0:
					headerStr += k1 + thisDelim
				if isinstance(v1, list) or isinstance(v1, tuple):
					valueStr += '['
					for tmpItem in v1:
						valueStr += str(tmpItem) + ' ' 
					valueStr += ']' + thisDelim
				else:
					valueStr += str(v1) + thisDelim
		valueStr += '\n' # weird
		
	'''
	print('headerStr:', headerStr)
	print('valueStr:', valueStr)
	'''
	
	#df0 = pd.concat(dfList0, axis=1) 
	#df0.to_excel('xxx0.xlsx', index=True, header=True)
	
	#
	# save as .csv
	#resultsFile = os.path.splitext(path)[0] + '_results.csv'
	resultsFile = savePath + '_batch_results.csv'
	print('DONE: saving results file:', resultsFile)
	with open(resultsFile, "w") as f:
		print(headerStr, file=f)
		print(valueStr, file=f)
	
	#
	# save raw data values
	#d = {}
	#dfList = [] # not used?
	dfList_branchLength = []
	dfList_euclideanDistance = []
	dfList_branchType = []
	dfList_tortuosity = []
	maxNumBranchLength = 0
	maxNumEuclideanDistance = 0
	maxNumBranchType = 0
	maxNumTortuosity = 0
	for idx, result in enumerate(resultsList): # resultsList can be across files
		# scale
		xVoxel = resultsList[idx]['xVoxel']
		yVoxel = resultsList[idx]['yVoxel']
		zVoxel = resultsList[idx]['zVoxel']
		
		# code is becoming sloppey, I use these to get tortuosity in k == 'euclideanDistance'
		thisBranchLength = None
		thisEuclideanDistance = None
		
		for k,v in result['data'].items():
			colName = k + '_' + str(idx)
			#print(idx, 'colName:', colName, type(v))
			#d[colName] = v
			#
			'''
			dfDict = {colName: v}
			df = pd.DataFrame(dfDict, index=range(len(v)))
			dfList.append(df)
			'''
			#
			if k == 'branchLength':
				#print('branchlength shape:', v.shape)
				if len(v) > maxNumBranchLength:
					maxNumBranchLength = len(v)
				#v = v * xVoxel # assuming x/y voxel is the same
				vScaled = np.multiply(v, xVoxel)
				thisBranchLength = vScaled # used in k == 'euclideanDistance' to get tortuosity
				dfDict = {colName: vScaled}
				df = pd.DataFrame(dfDict, index=range(len(v)))
				dfList_branchLength.append(df)
			if k == 'euclideanDistance':
				#print('branchlength shape:', v.shape)
				if len(v) > maxNumEuclideanDistance:
					maxNumEuclideanDistance = len(v)
				#v = v * xVoxel # assuming x/y voxel is the same
				vScaled = np.multiply(v, xVoxel)
				thisEuclideanDistance = vScaled # used in k == 'euclideanDistance' to get tortuosity
				dfDict = {colName: vScaled}
				df = pd.DataFrame(dfDict, index=range(len(v)))
				dfList_euclideanDistance.append(df)
				
				# assuming we have both branchLength and euclideanDistance, calculate the scale tortuosity
				# this is REALLY bad style !!!
				if len(v) > maxNumTortuosity:
					maxNumTortuosity = len(v)
				# this will print 'divide by zero encountered in true_divide' and value will become inf
				vScaled = np.divide(thisBranchLength, thisEuclideanDistance) # might fail on divide by 0
				dfDict = {colName: vScaled}
				df = pd.DataFrame(dfDict, index=range(len(v)))
				dfList_tortuosity.append(df)
				
			if k == 'branchType':
				#print('branchlength shape:', v.shape)
				if len(v) > maxNumBranchType:
					maxNumBranchType = len(v)
				dfDict = {colName: v} # no scaling
				df = pd.DataFrame(dfDict, index=range(len(v)))
				dfList_branchType.append(df)
			# this won't work, we need to calculate tortuosity from (scaled length / scaled euclidean)
			# hold off on this and calculate after we collate branchLength and euclideanDistance
			'''
			if k == 'tortuosity':
				if len(v) > maxNumTortuosity:
					maxNumTortuosity = len(v)
				#v = v * xVoxel # assuming x/y voxel is the same
				vScaled = np.multiply(v, xVoxel)
				dfDict = {colName: vScaled}
				df = pd.DataFrame(dfDict, index=range(len(v)))
				dfList_tortuosity.append(df)
			'''
	
	#
	# make 2d array of branch length to save numpy .npy and reopen for analysis
	# todo: take each of these 'for idx, result' and put into one loop
	print('=== making numpy arrays to save .npy')
	nResult = len(resultsList)
	npBranchLength = np.empty((maxNumBranchLength,nResult))  
	npBranchLength[:] = np.nan
	npEuclideanDistance = np.empty((maxNumEuclideanDistance,nResult))  
	npEuclideanDistance[:] = np.nan
	npBranchType = np.empty((maxNumBranchType,nResult))  
	npBranchType[:] = np.nan
	npTortuosity = np.empty((maxNumTortuosity,nResult))  
	npTortuosity[:] = np.nan
	for idx, result in enumerate(resultsList): # resultsList can be across files
		for k,v in result['data'].items():
			if k == 'branchLength':
				#v = v * xVoxel # assuming x/y voxel is the same
				vScaled = np.multiply(v, xVoxel)
				npBranchLength[:len(v),idx] = vScaled
			if k == 'euclideanDistance':
				vScaled = np.multiply(v, xVoxel)
				npEuclideanDistance[:len(v),idx] = vScaled
			if k == 'branchType':
				npBranchType[:len(v),idx] = v # no scaling
			if k == 'tortuosity':
				npTortuosity[:len(v),idx] = np.divide(npBranchLength[:,idx] / npEuclideanDistance[:,idx])

	print('   these should all have same shape ...')
	# branch length
	print('    npBranchLength.shape:', npBranchLength.shape)
	branchLengthFile = savePath + '_branchLength.npy'
	np.save(branchLengthFile, npBranchLength)

	# branch euclidean distance
	print('    npEuclideanDistance.shape:', npEuclideanDistance.shape)
	euclideanDistanceFile = savePath + '_euclideanDistance.npy'
	np.save(euclideanDistanceFile, npEuclideanDistance)
	# branch type

	print('    npBranchType.shape:', npBranchType.shape)
	branchTypeFile = savePath + '_branchType.npy'
	np.save(branchTypeFile, npBranchType)

	# tortuosity
	print('    npTortuosity.shape:', npTortuosity.shape)
	tortuosityFile = savePath + '_tortuosity.npy'
	np.save(tortuosityFile, npTortuosity)
			
	#
	# this works but replaced by single file with different sheets, see _summary.xlsx below
	# when processing one file
	#rawFile = os.path.splitext(path)[0] + '_results_raw.xlsx'
	'''
	rawFile = savePath + '_batch_results_raw.xlsx'
	print('    saving raw file:', rawFile)
	df = pd.concat(dfList, axis=1) 
	df.to_excel(rawFile, index=True, header=True)
	'''
	
	# save just branch length
	branchLengthFile = 'branchLength_summary.xlsx'
	df = pd.concat(dfList_branchLength, axis=1) 
	sheet_name = sheetName
	if os.path.isfile(branchLengthFile):
		mode = 'a'
	else:
		mode = 'w'
	with pd.ExcelWriter(branchLengthFile, mode=mode) as writer:
		df.to_excel(writer, sheet_name=sheet_name, index=True, header=True)
	
	# save just branch euclidean
	branchEuclideanFile = 'branchEuclidean_summary.xlsx'
	df = pd.concat(dfList_euclideanDistance, axis=1) 
	sheet_name = sheetName
	if os.path.isfile(branchEuclideanFile):
		mode = 'a'
	else:
		mode = 'w'
	with pd.ExcelWriter(branchEuclideanFile, mode=mode) as writer:
		df.to_excel(writer, sheet_name=sheet_name, index=True, header=True)
	
	# save just branch type
	branchTypeFile = 'branchType_summary.xlsx'
	df = pd.concat(dfList_branchType, axis=1) 
	sheet_name = sheetName
	if os.path.isfile(branchTypeFile):
		mode = 'a'
	else:
		mode = 'w'
	with pd.ExcelWriter(branchTypeFile, mode=mode) as writer:
		df.to_excel(writer, sheet_name=sheet_name, index=True, header=True)
	
	# save just tortuosity
	tortuosityFile = 'branchTortuosity_summary.xlsx'
	df = pd.concat(dfList_tortuosity, axis=1) 
	sheet_name = sheetName
	if os.path.isfile(tortuosityFile):
		mode = 'a'
	else:
		mode = 'w'
	with pd.ExcelWriter(tortuosityFile, mode=mode) as writer:
		df.to_excel(writer, sheet_name=sheet_name, index=True, header=True)
	