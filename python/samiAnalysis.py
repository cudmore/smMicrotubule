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
from datetime import datetime

import numpy as np
import pandas as pd
import tifffile

from IPython.display import display # to show pandas dataframe as a table

# statistical tests between analysis parameters
from scipy.stats import ttest_ind # used test sig difference between (length, euclidean distance)

# package for io, alternative is to just use tifffile
from aicsimageio import AICSImage, omeTifWriter							

# function for core algorithm
from aicssegmentation.core.vessel import filament_3d_wrapper
from aicssegmentation.core.pre_processing_utils import intensity_normalization, edge_preserving_smoothing_3d, image_smoothing_gaussian_3d
from skimage.morphology import remove_small_objects	

from skimage import morphology # to convert 3d mask to 1-pixel skeleton, using morphology.skeletonize_3d
import skan

# my code
from my_suggest_normalization_param import my_suggest_normalization_param # clone of aics segmentation
from readVoxelSize import readVoxelSize # to read x/y/z scale from Fiji save .tif

def myAnalyzeSkeleton(out=None, maskPath=None, imagePath=None):
	"""
	out: numpy array with 1-pixel skeleton
	maskPath : full path to _dvMask.tif file (can include appended _0.tif
	
	returns:
	    dict of results, 3d skeleton
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
	
	print('    === myAnalyzeSkeleton() maskData.shape:', maskData.shape)
	
	# make a 1-pixel skeleton from volume mask (similar to Fiji Skeletonize)
	mySkeleton = morphology.skeletonize_3d(maskData)

	# analyze the skeleton (similar to Fiji Analyze Skeleton)
	mySkanSkel = skan.Skeleton(mySkeleton)

	# look at the results
	branch_data = skan.summarize(mySkanSkel) # branch_data is a pandas dataframe
	nBranches = branch_data.shape[0]
	
	#
	# convert everything to nump arrays
	branchType = branch_data['branch-type'].to_numpy()
	branchDistance = branch_data['branch-distance'].to_numpy()
	euclideanDistance = branch_data['euclidean-distance'].to_numpy()
	# don't do tortuosity here, we need to scale to um/pixel in x/y/z

	#
	# scale
	branchDistance = np.multiply(branchDistance, xVoxel)
	euclideanDistance = np.multiply(euclideanDistance, xVoxel)
	# this will print 'divide by zero encountered in true_divide' and value will become inf
	tortuosity = np.divide(branchDistance, euclideanDistance) # might fail on divide by 0
	
	#
	# organize a return dicitonary
	retDict = OrderedDict()

	retDict['data'] = OrderedDict()
	retDict['data']['branchType'] = branchType
	retDict['data']['branchLength'] = branchDistance
	retDict['data']['euclideanDistance'] = euclideanDistance
	retDict['data']['tortuosity'] = tortuosity

	# todo: search for 0 values in (branchDistance, euclideanDistance)
	
	return retDict, mySkeleton # returning mySkeleton so we can save it
	
def run(path, f3_param=[[1, 0.01]], minArea=20, saveNumber=0):
	"""
	use aicssegmentation to pre-process raw data and then make/save a 3D mask
	"""
	print('    === run() path:', path, 'f3_param:', f3_param, 'minArea:', minArea, 'saveNumber:', saveNumber)
		
	# load the data
	reader = AICSImage(path) 
	IMG = reader.data.astype(np.float32)

	# load x/y/z voxel size (assumes .tif was saved with Fiji
	xVoxel, yVoxel, zVoxel = readVoxelSize(path)
	print('    file:', os.path.basename(path), 'has shape:', IMG.shape, 'xVoxel:', xVoxel, 'yVoxel:', yVoxel, 'zVoxel:', zVoxel)

	structure_channel = 0
	struct_img0 = IMG[0,structure_channel,:,:,:].copy()

	# give us a guess for our intensity_scaling_param parameters
	low_ratio, high_ratio = my_suggest_normalization_param(struct_img0)

	#intensity_scaling_param = [0.0, 22.5]
	intensity_scaling_param = [low_ratio, high_ratio]
	print('    === intensity_normalization() intensity_scaling_param:', intensity_scaling_param)
	
	# intensity normalization
	print('    === calling intensity_normalization()')
	struct_img = intensity_normalization(struct_img0, scaling_param=intensity_scaling_param)

	# smoothing with edge preserving smoothing 
	print('    === calling edge_preserving_smoothing_3d()')
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
	print('    === calling filament_3d_wrapper() f3_param:', f3_param)
	bw = filament_3d_wrapper(structure_img_smooth, f3_param)
		
	#
	#minArea = 20 # from recipe
	print('   === calling remove_small_objects() minArea:', minArea)
	seg = remove_small_objects(bw>0, min_size=minArea, connectivity=1, in_place=False)

	#
	# save original file again (with saveNumber
	saveNumberStr = ''
	if saveNumber>1:
		saveNumberStr = '_' + str(saveNumber)
		
	#
	# save _dvMask.tif
	seg = seg >0
	out=seg.astype(np.uint8)
	out[out>0]=255	
	#
	maskPath = os.path.splitext(path)[0] + '_dvMask' + saveNumberStr + '.tif'
	print('   === saving 3D mask [WILL FAIL IF FILE EXISTS] as maskPath:', maskPath)
	try:
		writer = omeTifWriter.OmeTifWriter(maskPath)
		writer.save(out)
	except(OSError) as e:
		print('    ******** ERROR: file already exists, did not resave, maskPath:', maskPath)
		
	#
	# ################
	# analyze skeleton, take a 3d mask and analyze as a 1-pixel skeleton
	retDict0, mySkeleton = myAnalyzeSkeleton(out=out, imagePath=path) 
	# ################
	retDict = OrderedDict()
	retDict['analysisDate'] = datetime.today().strftime('%Y%m%d')
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

	#
	retDict.update( retDict0 ) # this has 'keys' that are lists of ('len3d', 'eLen', 'branchType')

	# save 1-pixel skeleton: mySkeleton
	# save _dvSkel
	skelPath = os.path.splitext(path)[0] + '_dvSkel' + saveNumberStr + '.tif'
	print('    === saving 3D skel [WILL FAIL IF FILE EXISTS] as maskPath:', skelPath)
	try:
		writer = omeTifWriter.OmeTifWriter(skelPath)
		writer.save(mySkeleton)
	except(OSError) as e:
		print('    ******** ERROR: file already exists, did not resave, skelPath:', skelPath)
			
	return retDict
	
def parseBatchFile(path):
	"""
	given a file with one path per line, return a list paths
	"""
	with open(path) as f:
		content = f.readlines()
	# remove whitespace characters like `\n` at the end of each line
	content = [x.strip() for x in content] 
	return content
	
if __name__ == '__main__':
		
	#
	# parse command line arguments
	path = None
	pathList = None
	savePath = ''
	doBatchFile = None
	sheetName = None
	genotype = ''
	sex = ''
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
			#
			batchFileName = os.path.basename(batchFilePath) #wt-female.txt
			batchFileName = batchFileName.split('.')[0]
			genotype = batchFileName.split('-')[0]
			sex = batchFileName.split('-')[1]
			print('genotype:', genotype)
			print('sex:', sex)
			
	if pathList is None:
		pathList = [path]
	
	'''
	folderPath = '/Users/cudmore/box/data/sami/200414'
	pathList = []
	for dirpath, dirnames, files in os.walk(folderPath):
		for name in files:
			if name.endswith('_ch2.tif'):
				filePath = os.path.join(dirpath, name)
				#x,y,z = readVoxelSize(filePath, verbose=False)
				print('    ', filePath)
				pathList.append(filePath)
				#myConvert(path=filePath, imp=None)
	'''

	# put all results in pandas df, like ('filename', 'path', 'genotype', 'sex', 'branchType', 'len3d', 'euclideanDist', 'tort')
	column_names = ['myCellNumber', 'filename', 'genotype', 'sex', 'branchType', 'len3d', 'euclideanDist', 'tortuosity', 'path']
	dfMaster = pd.DataFrame(columns = column_names)
	
	resultsList = []
	#pathList = [path]
	for myCellNumber, path in enumerate(pathList):	
	
		print('\n***')
		print('*** path:', path)
		#print('    fileName:', fileName)
		#print('    fileNameNoExtension:', fileNameNoExtension)
	
		#1
		saveNumber = 1
		f3_param=[[1, 0.01]]
		f3_param=[[0.5, 0.005]]
		minArea=20
		# run the masking and skeletonization and append to list
		retDict = run(path, f3_param, minArea=minArea, saveNumber=saveNumber)
		resultsList.append(retDict)
		
		# this allows us to use different detection parameters and save _mask_<n> and _skel_<n>
		'''
		#2
		saveNumber = 2
		f3_param=[[1, 0.01]]
		#f3_param=[1, 0.01]
		#minArea=40
		minArea=20
		# run the masking and skeletonization and append to list
		retDict = run(path, f3_param, minArea=minArea, saveNumber=saveNumber)
		resultsList.append(retDict)
		'''
		
		#
		# put these in df
		df2 = pd.DataFrame(columns = column_names)
		# lists
		df2['branchType'] = retDict['data']['branchType']
		df2['len3d'] = retDict['data']['branchLength']
		df2['euclideanDist'] = retDict['data']['euclideanDistance']
		df2['tortuosity'] = retDict['data']['tortuosity']
		# single values (all same)
		df2['myCellNumber'] = myCellNumber # fills in all rows
		df2['filename'] = os.path.basename(path) # fills in all rows
		df2['genotype'] = genotype # fills in all rows, from name of batch=batFilePath, e.g. wt-female.txt
		df2['sex'] = sex # fills in all rows
		df2['path'] = path # fills in all rows
		
		# append to master
		dfMaster = dfMaster.append(df2, ignore_index = True) 
		
	# save master, across cells in this batch/cond, generally (genotype, sex)
	#pandasPath = '/Users/cudmore/Desktop/masterDF.csv'
	pandasPath = savePath + '_results.csv'
	print('saving master pandas dataframe pandasPath:', pandasPath)
	dfMaster.to_csv(pandasPath)
	
	#
	# do some quick stats on our master dfMaster
	
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
		# result is for one cell
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
		
	#
	# save as .csv from the flat output
	#resultsFile = os.path.splitext(path)[0] + '_results.csv'
	resultsFile = savePath + '_batch_results.csv'
	print('DONE: saving results file:', resultsFile)
	with open(resultsFile, "w") as f:
		print(headerStr, file=f)
		print(valueStr, file=f)
	
	##
	##
	## see old_code.py
	##
	##
	