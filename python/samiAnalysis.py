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
		#cd /Users/cudmore/box/data/sami/
		cd /Users/cudmore/Sites/smMicrotubule/python
	
	run it:
		#python samiAnalysis.py
		python samiAnalysis.py batch=../analysis/wt-female.txt 
	
	remember, I was using original recipe from jupyter notebook in:
		notebooks/playground_filament3d.ipynb
"""

import os, sys, time
from collections import OrderedDict
from datetime import datetime

import numpy as np
import pandas as pd
#import tifffile

import multiprocessing as mp

#from IPython.display import display # to show pandas dataframe as a table

# statistical tests between analysis parameters
#from scipy.stats import ttest_ind # used test sig difference between (length, euclidean distance)

# package for io, alternative is to just use tifffile
#from aicsimageio import AICSImage, omeTifWriter							

# function for core algorithm
from aicssegmentation.core.vessel import filament_3d_wrapper
from aicssegmentation.core.pre_processing_utils import edge_preserving_smoothing_3d, image_smoothing_gaussian_3d
from skimage.morphology import remove_small_objects	

from skimage import morphology # to convert 3d mask to 1-pixel skeleton, using morphology.skeletonize_3d
import skan

# my code
# I rewrote (my_suggest_normalization_param, my_intensity_normalization) so they do not print loads of crap
from my_suggest_normalization_param import my_suggest_normalization_param # clone of aics segmentation
from my_intensity_normalization import my_intensity_normalization

from readVoxelSize import readVoxelSize # to read x/y/z scale from Fiji save .tif

import bimpy.util.bTiffFile

def myAnalyzeSkeleton(out=None, maskPath=None, imagePath=None, saveBase=None, verbose=False):
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
		#maskData = tifffile.imread(maskPath)
		maskData, maskHeader = bimpy.util.bTiffFile.imread(maskPath)
		
	# was used by shape_index
	#imageData = tifffile.imread(imagePath)
	
	if verbose: print('    === myAnalyzeSkeleton() maskData.shape:', maskData.shape)
	
	##
	##
	# make a 1-pixel skeleton from volume mask (similar to Fiji Skeletonize)
	mySkeleton = morphology.skeletonize_3d(maskData)
	##
	##
	
	##
	##
	# analyze the skeleton (similar to Fiji Analyze Skeleton)
	## BE SURE TO INCLUDE VOXEL SIZE HERE !!!!!! 20200503
	mySpacing_ = (zVoxel, xVoxel, yVoxel)
	mySkanSkel = skan.Skeleton(mySkeleton, spacing=mySpacing_)
	##
	##
	
	# look at the results
	branch_data = skan.summarize(mySkanSkel) # branch_data is a pandas dataframe
	nBranches = branch_data.shape[0]
	
	# working on eroded/ring density 20200501
	# save entire skan analysis as csv
	# 20200713, was this
	'''
	tmpFolder, tmpFileName = os.path.split(imagePath)
	tmpFileNameNoExtension, tmpExtension = tmpFileName.split('.')
	saveSkelPath = os.path.join(tmpFolder, tmpFileNameNoExtension + '_skel.csv')
	if verbose: print('saving skan results to saveSkelPath:', saveSkelPath)
	'''
	
	saveSkelPath = saveBase + '_skel.csv'
	print('    myAnalyzeSkeleton() saving saveSkelPath:', saveSkelPath)
	
	branch_data.to_csv(saveSkelPath)
	
	#
	# convert everything to nump arrays
	branchType = branch_data['branch-type'].to_numpy()
	branchDistance = branch_data['branch-distance'].to_numpy()
	euclideanDistance = branch_data['euclidean-distance'].to_numpy()
	# don't do tortuosity here, we need to scale to um/pixel in x/y/z

	#
	# scale
	# 20200503 TRYING TO DO THIS WHEN CALLING skan.Skeleton(mySkeleton, spacing=mySpacing_) !!!!!!!!!!!!!!!!!!!!!!!
	'''
	branchDistance = np.multiply(branchDistance, xVoxel)
	euclideanDistance = np.multiply(euclideanDistance, xVoxel)
	'''
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
	
	# 20200503 working on samiPostAnalysis density
	# we need all the src/dst point so we can quickly determine if they are in mask (full, eroded, ring)
	# 'image-coord-src-0', 'image-coord-src-1', 'image-coord-src-2', 'image-coord-dst-0', 'image-coord-dst-1', 'image-coord-dst-2'
	image_coord_src_0 = branch_data['image-coord-src-0'].to_numpy()
	image_coord_src_1 = branch_data['image-coord-src-1'].to_numpy()
	image_coord_src_2 = branch_data['image-coord-src-2'].to_numpy()
	image_coord_dst_0 = branch_data['image-coord-dst-0'].to_numpy()
	image_coord_dst_1 = branch_data['image-coord-dst-1'].to_numpy()
	image_coord_dst_2 = branch_data['image-coord-dst-2'].to_numpy()
	retDict['data']['image_coord_src_0'] = image_coord_src_0
	retDict['data']['image_coord_src_1'] = image_coord_src_1
	retDict['data']['image_coord_src_2'] = image_coord_src_2
	retDict['data']['image_coord_dst_0'] = image_coord_dst_0
	retDict['data']['image_coord_dst_1'] = image_coord_dst_1
	retDict['data']['image_coord_dst_2'] = image_coord_dst_2
	
	return retDict, mySkeleton # returning mySkeleton so we can save it
	
def myRun(path, myCellNumber, genotype, sex, saveBase = '/Users/cudmore/Desktop/samiVolume2', f3_param=[1, 0.01], minArea=20, verbose=False): #, saveNumber=0):
	"""
	use aicssegmentation to pre-process raw data and then make/save a 3D mask
	
	path: path to raw tif, the _ch2.tif
	myCellNumber: cell number from batch file (per genotype/sex), NOT unique across different (genotype, sex)
	saveBase: full path to folder to save to (e.g. /Users/cudmore/Desktop/samiVolume2), must exist
	...
	saveNumber: not used
	"""
		
	print('  === myRun() path:', path, 'saveBase:', saveBase, 'f3_param:', f3_param, 'minArea:', minArea) #, 'saveNumber:', saveNumber)

	#20200608
	#saveBase = '/Users/cudmore/Desktop/samiVolume2'
	tmpPath, tmpFileName = os.path.split(path)
	tmpPath = tmpPath.replace('../data/', '')
	tmpFileNameNoExtension, tmpExtension = tmpFileName.split('.')
	saveBase = os.path.join(saveBase, tmpPath)
	if not os.path.isdir(saveBase):
		print('    making output dir:', saveBase)
		os.makedirs(saveBase)
	#saveBase = os.path.join(saveBase, tmpPath, tmpFileNameNoExtension)
	saveBase = os.path.join(saveBase, tmpFileNameNoExtension)
	#saveBase = os.path.join(saveBase, os.path.splitext(path)[0].replace('../data/', ''))
	if verbose: print('    saveBase:', saveBase)
			
	if not os.path.isfile(path):
		print('ERROR: myRun() did not find file:', path)
		return None
	
	# load the data
	IMG, tifHeader = bimpy.util.bTiffFile.imread(path)
	
	saveDataPath = saveBase + '.tif'
	print('   === saving raw data to saveDataPath:', saveDataPath)
	bimpy.util.bTiffFile.imsave(saveDataPath, IMG, tifHeader=tifHeader, overwriteExisting=True)

	IMG = IMG.astype(np.float32)

	# channel 1: load then save (do nothing else with it)
	channelOnePath = path.replace('_ch2.tif', '_ch1.tif')
	channelOneData, channelOneTiffHeader = bimpy.util.bTiffFile.imread(channelOnePath)
	saveChannelOnePath = saveBase.replace('_ch2', '_ch1.tif')
	bimpy.util.bTiffFile.imsave(saveChannelOnePath, channelOneData, tifHeader=tifHeader, overwriteExisting=True)
	
	# load x/y/z voxel size (assumes .tif was saved with Fiji
	xVoxel, yVoxel, zVoxel = readVoxelSize(path)
	print('    file:', os.path.basename(path), 'has shape:', IMG.shape, 'xVoxel:', xVoxel, 'yVoxel:', yVoxel, 'zVoxel:', zVoxel)

	# give us a guess for our intensity_scaling_param parameters
	#low_ratio, high_ratio = my_suggest_normalization_param(struct_img0)
	low_ratio, high_ratio = my_suggest_normalization_param(IMG)

	#intensity_scaling_param = [0.0, 22.5]
	intensity_scaling_param = [low_ratio, high_ratio]
	if verbose: print('    === my_intensity_normalization() intensity_scaling_param:', intensity_scaling_param)
	
	# intensity normalization
	if verbose: print('    === calling my_intensity_normalization()')
	#struct_img = my_intensity_normalization(struct_img0, scaling_param=intensity_scaling_param)
	struct_img = my_intensity_normalization(IMG, scaling_param=intensity_scaling_param)

	# smoothing with edge preserving smoothing 
	if verbose: print('    === calling edge_preserving_smoothing_3d()')
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
	if verbose: print('    === calling filament_3d_wrapper() f3_param:', f3_param)
	bw = filament_3d_wrapper(structure_img_smooth, [f3_param]) # f3_param is a list of a list
		
	#
	#minArea = 20 # from recipe
	if verbose: print('   === calling remove_small_objects() minArea:', minArea)
	seg = remove_small_objects(bw>0, min_size=minArea, connectivity=1, in_place=False)

	#
	# save original file again (with saveNumber
	saveNumberStr = ''
	#if saveNumber>1:
	#	saveNumberStr = '_' + str(saveNumber)
		
	#
	# save _dvMask.tif
	seg = seg >0
	out=seg.astype(np.uint8)
	out[out>0]=255	
	#
	#maskPath = os.path.splitext(path)[0] + '_dvMask' + saveNumberStr + '.tif'
	maskPath = saveBase + '_dvMask' + saveNumberStr + '.tif'
	print('   === saving 3D mask as maskPath:', maskPath)
	#tifffile.imsave(maskPath, out)
	bimpy.util.bTiffFile.imsave(maskPath, out, tifHeader=tifHeader, overwriteExisting=True)
	
	'''
	try:
		writer = omeTifWriter.OmeTifWriter(maskPath)
		writer.save(out)
	except(OSError) as e:
		print('    ******** ERROR: file already exists, did not resave, maskPath:', maskPath)
	'''
		
	#
	# ################
	# analyze skeleton, take a 3d mask and analyze as a 1-pixel skeleton
	retDict0, mySkeleton = myAnalyzeSkeleton(out=out, imagePath=path, saveBase=saveBase) 
	# ################
	retDict = OrderedDict()
	retDict['analysisDate'] = datetime.today().strftime('%Y%m%d')
	#
	retDict['saveBase'] = saveBase
	retDict['myCellNumber'] = myCellNumber
	retDict['genotype'] = genotype
	retDict['sex'] = sex
	#
	retDict['path'] = path # 20200713 working on parallel
	retDict['tifPath'] = path
	retDict['maskPath'] = maskPath
	retDict['tifFile'] = os.path.basename(path)
	retDict['xVoxel'] = xVoxel
	retDict['yVoxel'] = yVoxel
	retDict['zVoxel'] = zVoxel
	#
	retDict['params'] = OrderedDict()
	#retDict['params']['saveNumber'] = saveNumber
	retDict['params']['intensity_scaling_param'] = intensity_scaling_param # calculated in my_suggest_normalization_param
	#retDict['params']['f3_param'] = f3_param[0] # cludge, not sure where to put this. f3_param is a list of list but screws up my .csv output !!!
	retDict['params']['f3_param'] = f3_param # cludge, not sure where to put this. f3_param is a list of list but screws up my .csv output !!!
	retDict['params']['minArea'] = minArea

	#
	retDict.update( retDict0 ) # this has 'keys' that are lists of ('len3d', 'eLen', 'branchType')

	# save 1-pixel skeleton: mySkeleton
	# save _dvSkel
	#skelPath = os.path.splitext(path)[0] + '_dvSkel' + saveNumberStr + '.tif'
	skelPath = saveBase + '_dvSkel' + saveNumberStr + '.tif'
	print('    === saving 3D skel as maskPath:', skelPath)
	#tifffile.imsave(skelPath, mySkeleton)
	bimpy.util.bTiffFile.imsave(skelPath, mySkeleton, tifHeader=tifHeader, overwriteExisting=True)

	return retDict
	
def saveFlatResults(resultsList, flatResultPath):
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
	#resultsFile = savePath + '_batch_results_p.csv'
	print('DONE: saving results file:', flatResultPath)
	with open(flatResultPath, "w") as f:
		print(headerStr, file=f)
		print(valueStr, file=f)

def getColumnNames():
	theRet = ['myCellNumber', 'saveBase', 'filename', 'genotype', 'sex', 'xVoxel', 'yVoxel', 'zVoxel', 'branchType', 'len3d', 'euclideanDist', 'tortuosity', 'path']
	return theRet
	
def oneCellDataFrame(retDict):
	column_names = getColumnNames()
	#
	# put these in df
	df2 = pd.DataFrame(columns = column_names)
	# lists
	df2['branchType'] = retDict['data']['branchType']
	df2['len3d'] = retDict['data']['branchLength']
	df2['euclideanDist'] = retDict['data']['euclideanDistance']
	df2['tortuosity'] = retDict['data']['tortuosity']
	# single values (all same)
	df2['myCellNumber'] = retDict['myCellNumber'] # fills in all rows
	df2['filename'] = os.path.basename(retDict['path']) # fills in all rows
	df2['genotype'] = retDict['genotype'] # fills in all rows, from name of batch=batFilePath, e.g. wt-female.txt
	df2['sex'] = retDict['sex'] # fills in all rows
	df2['xVoxel'] = retDict['xVoxel'] # fills in all rows
	df2['yVoxel'] = retDict['yVoxel'] # fills in all rows
	df2['zVoxel'] = retDict['zVoxel'] # fills in all rows
	df2['path'] = retDict['path'] # fills in all rows
	df2['saveBase'] = retDict['saveBase'] # fills in all rows

	# 20200503 working on density analysis
	# ['data']['image_coord_src_0']
	df2['image_coord_src_0'] = retDict['data']['image_coord_src_0']
	df2['image_coord_src_1'] = retDict['data']['image_coord_src_1']
	df2['image_coord_src_2'] = retDict['data']['image_coord_src_2']
	df2['image_coord_dst_0'] = retDict['data']['image_coord_dst_0']
	df2['image_coord_dst_1'] = retDict['data']['image_coord_dst_1']
	df2['image_coord_dst_2'] = retDict['data']['image_coord_dst_2']

	return df2
	
def parseBatchFile(path):
	"""
	given a file with one path per line, return a list paths
	"""
	retList = []
	with open(path) as f:
		content = f.readlines()
	# remove whitespace characters like `\n` at the end of each line
	content = [x.strip() for x in content] 
	for line in content:
		if line.startswith('#'):
			continue
		else:
			retList.append(line)
	return retList
	
if __name__ == '__main__':
		
	startTime = time.time()

	#
	# parse command line arguments
	path = None
	pathList = None
	savePath = ''
	doBatchFile = None
	#sheetName = None
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
			batchFileName = fileNameNoExtension
			#sheetName = fileNameNoExtension
		elif arg.startswith('batch='):
			doBatchFile = True
			# open specified file and make a list of path
			#os.path.splitext(path)[0]
			tmp, batchFilePath = arg.split('=')
			pathList = parseBatchFile(batchFilePath) ## LOCAL FUNCTION
			savePath = os.path.splitext(batchFilePath)[0]
			batchFileNoExtension = os.path.splitext(os.path.basename(batchFilePath))[0]
			#sheetName = batchFileNoExtension
			#
			batchFileName = os.path.basename(batchFilePath) #wt-female.txt
			batchFileName = batchFileName.split('.')[0]
			genotype = batchFileName.split('-')[0]
			sex = batchFileName.split('-')[1]
			print('    genotype:', genotype)
			print('    sex:', sex)
			
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
	column_names = ['myCellNumber', 'filename', 'genotype', 'sex', 'xVoxel', 'yVoxel', 'zVoxel', 'branchType', 'len3d', 'euclideanDist', 'tortuosity', 'path']
	dfMaster = pd.DataFrame(columns = column_names)
	
	resultsList = []
	#pathList = [path]
	nPathList = len(pathList)
	for myCellNumber, path in enumerate(pathList):	
	
		print('\n*** myCellNumber:', myCellNumber, 'of', nPathList)
		print('*** path:', path)
		#print('    fileName:', fileName)
		#print('    fileNameNoExtension:', fileNameNoExtension)
	
		#
		# series
		#1
		#saveNumber = 1
		saveBase = '/Users/cudmore/Desktop/samiVolume1'
		if not os.path.isdir(saveBase):
			print('samiAnalysis is making output dir:', saveBase)
			os.mkdir(saveBase)
		f3_param=[[1, 0.01]]
		f3_param=[[0.5, 0.005]]
		minArea=20
		# run the masking and skeletonization and append to list
		args=(path, myCellNumber, genotype, sex, saveBase, f3_param, minArea) #, saveNumber)
		# was working
		#retDict = myRun(path, f3_param, minArea=minArea)
		retDict = myRun(args)
		if retDict is None:
			print('\n\n\n                            ERROR: File not found path:', path, '\n\n\n')
			continue
		
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
		retDict = myRun(path, f3_param, minArea=minArea, saveNumber=saveNumber)
		resultsList.append(retDict)
		'''
		
		#
		df2 = oneCellDataFrame(retDict)

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
		df2['xVoxel'] = retDict['xVoxel'] # fills in all rows
		df2['yVoxel'] = retDict['yVoxel'] # fills in all rows
		df2['zVoxel'] = retDict['zVoxel'] # fills in all rows
		df2['path'] = path # fills in all rows
		
		# 20200503 working on density analysis
		# ['data']['image_coord_src_0']
		df2['image_coord_src_0'] = retDict['data']['image_coord_src_0']
		df2['image_coord_src_1'] = retDict['data']['image_coord_src_1']
		df2['image_coord_src_2'] = retDict['data']['image_coord_src_2']
		df2['image_coord_dst_0'] = retDict['data']['image_coord_dst_0']
		df2['image_coord_dst_1'] = retDict['data']['image_coord_dst_1']
		df2['image_coord_dst_2'] = retDict['data']['image_coord_dst_2']
		'''
		
		# append to master
		dfMaster = dfMaster.append(df2, ignore_index = True) 
		
	# save master, across cells in this batch/cond, generally (genotype, sex)
	#pandasPath = '/Users/cudmore/Desktop/masterDF.csv'
	#pandasPath = savePath + '_results.csv'
	pandasPath = os.path.join(saveBase, batchFileName + '_results.csv')
	print('saving master pandas dataframe pandasPath:', pandasPath)
	dfMaster.to_csv(pandasPath)
	
	#
	# do some quick stats on our master dfMaster
	
	#
	#
	#flatResultPath = savePath + '_batch_results.csv'
	flatResultPath = os.path.join(saveBase, batchFileName + '_results_summary.csv')
	saveFlatResults(resultsList, flatResultPath)

	"""
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
	"""
	
	stopTime = time.time()
	print('  took', round(stopTime-startTime,2), 'seconds')

	##
	##
	## see old_code.py
	##
	##
	