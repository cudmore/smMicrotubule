"""
Robert Cudmore
20200713, 20200330

parallel version of samiAnalysis

wt-female (n=22)
	53 sec (cpu-2)
	53 sec (cpu-1), Nominal improvement
	in series takes 124 seconds (2.3 x speed)

"""

import os, sys, time

import numpy as np
import pandas as pd

import multiprocessing as mp

from samiAnalysis import myRun
from samiAnalysis import parseBatchFile
from samiAnalysis import saveFlatResults
from samiAnalysis import getColumnNames
from samiAnalysis import oneCellDataFrame

if __name__ == '__main__':
		
	startTime = time.time()
	
	#
	# parse command line arguments
	path = None
	pathList = None
	savePath = ''
	#doBatchFile = None
	#sheetName = None
	genotype = ''
	sex = ''
	for i, arg in enumerate(sys.argv):
		if i == 0:
			# my .py filename (not full path)
			continue
		#print('i:', i, 'arg:', arg)
		if arg.startswith('file='):
			#doBatchFile = False
			tmp, filePath = arg.split('=')
			path = filePath
			savePath = os.path.splitext(path)[0]
			fileNameNoExtension = os.path.splitext(os.path.basename(path))[0]
			#sheetName = fileNameNoExtension
		
		elif arg.startswith('batch='):
			#doBatchFile = True
			# open specified file and make a list of path
			#os.path.splitext(path)[0]
			tmp, batchFilePath = arg.split('=')
			pathList = parseBatchFile(batchFilePath) ## LOCAL FUNCTION
			savePath = os.path.splitext(batchFilePath)[0] # get path with filename without extension
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
	
	# put all results in pandas df, like ('filename', 'path', 'genotype', 'sex', 'branchType', 'len3d', 'euclideanDist', 'tort')
	#column_names = ['myCellNumber', 'filename', 'genotype', 'sex', 'xVoxel', 'yVoxel', 'zVoxel', 'branchType', 'len3d', 'euclideanDist', 'tortuosity', 'path']
	column_names = getColumnNames()
	dfMaster = pd.DataFrame(columns = column_names)
	
	resultsList = []
	#pathList = [path]
	nPathList = len(pathList)

	'''
	saveBase = '/Users/cudmore/Desktop/samiVolume2'
	saveBase = '/Users/cudmore/Desktop/samiVolume3'
	'''
	
	from gAnalysisPath import gSami_Params
	saveBase = gSami_Params['gAnalysisPath']
	
	if not os.path.isdir(saveBase):
		print('samiAnalysisParallel is making output dir:', saveBase)
		os.mkdir(saveBase)
	
	f3_param = gSami_Params['samiAnalysis_params']['f3_param']
	minArea = gSami_Params['samiAnalysis_params']['minArea']
	
	#f3_param=[[1, 0.01]] #saveNumber = 1
	'''
	f3_param=[[0.5, 0.005]] # samiVolume3
	minArea=20
	'''
	#
	# parallel
	cpuCount = mp.cpu_count()
	cpuCount -= 2
	pool = mp.Pool(processes=cpuCount)
	results = []

	for myCellNumber, path in enumerate(pathList):	
		#results = [pool.apply_async(myRun, args=(file,myIdx+1)) for myIdx, file in enumerate(filenames)]
		args=(path, myCellNumber, genotype, sex, saveBase, f3_param, minArea) #, saveNumber)
		oneAnswer = pool.apply_async(myRun, args=args)
		results.append(oneAnswer)
		
	# parallel, get results
	for myCellNumber, result in enumerate(results):
		path = pathList[myCellNumber]
		
		retDict = result.get() ### get results from parallel

		if retDict is None:
			print('\n\n\n                            ERROR: File not found idx:', idx, '\n\n\n')
			continue

		resultsList.append(retDict)

		#
		df2 = oneCellDataFrame(retDict)
		
		# append to master
		dfMaster = dfMaster.append(df2, ignore_index = True) 
		
	# save master, across cells in this batch/cond, generally (genotype, sex)
	#pandasPath = '/Users/cudmore/Desktop/masterDF.csv'
	#pandasPath = savePath + '_results_p.csv'
	pandasPath = os.path.join(saveBase, batchFileName + '_results.csv')
	print('saving master pandas dataframe pandasPath:', pandasPath)
	dfMaster.to_csv(pandasPath)
	
	#
	#
	#flatResultPath = savePath + '_batch_results_p.csv'
	flatResultPath = os.path.join(saveBase, batchFileName + '_results_summary.csv')
	saveFlatResults(resultsList, flatResultPath)
		
	stopTime = time.time()
	print('  \nsamiAnalysisParallel took', round(stopTime-startTime,2), 'seconds for genotype:', genotype, 'sex:', sex, 'n=', len(results), '\n\n')
	##
	##
	## see old_code.py
	##
	##
	