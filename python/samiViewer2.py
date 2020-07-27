"""
20200604

given one batch file list of all cells, plot (raw, eroded, ring) in a grid

"""

import os, sys, math

import numpy as np

import tifffile

import napari

import bimpy

def samiDaskNapari(gAnalysisPath, batchFilePath):

	print('samiDaskNapari() batchFilePath:', batchFilePath)
	
	# path to _dvMask and dvSekl
	#smMicrotubulePath = '/Users/cudmore/Sites/smMicrotubule/data'
	
	# path to (full, eroded, ring) masks
	#dataPath = '/Users/cudmore/Desktop/samiVolume' # old
	#dataPath = '/Users/cudmore/Desktop/samiVolume2' # new 20200605, reduced median and added initial 1 iteration erosion
	
	dataPath = gAnalysisPath

	batchFile = bimpy.util.getFileListFromFile(batchFilePath)
	
	maxRows = 0
	maxCols = 0
	maxSlices = 0
	
	dataList = []
	fullMaskList = []
	erodedList = []
	ringList = []
	#
	dvMaskList = []
	dvSkelList = []
	rawCh1List = []
	
	numFiles = len(batchFile)
	
	print('  loading', numFiles, 'files')
	
	for idx, file in enumerate(batchFile):
		file = file.replace('../data/', '')
		
		fileNameNoExtension, _ = os.path.splitext(file)
		
		'''
		isBad = False
		if batchFilePath.endswith('ko-female.txt') and idx in [16]:
			print('  bad', idx, file)
			isBad = True
		'''
			
		#
		# for (_ch1.tif, _dvMask, _dvSkel)
		#filePathNoExtension = os.path.join(smMicrotubulePath, fileNameNoExtension)
		filePathNoExtension = os.path.join(dataPath, fileNameNoExtension)
		dvMaskPath = filePathNoExtension + '_dvMask.tif'
		dvSkelPath = filePathNoExtension + '_dvSkel.tif'
		rawCh1Path = filePathNoExtension.replace('_ch2', '')  + '_ch1.tif'
		
		dvMaskData = tifffile.imread(dvMaskPath)
		dvSkelData = tifffile.imread(dvSkelPath)
		rawCh1Data = tifffile.imread(rawCh1Path)

		# 20200707
		'''
		if isBad:
			rawCh1Data[:] = 0
		'''
		
		dvMaskList.append(dvMaskData)
		dvSkelList.append(dvSkelData)
		rawCh1List.append(rawCh1Data)
		
		#
		# for (image, full mask, eroded mask, ring mask) of _ch2
		filePathNoExtension = os.path.join(dataPath, fileNameNoExtension)
		
		rawPath = filePathNoExtension + '.tif'
		fullMask = filePathNoExtension + '_filledHolesMask.tif'
		erodedPath = filePathNoExtension + '_erodedMask.tif'
		ringPath = filePathNoExtension + '_ringMask.tif'
		
		if not os.path.isfile(rawPath):
			print('  error: samiDaskNapari() did not find file:', rawPath)
			continue
		
		imageData = tifffile.imread(rawPath)
		fullMaskData = tifffile.imread(fullMask)
		erodedData = tifffile.imread(erodedPath)
		ringData = tifffile.imread(ringPath)
				
		# 20200707
		'''
		if isBad:
			imageData[:] = 0
		'''
		
		dataList.append(imageData)
		fullMaskList.append(fullMaskData)
		erodedList.append(erodedData)
		ringList.append(ringData)
		
		#print('  ', filePathNoExtension)
		#print('  ', imageData.shape, erodedData.shape, ringData.shape)

		numSlices = imageData.shape[0]
		if numSlices> maxSlices:
			maxSlices = numSlices
		numRows = imageData.shape[1]
		if numRows> maxRows:
			maxRows = numRows
		numCols = imageData.shape[2]
		if numCols> maxCols:
			maxCols = numCols
		
	
	nMatrix = 5 # choose number of columns
	mMatrix = math.ceil(numFiles / nMatrix) # make sure number of rows fit numFiles
	
	print('  ', 'num files:', numFiles, mMatrix, 'X', nMatrix, 'maxSlices:', maxSlices, 'maxRows:', maxRows, 'maxCols:', maxCols)

	data = np.ndarray((maxSlices, maxRows*mMatrix, maxCols*nMatrix), dtype=np.uint8)
	fullMask = np.ndarray((maxSlices, maxRows*mMatrix, maxCols*nMatrix), dtype=np.uint8)
	erodedMask = np.ndarray((maxSlices, maxRows*mMatrix, maxCols*nMatrix), dtype=np.uint8)
	ringMask = np.ndarray((maxSlices, maxRows*mMatrix, maxCols*nMatrix), dtype=np.uint8)
	#
	dvMask = np.ndarray((maxSlices, maxRows*mMatrix, maxCols*nMatrix), dtype=np.uint8)
	dvSkel = np.ndarray((maxSlices, maxRows*mMatrix, maxCols*nMatrix), dtype=np.uint8)
	rawCh1 = np.ndarray((maxSlices, maxRows*mMatrix, maxCols*nMatrix), dtype=np.uint8)
	
	print('  final data.shape:', data.shape)
	
	startRow = 0
	currentIdx = 0
	for i in np.arange(mMatrix):
		startCol = 0
		for j in np.arange(nMatrix):
			if currentIdx > len(dataList)-1:
				continue
			currentStack = dataList[currentIdx]
			
			startSlice = 0 # always
			stopSlice = startSlice + currentStack.shape[0]
			stopRow = startRow + currentStack.shape[1]
			stopCol = startCol + currentStack.shape[2]

			#print('i:', i, 'j:', j, 'slices:', startSlice, stopSlice, 'rows:', startRow, stopRow, 'cols:', startCol, stopCol, 'currentStack:', currentStack.shape)
			data[startSlice:stopSlice, startRow:stopRow, startCol:stopCol] = currentStack
			erodedMask[startSlice:stopSlice, startRow:stopRow, startCol:stopCol] = erodedList[currentIdx]
			ringMask[startSlice:stopSlice, startRow:stopRow, startCol:stopCol] = ringList[currentIdx]
			#
			fullMask[startSlice:stopSlice, startRow:stopRow, startCol:stopCol] = fullMaskList[currentIdx]
			dvMask[startSlice:stopSlice, startRow:stopRow, startCol:stopCol] = dvMaskList[currentIdx]
			dvSkel[startSlice:stopSlice, startRow:stopRow, startCol:stopCol] = dvSkelList[currentIdx]
			rawCh1[startSlice:stopSlice, startRow:stopRow, startCol:stopCol] = rawCh1List[currentIdx]
			
			# increment
			startCol += maxCols
			currentIdx += 1

		#
		startRow += maxRows
			
			
	with napari.gui_qt():
		title = gAnalysisPath + ":" + batchFilePath + '(n=' + str(numFiles) + ')'
		viewer = napari.Viewer(title=title)
		
		viewer.add_image(data=fullMask, name='full mask', contrast_limits=[0,1], opacity=0.6, colormap='gray')

		viewer.add_image(data=dvMask, name='dv mask', contrast_limits=[0,1], opacity=0.6, colormap='yellow', visible=False)
		viewer.add_image(data=dvSkel, name='dv skel', contrast_limits=[0,1], opacity=0.6, colormap='cyan', visible=False)
		
		viewer.add_image(data=erodedMask, name='eroded mask', contrast_limits=[0,1], opacity=0.6, colormap='inferno', visible=False)
		viewer.add_image(data=ringMask, name='ring mask', contrast_limits=[0,1], opacity=0.6, colormap='magenta', visible=False)

		viewer.add_image(data=data, name='raw ch2', contrast_limits=[0,110], opacity=0.6, colormap='green')
		viewer.add_image(data=rawCh1, name='raw ch1', contrast_limits=[0,110], opacity=0.6, colormap='red')

def samiNapari2(path):
	"""
	path: full patch to _ch2 stack
	"""
		
	fileNameNoExtension, _ = os.path.splitext(path)
	
	with napari.gui_qt():
		
		erodedPath = fileNameNoExtension + '_erodedMask.tif'
		ringPath = fileNameNoExtension + '_ringMask.tif'
		
		imageData = tifffile.imread(path)
		erodedData = tifffile.imread(erodedPath)
		ringData = tifffile.imread(ringPath)
		
		viewer = napari.Viewer(title=fileNameNoExtension, ndisplay=3)
		
		viewer.add_image(data=imageData, name='raw', colormap='green')
		viewer.add_image(data=erodedData, name='eroded', contrast_limits=[0,1], colormap='blue')
		viewer.add_image(data=ringData, name='ring', contrast_limits=[0,1], colormap='cyan')

if __name__ == '__main__':
	
	from gAnalysisPath import gSami_Params
	gAnalysisPath = gSami_Params['gAnalysisPath']
	
	gAnalysisPath = '/Users/cudmore/Desktop/samiVolume3'
	
	batchFilePath = ''
	for i, arg in enumerate(sys.argv):
		if i == 0:
			# my .py filename
			continue
		#print('i:', i, 'arg:', arg)
		if i==1:
			batchFilePath = arg

	#path = '/Users/cudmore/Desktop/samiVolume/200108/WT_Female/Cell_1/1_5ADVMLEG1L1_ch2.tif'
	#samiNapari2(path)
	
	if not batchFilePath:
		batchFilePath = '../analysis/wt-female.txt'
		#batchFilePath = '../analysis/ko-female.txt'
	
		#batchFilePath = '../analysis/wt-male.txt'
		#batchFilePath = '../analysis/ko-male.txt'
	
	samiDaskNapari(gAnalysisPath, batchFilePath)