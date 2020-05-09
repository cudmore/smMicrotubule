"""
	Author: Robert Cudmore
	Date: 20200411

	Purpose:
		Take a raw image stack and produce a segmentation to count the number of pixels in the largest segmented object

	Input:
		8-bit image
	
	Output:
		comma-seperated list of detection parameters and largest 3 segmented objects
	
	Algorithm:
		1) Gaussian blur
		2) Threshold
		3) Fill holes
		4) label connected component 1, 2, 3, ...
		5) calculate volume of largest component (v = (number of pixels * (width*height*depth))
		
	Requirements:
		Using MorphoLibJ, see https://imagej.net/MorphoLibJ
		
	To install MorphoLibJ in Fiji/ImageJ 2
		See: https://imagej.net/MorphoLibJ#Installation
		1) In ImageJ2 (including Fiji), you just need to add the IJPB-plugins site to your list of update sites:
		2) Select Help    Update... from the menu to start the updater.
		3) Click on Manage update sites. This brings up a dialog where you can activate additional update sites.
		4) Activate the IJPB-plugins update site and close the dialog. Now you should see an additional jar file for download.
		5) Click Apply changes and restart ImageJ.

	Running the script:
		(1) Run from within Fiji
  			Drag and drop this .py file into Fiji (will open in Fiji text editor)
  			Open a .tif to process
  			Go back into the Fiji text editor and press keyboard command+r (or select run-run menu)
  	
		(2) Run from the macOS terminal/command-line
			# assuming you have Fiji installed in your /applications folder
			/Applications/Fiji.app/Contents/MacOS/ImageJ-macosx fiji/samiVolume.py batch=/Users/cudmore/box/data/sami/analysis/wt-male.txt

			Where 'wt-male.txt' is a text file with a list of full .tif file paths
  	
	Notes:
		See:
			https://imagej.nih.gov/ij/developer/api/ij/ImagePlus.html
		
		Always look at these for hints:
			https://imagej.net/Jython_Scripting
			https://ilovesymposia.com/2014/02/26/fiji-jython/

		Post that hinted how to do this with MorphoLibJ
			https://forum.image.sc/t/selecting-connected-pixels-in-a-binary-mask-in-3d/4142/2

"""

from __future__ import print_function

import os, sys
from collections import OrderedDict

from ij import IJ, WindowManager
from ij import ImagePlus # to query ImagePlus.GRAY8
from ij.measure import ResultsTable

# see: https://javadoc.scijava.org/MorphoLibJ/index.html?inra/ijpb/label/LabelImages.html
from inra.ijpb.label import LabelImages

#####################################################################
def myRun(path=None, imp=None):

	theRet = OrderedDict()
	theRet['path'] = ''
	theRet['file'] = ''
	theRet['cond'] = '' # calling function needs to fill this in
	theRet['xVoxel'] = None
	theRet['yVoxel'] = None
	theRet['zVoxel'] = None
	
	if path is not None:
		# load from a file path
		print(' .  myRun() path:', path)
		theRet['path'] = path
		theRet['file'] = os.path.split(path)[1]
		pathNoExtension = os.path.splitext(path)[0]
		imp = IJ.openImage(path)
		imp.show()
	elif imp is not None:
		# use open image in Fiji
		imageTitle = imp.getTitle() # title of the window, usually the filename
		print(' .  myRun() imageTitle:', imageTitle)
		myFileInfo = imp.getOriginalFileInfo() # can be None when file was not dragged/dropped
		if myFileInfo is not None:
			#print('myFileFinfo:', myFileInfo.directory, myFileInfo.fileName)
			tmpPath = os.path.join(myFileInfo.directory, myFileInfo.fileName)
			pathNoExtension = os.path.splitext(tmpPath)[0]
			theRet['path'] = tmpPath
			theRet['file'] = myFileInfo.fileName
		else:
			# happens when .tif was not dragged and dropped by the user, there is no file or path !!!
			theRet['file'] = imageTitle
	else:
		print('samiVolume.py does not know what to do: path:', path, 'imp:', imp)
		return None
		
	imageType = imp.getType()
	if imageType != ImagePlus.GRAY8:
		# abort
		print('error: samiVolume() expects an 8-bit image. path:', path)
		return None
		
	
	cal = imp.getCalibration() # cal is type 'ij.measure.Calibration'
	#print(' .   width:', cal.pixelWidth, 'height:', cal.pixelHeight, 'depth:', cal.pixelDepth)
	theRet['xVoxel'] = cal.pixelWidth
	theRet['yVoxel'] = cal.pixelHeight
	theRet['zVoxel'] = cal.pixelDepth

	# blur
	xGaussian = 2
	yGaussian = 2
	zGaussian = 2
	theRet['xGaussian'] = xGaussian
	theRet['yGaussian'] = yGaussian
	theRet['zGaussian'] = zGaussian
	paramStr = "x=" + str(xGaussian) + ' y=' + str(yGaussian) + ' z=' + str(zGaussian)
	#IJ.run("Gaussian Blur 3D...", "x=2 y=2 z=2");
	IJ.run("Gaussian Blur 3D...", paramStr);

	#
	# threshold
	IJ.setAutoThreshold(imp, "Default dark");
	# run("Threshold...");
	minThreshold = 2
	maxThreshold = 255
	#IJ.setThreshold(2, 255);
	theRet['minThreshold'] = minThreshold
	theRet['maxThresholdeshold'] = maxThreshold
	IJ.setThreshold(minThreshold, maxThreshold)
	IJ.run("Convert to Mask", "method=Default background=Dark list")

	#
	# fill holes
	# this make a new window -fillHoles, 8-bit
	IJ.run("Fill Holes (Binary/Gray)");

	#
	# label 3d connected components in binary mask, each connected component has intensity 1,2,3
	# this make a new window -lbl, 16-bit
	IJ.run("Connected Components Labeling", "connectivity=6 type=[16 bits]");

	#imp_lbl = WindowManager.getCurrentImage()
	imp_lbl = IJ.getImage()

	#
	# take the image with the connected component labels and step through each label and get pixel counts
	# Returns the set of unique labels existing in the given image, excluding the value zero (used for background).
	labels = LabelImages.findAllLabels( imp_lbl )

	pixelCountArray = LabelImages.pixelCount(imp_lbl, labels)
	#print('pixelCountArray:', type(pixelCountArray), pixelCountArray)

	sortedArray = sorted(pixelCountArray, reverse=True)
	#for idx, pixelCount in enumerate(pixelCountArray):
	theRet['pixelCount_1'] = None
	theRet['pixelCount_2'] = None
	theRet['pixelCount_3'] = None
	for idx, pixelCount in enumerate(sortedArray):
		labelIdx = idx + 1 # because label 0 (background is skipped)
		theRet['pixelCount_' + str(labelIdx)] = pixelCount
		#print(' .   label:', labelIdx, 'pixelCount:', pixelCount)
		if labelIdx >= 3:
			break
	
	#
	# todo: we need to make 2 copies of the mask, one that is reduced (intracellular) and the remaining shell (membrane)
	
	#
	# save the mask
	maskFilePath = os.path.join(pathNoExtension, '_dvVolMask_Full.tif')
	#IJ.save
	
	return theRet

#####################################################################
def parseBatchFile(path):
	"""
	given a .txt file with one path per line, return a list paths
	"""
	with open(path) as f:
		content = f.readlines()
	# remove whitespace characters like `\n` at the end of each line
	content = [x.strip() for x in content] 
	return content

#####################################################################
def closeAll():
	""" Close all open windows """
	if WindowManager.getIDList() is None:
		return
	for id in WindowManager.getIDList():
		imp = WindowManager.getImage(id)
		imp.close()
		
#####################################################################
if __name__ in ['__builtin__','__main__']:
	"""
	__main__ when called from command line
	__builtin__ when called from inside fiji
	"""
	
	commandLine = True
	if __name__ == '__builtin__':
		commandLine = False
	
	#print('__name__ =', __name__)
	myName = '' # the name of this file
	for i, arg in enumerate(sys.argv):
		#print(' .   arg:', i, arg)
		if i == 0:
			# my .py filename (not full path)
			myName = arg
			continue
		elif arg.startswith('batch='):
			tmp, batchFilePath = arg.split('=')
			print('processing batch file:', batchFilePath)
			pathList = parseBatchFile(batchFilePath)
			savePath = os.path.splitext(batchFilePath)[0]
			print(' .   savePath:', savePath)
			batchFile = os.path.basename(batchFilePath)
			print(' .   batchFile:', batchFile)
			batchFileNoExtension = os.path.splitext(os.path.basename(batchFilePath))[0]
			print(' .   batchFileNoExtension:', batchFileNoExtension)
			'''
			for path in pathList:
				print(' .   path:', path)
			print(' .   savePath:', savePath)
			print(' .   batchFileNoExtension:', batchFileNoExtension)
			'''
			
	resultList = []
	if commandLine:
		nPath = len(pathList)
		for idx, path in enumerate(pathList):
			print('[', idx+1, 'of', nPath, ']')
			result = myRun(path=path)
			if result is not None:
				result['cond'] = batchFile
				resultList.append(result)
		closeAll()
		
	else:
		# assume we have a window with an imp
		imp = IJ.getImage()
		result = myRun(imp=imp)
		if result is not None:
			resultList.append(result)
		
	#
	# print results
	# This is a mixture of header/results str in python/jython and a ImageJ ResultsTable() ... confusing
	headerStr = '' # header using first result (assuming all results have same keys)
	resultStr = ''
	nResult = len(resultList)
	if nResult > 0:
		myResultsTable = ResultsTable()
	for idx, result in enumerate(resultList):
		myResultsTable.incrementCounter() # add a row to results table
		for k,v in result.items():
			if idx==0:
				headerStr += k + ','
			resultStr += str(v) + ','
			if v is None:
				# need this because on save because Python None get converted to null and generates and error
				v = 'None'
			myResultsTable.addValue(k, v)
		if idx==0:
			headerStr += '\n'
		resultStr += '\n'

	if nResult > 0:
		myResultsTable.setPrecision(6)
		myResultsTable.show('samiVolume Results')
		#myResultsTable.save('/Users/cudmore/box/data/sami/volumeResults.csv')
		myResultsTable.save(savePath + '-volume.csv')
		#myResultsTable.close()
		IJ.run("Close");
		
	if nResult > 0:
		print('=== results headers followed by values ===')
		print(headerStr)
		print(resultStr)
	else:
		print('No Results ???')
	#
	# put results in a table (only useful for running from Fiji with single image
	# done
	print('Finished', myName)