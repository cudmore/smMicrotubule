"""
Author: RObert Cudmore
Date: 20200424

	Fiji API reference:
		ImagePlus: https://imagej.nih.gov/ij/developer/api/ij/ImagePlus.html
		ChannelSplitter: https://imagej.nih.gov/ij/developer/api/ij/plugin/ChannelSplitter.html
"""

from __future__ import print_function  # Only needed for Python 2, MUST be first import
import os, sys
import json

from ij import IJ
from ij.io import DirectoryChooser, OpenDialog
from ij.plugin import ChannelSplitter
from ij.process import ImageConverter

# lifeline version is 1.51n99, newer (April 2020) version is 1.52p99
fijiVersion = IJ.getFullVersion()
print('ImageJ/Fiji Version:', fijiVersion)

# allows us to import our own .py files, e.g. 'import simpleFileInfo
import inspect
print(' . ', inspect.getsourcefile(lambda:0))
myFilePath = os.path.dirname(os.path.abspath(inspect.getsourcefile(lambda:0)))
print('    myFilePath:', myFilePath)
sys.path.append(myFilePath)

# must be in same folder as this file 
# and folder needs an empty __init__.py file
import simpleFileInfo 

def convertAndSaveFile(fullFilePath, convertTo8Bit=False):
	"""
	"""
	print('  convertAndSaveFile() fullFilePath:', fullFilePath)
	
	folderPath, fileName = os.path.split(fullFilePath)
	fileNameNoExtension = os.path.splitext(fileName)[0]
	saveFileNoExtension = os.path.join(folderPath, fileNameNoExtension)
	
	#
	# load file and build a dict of parameters like (channels, pixels, voxels)
	myFileInfo, imp = simpleFileInfo.myLoadFile(fullFilePath, doShow=True)
	
	if convertTo8Bit:
		# from (ImagePlus.GRAY8, ImagePlus.GRAY16, ImagePlus.GRAY32, ImagePlus.COLOR_256 or ImagePlus.COLOR_RGB)
		myType = imp.getType()
		if myType in [1,2]:
			ic = ImageConverter(imp) # converts channelImp in place
			scaleConversions = True
			ic.setDoScaling(scaleConversions) # scales inensities to fill 0...255
			ic.convertToGray8()		

	#
	# split channels
	channelArray = ChannelSplitter.split(imp) # returns an array of imp
	
	for channelIdx, channelImp in enumerate(channelArray):
		# channelImp is an imp, this will NOT have fileinfo (we just created it)
		#channelImp.show()
		
		saveFilePath = saveFileNoExtension + '_ch' + str(channelIdx+1) + '.tif'
		print('    ch:', channelIdx+1, 'saveFilePath:', saveFilePath)
		
		IJ.save(channelImp, saveFilePath)

	return myFileInfo
	
def convertFolder(folderPath, theseFileExtensions=['.oir'], thisFileEnding=None, convertTo8Bit=False):
	"""
	Recursively traverse down folderPath
	and open/process all files ending in one of [theseFileExtensions]

	save a _db.csv file in folderPath with image acquisition parameters 
	"""
	
	if not os.path.isdir(folderPath):
		print('ERROR: convertFolder() did not get a folder path folderPath:', folderPath)
		return
	
	print('convertFolder() folderPath:', folderPath)
	
	myInfoDictList = []
	
	for dirpath, dirnames, files in os.walk(folderPath):
		for name in files:
			fileNameNoExtension, fileExtension = os.path.splitext(name)
			if fileExtension in theseFileExtensions:
				if thisFileEnding is not None:
					if name.endswith(thisFileEnding):
						pass
					else:
						continue
				filePath = os.path.join(dirpath, name)
				#
				myInfoDict = convertAndSaveFile(filePath, convertTo8Bit=convertTo8Bit)
				#
				myInfoDictList.append(myInfoDict)
	
	# number of files we processed
	numFiles = len(myInfoDictList)
	
	#
	# save
	enclosingFolderName = os.path.split(folderPath)[1]
	saveDatabaseFilePath = os.path.join(folderPath, enclosingFolderName + '_db.csv')
	
	print('    processed', numFiles, 'files, saving database in', saveDatabaseFilePath)
	
	with open(saveDatabaseFilePath, 'w') as f:
		for rowIdx, myInfoDict in enumerate(myInfoDictList):
			if rowIdx==0:
				# headers
				headerStr = ''
				for k,v in myInfoDict.items():
					headerStr += k + ','
				print(headerStr, file=f)
			rowStr = ''
			for k,v in myInfoDict.items():
				rowStr += str(v) + ','
			print(rowStr, file=f)

if __name__ in ['__main__', '__builtin__']: 

	###
	###
	### PARAMETERS
	###
	###
	thisFileEnding = None # to process all files with extension in theseFileExtension
	#thisFileEnding = 'ADVMLEG1L1.oir' # deconvolved is like 12_5ADVMLEG1L1.oir
	#theseFileExtensions=['.tif']
	theseFileExtensions=['.oir']
	convertTo8Bit = True

	print(' . thisFileEnding:', thisFileEnding)
	print(' . theseFileExtensions:', theseFileExtensions)
	print(' . convertTo8Bit:', convertTo8Bit)
	
	#
	# ask user for a folder path
	userFolder = DirectoryChooser('Choose a folder')
	folderPath = userFolder.getDirectory()
	if folderPath is None:
		print('Warning: User did not choose a folder')	
	else:
		# remove trailing /  or \
		if folderPath.endswith(os.sep):
			folderPath = folderPath[:-1]		
		#
		convertFolder(folderPath, theseFileExtensions=theseFileExtensions, thisFileEnding=thisFileEnding, convertTo8Bit=convertTo8Bit)
		#

	print('simpleConvert is done !!!')
