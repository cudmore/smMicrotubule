
"""
This Fiji Jython script will convert a two channel .oir xxx file into 2 (two) channel .tif files

Run this inside Fiji by selecting menu:
	Plugins - Macro - Run...

Assumes:
	The file you select is a proper .oir and has two channels

Caveats:
	1) This assumes all paths including folder names and file names do NOT have spaces !!!
	2) if you have a finder window open, the new files will not immediately appear
	3) This will overwite _ch1/_ch2 .tif file if they already exist
	
Recipe:
=======
run("Bio-Formats Importer", "open=/Users/cudmore/box/data/sami/Cell_1/1_5ADVMLEG1L1.oir color_mode=Default rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT");
run("8-bit");
run("Split Channels");
selectWindow("C1-1_5ADVMLEG1L1.oir");
run("Save", "save=/Users/cudmore/box/data/sami/Cell_1/c12_ch1.tif");
selectWindow("C2-1_5ADVMLEG1L1.oir");
run("Save", "save=/Users/cudmore/box/data/sami/Cell_1/c12_ch2.tif");
"""

# has to be the first line
from __future__ import print_function

import os, sys

from ij import IJ, WindowManager
from ij.io import OpenDialog
from ij.process import ImageConverter # by default convert to 8-bit will scale, i need to turn ot off. See: https://ilovesymposia.com/2014/02/26/fiji-jython/
# see: https://imagej.nih.gov/ij/developer/api/ij/process/ImageConverter.html
# Set true to scale to 0-255 when converting short to byte or float to byte and to 0-65535 when converting float to short.

from loci.plugins import BF

#path = '/Users/cudmore/box/data/sami/Cell_1/1_5ADVMLEG1L1.oir'

def myConvert(path='', imp=None):
	if path == '' and imp is None:
		# ask user for file. I do not know how to handle when users hits cancel??? Script will just fail
		notUsedBy_oxs = 'Open a two-channel deconvoluted .oir file'
		path = OpenDialog(notUsedBy_oxs).getPath()

	if len(path)>0:
		print('    user selected path:', path)
		fileName = os.path.basename(path)

		# load
		print(' .   loading:', path)
		imps = BF.openImagePlus(path)
		for imp in imps:
		    imp.show()

	# convert to 8-bit
	ImageConverter.setDoScaling(True)
	IJ.run('8-bit')
	
	# split channels
	IJ.run("Split Channels");
	
	
	# save channel 1
	windowName = 'C1-' + fileName
	print(' .   ch1 windowName:', windowName)
	IJ.selectWindow(windowName)
	windowName_notSure = WindowManager.getImage(windowName)
	saveFilePath = os.path.splitext(path)[0] + '_ch1.tif' # this is getting garbled when we have spaces (' ') in path !!!!!!!!!
	# remove spaces
	#saveFilePath = saveFilePath.replace(' ', '_')
	print('    saving', saveFilePath)
	IJ.save(windowName_notSure, saveFilePath)
	#IJ.run("Save", 'save=' + saveFilePath)
	#
	#IJ.close() # close channel 1 ... this does not work

	windowName_notSure.close()
	
	# save channel 2
	windowName = 'C2-' + fileName
	print(' .   ch2 windowName:', windowName)
	IJ.selectWindow(windowName)
	windowName_notSure = WindowManager.getImage(windowName)
	saveFilePath = os.path.splitext(path)[0] + '_ch2.tif'
	# remove spaces
	#saveFilePath = saveFilePath.replace(' ', '_')
	print('    saving', saveFilePath)
	IJ.save(windowName_notSure, saveFilePath)
	#IJ.run("Save", 'save=' + saveFilePath)
	#IJ.close() # close channel 2 ... this does not work
	windowName_notSure.close()
	
	
	# close all the windows ... this does not work
	'''
	print('    closing windows')
	closeAll()
	'''
	
	'''
	imps = BF.openImagePlus(path)
	for imp in imps:
	    imp.changes= False
	    imp.close()
	'''
	
	print('    done')

#####################################################################
def closeAll():
	for imp in map(WindowManager.getImage, WindowManager.getIDList()):
		imp.close()

#####################################################################
def parseBatchFile(path):
	"""
	given a .txt file with one path per line, return a list paths
	"""
	with open(path) as f:
		content = f.readlines()
	print('content:', content)
	# remove whitespace characters like `\n` at the end of each line
	content = [x.strip() for x in content] 
	return content

#####################################################################
if __name__ in ['__builtin__','__main__']:

	#
	# Fiji 2017 does not have '__builtin__'
	
	commandLine = False # we will never call this from the command line
	if __name__ == '__builtin__':
		commandLine = False
	
	print('__name__ =', __name__)
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

		resultList = []
	if commandLine:
		nPath = len(pathList)
		for idx, path in enumerate(pathList):
			print('[', idx+1, 'of', nPath, ']')
			#
			result = myConvert(path=path)
			#
			if result is not None:
				result['cond'] = batchFile
				resultList.append(result)
		closeAll()
	else:
		# specify path here
		# .endswith('_3ADVMLEG1L1.oir')
		#folderPath = '/Users/cudmore/box/data/sami/data/200414'
		folderPath = '/Users/cudmore/box/data/sami/data/200421'
	
		print('folderPath:', folderPath)

		for dirpath, dirnames, files in os.walk(folderPath):
			for name in files:
				if name.endswith('_3ADVMLEG1L1.oir'):
					filePath = os.path.join(dirpath, name)
					#x,y,z = readVoxelSize(filePath, verbose=False)
					print('    ', filePath)
					myConvert(path=filePath, imp=None)
		
		'''
		batchFilePath = '/Users/cudmore/box/data/sami/batchConvert.txt'
		pathList = parseBatchFile(batchFilePath)
		for path in pathList:
			if os.path.isfile(path):
				myConvert(path=path, imp=None)
			else:
				print('ERROR: file not found:', path)
		'''