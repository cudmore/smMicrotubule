"""
Given a full path to original .tif file, open viewer with 3 volumes
	1) original
	2) _dvMask
	3) _dvSkel

On macOS, you can get the full path to a file by:
	1) drag and drop the file into a command prompt
	2) right-click the file, press option key, and select menu 'Copy ... As Pathname'

Usage:

	# to view analysis output 1 (one)
	python samiViewer /Users/cudmore/box/data/sami/Cell_1/1_5ADVMLEG1L1_ch2.tif

	python samiViewer.py /Users/cudmore/box/data/sami/191230/BIN1_smKO_Male/Cell_12/12_5ADVMLEG1L1_ch2.tif 
	
	# to view analysis output 2, you need to specify '2'
	python samiViewer /Users/cudmore/box/data/sami/Cell_1/1_5ADVMLEG1L1_ch2.tif 2

"""

import os, sys

#from skimage import data
import napari

#path = '/Volumes/fourt/2020/resize/cell1.tif'

if __name__ == '__main__':

	# get (path, analysis number) from command line
	path = None
	saveNumber = None
	for i, arg in enumerate(sys.argv):
		if i == 0:
			# my .py filename
			continue
		#print('i:', i, 'arg:', arg)
		if i==1:
			path = arg
		elif i == 2:
			saveNumber = arg

	if saveNumber is None:
		saveNumber = 1
		
	with napari.gui_qt():
		#viewer = napari.view_image(data.astronaut(), rgb=True)
		
		if path is not None:
			# open viewer with (path, _dvMask, dvSkel)
			pathNoExtension = os.path.splitext(path)[0]
			maskPath = pathNoExtension + '_dvMask_' + str(saveNumber) + '.tif'
			skelPath = pathNoExtension + '_dvSkel_' + str(saveNumber) + '.tif'
			
			title = pathNoExtension + ' ' + str(saveNumber)
			viewer = napari.Viewer(title=title, ndisplay=3)
			viewer.add_image(path=path, name='raw', colormap='green')
			viewer.add_image(path=maskPath, name='dvMask', opacity=0.6, visible=False, colormap='gray')
			viewer.add_image(path=skelPath, name='dvSkel', opacity=0.6, colormap='red')
			
		else:
			# open empty viewer
			viewer = napari.Viewer()

