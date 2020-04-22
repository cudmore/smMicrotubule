"""
	Robert Cudmore
	20200412
	
	See:
		https://imagej.net/AnalyzeSkeleton
		https://javadoc.scijava.org/Fiji/sc/fiji/analyzeSkeleton/package-summary.html
		https://javadoc.scijava.org/Fiji/sc/fiji/analyzeSkeleton/SkeletonResult.html
		
"""

import os, sys, math
from collections import OrderedDict

#from skeleton_analysis import AnalyzeSkeleton_,Graph,Edge,Vertex

from ij import IJ, WindowManager
from ij import ImagePlus # to query ImagePlus.GRAY8
from ij.measure import ResultsTable

import sc.fiji.analyzeSkeleton #.AnalyzeSkeleton_

#####################################################################
def myRun(path=None, imp=None):
	"""
	"""
	theRet = OrderedDict()
	
	if path is not None:
		# convert full path to raw data to _dvSkel.tif
		savePath, fileName = os.path.split(path)
		pathFileNoExtension, ext= os.path.splitext(fileName)
		dvSkelPath = os.path.join(savePath, pathFileNoExtension+'_dvSkel.tif')
		#
		#imp = IJ.openImage(path)
		print('dvSkelPath() loading:', dvSkelPath)
		imp = IJ.openImage(dvSkelPath)
		imp.show()
		
	cal = imp.getCalibration() # cal is type 'ij.measure.Calibration'
	#print(' .   width:', cal.pixelWidth, 'height:', cal.pixelHeight, 'depth:', cal.pixelDepth)
	pixelWidth = cal.pixelWidth
	pixelHeight = cal.pixelHeight
	pixelDepth = cal.pixelDepth

	'''
	skelImp = imp.duplicate() # makes DUP_ window
	skelImp.show()
	'''
	skelImp = imp
	
	doOnRaw = False
	if doOnRaw:
		# median filter
		#IJ.run(imp,"Gaussian Blur...","sigma=5")
		IJ.run("Median 3D...", "x=2 y=2 z=2")
		
		# threshold and make binary
		#IJ.setAutoThreshold("Otsu dark")
		IJ.run("Convert to Mask", "method=Otsu background=Dark calculate");

	
	# skeletonize
	IJ.run("Skeletonize (2D/3D)");
	 
	skelImp.changes= False # will not prompt to save when closing imp window

	# analyze skeleton
	# makes 2x tables 'Results' and 'Branch Information'
	# makes 2x window 'Tagged skeleton' and <orig image>-labeled-skeletons
	'''
	paramStr = 'prune=[shortest branch] show display'
	IJ.run("Analyze Skeleton (2D/3D)", paramStr)
	'''
	skel = sc.fiji.analyzeSkeleton.AnalyzeSkeleton_()
	skel.setup("", skelImp)
	
	# .NONE : no pruning mode index
	# .LOWEST_INTENSITY_BRANCH : lowest intensity branch
	pruneIdx = sc.fiji.analyzeSkeleton.AnalyzeSkeleton_.NONE # The pruneIndex, as asked by the initial gui dialog.
	pruneEnds = False # flag to prune end-point-ending branches
	shortPath = False # flag to calculate the longest shortest path
	silent = False # 
	verbose = True
	
	# run(int pruneIndex, boolean pruneEnds, boolean shortPath, ij.ImagePlus origIP, boolean silent, boolean verbose)
	skelResult = skel.run(pruneIdx, pruneEnds, shortPath, skelImp, silent, verbose);
	print('skelResult:', type(skelResult), skelResult)
	
	'''
	# run AnalyzeSkeleton
	# (see https://fiji.sc/AnalyzeSkeleton 
	# and https://fiji.sc/javadoc/skeleton_analysis/package-summary.html)
	skel = AnalyzeSkeleton_()
	skel.setup("",imp)
	skelResult = skel.run(skel.NONE, False, True, None, True, True)
 	'''
 	
	# get the separate skeletons
	graph = skelResult.getGraph()
	print len(graph)
	print skelResult.getNumOfTrees()
 
 	myResults(skelResult)
 	
	# find the longest graph
	graph = sorted(graph, key=lambda g: getGraphLength(g), reverse=True)
	longestGraph = graph[0]
 	print('longestGraph:', longestGraph)
 	
	# find the longest edge
	edges = longestGraph.getEdges()
	edges = sorted(edges, key=lambda edge: edge.getLength(), reverse=True)
	longestEdge = edges[0]
	print('longestEdge:', longestEdge)
	
	
#####################################################################
def getGraphLength(graph):
	length = 0
	for g in graph.getEdges():
		length = length + g.getLength()
	return length

#####################################################################
def myResults(results):
	myResultsTable = ResultsTable()
	for idx, graph in enumerate(results.getGraph()):
		for edge in graph.getEdges():
			edgeLength = edge.getLength()
			v1 = edge.getV1()
			v2 = edge.getV2()
			dist = euclideanDistance(v1, v2)
			#print('v1:', type(v1), v1.getPoints())
			#
			myResultsTable.incrementCounter() # add a row to results table
			myResultsTable.addValue('graphID', idx)
			myResultsTable.addValue('length_3d', edgeLength)
			myResultsTable.addValue('dist', dist)
			if dist > 0:
				myResultsTable.addValue('tort', edgeLength / dist)
			else:
				myResultsTable.addValue('tort', 'inf')

	myResultsTable.setPrecision(6)
	myResultsTable.show('samiSkel_results')
	#myResultsTable.save('/Users/cudmore/box/data/sami/volumeResults.csv')

#####################################################################
def euclideanDistance(v1, v2):
	"""
	v1 and v2: Vertex
	"""
	x1 = v1.getPoints()[0].x
	y1 = v1.getPoints()[0].y
	z1 = v1.getPoints()[0].z

	x2 = v2.getPoints()[0].x
	y2 = v2.getPoints()[0].y
	z2 = v2.getPoints()[0].z

	dist = math.sqrt((x2-x1)**2 + (y2-y1)**2 + (z2-z1)**2)

	return dist
	
#####################################################################
def closeAll():
	for imp in map(WindowManager.getImage, WindowManager.getIDList()):
		imp.close()
	
#####################################################################
if __name__ in ['__builtin__','__main__']:

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
		path = '/Users/cudmore/box/data/sami/191230/WT_Male/Cell_4_/4_1_5ADVMLEG1L1_ch2.tif'
		'''
		imp = IJ.openImage(path)
		imp.show()
		'''
		
		# debug, remove and use user specified front window
		result = myRun(path=path)

		# assume we have a window with an imp
		'''
		imp = IJ.getImage()
		result = myRun(imp=imp)
		'''

		if result is not None:
			resultList.append(result)
