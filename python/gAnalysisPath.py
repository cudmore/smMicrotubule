"""
Script with parameters to control Sami analysis

This specifies the save path and parameters for analysis to be used in:
	samiAnalysis.py
	samiVolume2.py
	samiDensity.py (no parameters)

To use this, specify a block like 'case 1', make sure
	- it has a unique output folder gSami_Params['gAnalysisPath']
	- specify detection parameters

To troubleshoot and run on just a few files, set batchFileList
	batchFileList = ['/Users/cudmore/Sites/smMicrotubule/analysis/two_wt-female.txt']

"""

#gAnalysisPath = '/Users/cudmore/Desktop/samiVolume3' # pre erode 3
#gAnalysisPath = '/Users/cudmore/Desktop/samiVolume4' # pre erode 1
#gAnalysisPath = '/Users/cudmore/Desktop/samiVolume5' # pre erode 2

from collections import OrderedDict

gSami_Params = OrderedDict()
gSami_Params['gAnalysisPath'] = ''
gSami_Params['gBatchFileList'] = []
gSami_Params['samiAnalysis_params'] = OrderedDict()
gSami_Params['samiVolume2_params'] = OrderedDict()
gSami_Params['samiDensity_params'] = OrderedDict()

#
# todo: make samiAnalysis.py use batch file list !!!
#

#
batchFileList = []

# normally all 4 (genotype, sex)
batchFile = '/Users/cudmore/Sites/smMicrotubule/analysis/wt-female.txt'
batchFileList.append(batchFile)
batchFile = '/Users/cudmore/Sites/smMicrotubule/analysis/wt-male.txt'
batchFileList.append(batchFile)
batchFile = '/Users/cudmore/Sites/smMicrotubule/analysis/ko-female.txt'
batchFileList.append(batchFile)
batchFile = '/Users/cudmore/Sites/smMicrotubule/analysis/ko-male.txt'
batchFileList.append(batchFile)

# or debugging
#batchFileList = ['/Users/cudmore/Sites/smMicrotubule/analysis/two_wt-female.txt']

#
gSami_Params['gBatchFileList'] = batchFileList

#
# specify different 'gAnalysisPath' folder with different parameters
#

#
# case 1
gSami_Params['gAnalysisPath'] = '/Users/cudmore/Desktop/samiVolume1'
gSami_Params['samiAnalysis_params']['f3_param'] = [0.5, 0.005]
gSami_Params['samiAnalysis_params']['minArea'] = 20 # pixels

gSami_Params['samiVolume2_params']['medianParams'] = (3,5,5)	# for mediam filter
gSami_Params['samiVolume2_params']['thresholdMin'] = 2			# for myThreshold_min_max()
gSami_Params['samiVolume2_params']['thresholdMax'] = 255
gSami_Params['samiVolume2_params']['includeDistanceLessThan'] = 2 # distance map, includeDistanceLessThan = 2 * xVoxel
gSami_Params['samiVolume2_params']['erodedIterations0'] = 2 # after filling holes, filledHolesMask by 1 iteration 
gSami_Params['samiVolume2_params']['erodedIterations'] = 4 # erode to get final eroded

# gSami_Params['samiDensity_params'] # No parameters !!!!

#
# todo: before i do more cases, i need to remove nucleus from masks (eroded, ring)
#

# case 2
gSami_Params['gAnalysisPath'] = '/Users/cudmore/Desktop/samiVolume2'
gSami_Params['samiAnalysis_params']['f3_param'] = [1, 0.01]
gSami_Params['samiAnalysis_params']['minArea'] = 20 # pixels

gSami_Params['samiVolume2_params']['medianParams'] = (3,5,5)	# for mediam filter
gSami_Params['samiVolume2_params']['thresholdMin'] = 2			# for myThreshold_min_max()
gSami_Params['samiVolume2_params']['thresholdMax'] = 255
gSami_Params['samiVolume2_params']['includeDistanceLessThan'] = 2 # distance map, includeDistanceLessThan = 2 * xVoxel
gSami_Params['samiVolume2_params']['erodedIterations0'] = 2 # after filling holes, filledHolesMask by 1 iteration 
gSami_Params['samiVolume2_params']['erodedIterations'] = 4 # erode to get final eroded


#
#
#
# case 1
gSami_Params['gAnalysisPath'] = '/Users/cudmore/Desktop/samiVolume3'
gSami_Params['samiAnalysis_params']['f3_param'] = [0.5, 0.005]
gSami_Params['samiAnalysis_params']['minArea'] = 20 # pixels

gSami_Params['samiVolume2_params']['medianParams'] = (3,5,5)	# for mediam filter
gSami_Params['samiVolume2_params']['thresholdMin'] = 2			# for myThreshold_min_max()
gSami_Params['samiVolume2_params']['thresholdMax'] = 255
gSami_Params['samiVolume2_params']['includeDistanceLessThan'] = 2 # distance map, includeDistanceLessThan = 2 * xVoxel
gSami_Params['samiVolume2_params']['erodedIterations0'] = 3 # after filling holes, filledHolesMask by 1 iteration 
gSami_Params['samiVolume2_params']['erodedIterations'] = 5 # erode to get final eroded
