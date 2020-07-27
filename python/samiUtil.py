import numpy as np

def removeZSpreadFromMask(maskData, removeNumSlices=4, myName=''):
	"""
	remove top/bottom slices to fix z-spread
	
	removeNumSlices:
		0 to remove 1 slice on top/bottom
		1 to remove 2 slices on top/bottom ...
	myName: used to print name of stack (to keep track of which one)
	"""
	numSlices = maskData.shape[0]
	
	# from top to bottom
	firstNonZeroSlice = None
	for i in range(numSlices):
		theCount = np.count_nonzero(maskData[i,:,:] > 0) # count pixels >0
		if theCount > 0 and firstNonZeroSlice is None:
			firstNonZeroSlice = i
		if firstNonZeroSlice is not None:
			if i >= firstNonZeroSlice and i <= (firstNonZeroSlice+removeNumSlices):
				#print('removeZSpreadFromMask()', myName, 'top-bottom removing slice', i)
				maskData[i,:,:] = 0
			else:
				break
	
	# from bottom to top
	lastNonZeroSlice = None
	for i in reversed(range(numSlices)):
		theCount = np.count_nonzero(maskData[i,:,:] > 0) # count pixels >0
		if theCount > 0 and lastNonZeroSlice is None:
			lastNonZeroSlice = i
		if lastNonZeroSlice is not None:
			if i <= lastNonZeroSlice and i >= (lastNonZeroSlice-removeNumSlices):
				#print('removeZSpreadFromMask()', myName, 'bottom-top removing slice', i)
				maskData[i,:,:] = 0
			else:
				break
	
	return maskData
	
