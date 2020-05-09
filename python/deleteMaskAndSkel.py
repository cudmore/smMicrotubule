"""
delete all _dvMask and _dvSkel
We need this because save and replace of .tif is tricky (for good reason)
"""

import os

def deleteMaskAndSkel(path):
	if not os.path.isdir(path):
		print('ERROR: deleteMaskAndSkel() did not find path:', path)
	for dirpath, dirnames, files in os.walk(path):
		for name in files:
			if name.endswith('_dvMask_1.tif'):
				filePath = os.path.join(dirpath, name)
				print('    removing:', filePath)
				os.remove(filePath)
			if name.endswith('_dvSkel_1.tif'):
				filePath = os.path.join(dirpath, name)
				print('    removing:', filePath)
				os.remove(filePath)
		
			if name.endswith('_dvMask.tif'):
				filePath = os.path.join(dirpath, name)
				print('    removing:', filePath)
				os.remove(filePath)
			if name.endswith('_dvSkel.tif'):
				filePath = os.path.join(dirpath, name)
				print('    removing:', filePath)
				os.remove(filePath)
		
			if name.endswith('_ch1.tif'):
				filePath = os.path.join(dirpath, name)
				print('    removing:', filePath)
				os.remove(filePath)
			if name.endswith('_ch2.tif'):
				filePath = os.path.join(dirpath, name)
				print('    removing:', filePath)
				os.remove(filePath)
				
path = '../data/191230'
deleteMaskAndSkel(path)

path = '../data/200108'
deleteMaskAndSkel(path)

path = '../data/200414'
deleteMaskAndSkel(path)

path = '../data/200421'
deleteMaskAndSkel(path)

