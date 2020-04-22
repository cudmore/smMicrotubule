"""
delete all _dvMask and _dvSkel
We need this because save and replace of .tif is tricky (for good reason)
"""

import os

def deleteMaskAndSkel(path):
	for dirpath, dirnames, files in os.walk(path):
		for name in files:
			if name.endswith('_dvMask.tif'):
				filePath = os.path.join(dirpath, name)
				print('removing:', filePath)
				os.remove(filePath)
			if name.endswith('_dvSkel.tif'):
				filePath = os.path.join(dirpath, name)
				print('removing:', filePath)
				os.remove(filePath)
				
path = '/Users/cudmore/box/data/sami/191230'
deleteMaskAndSkel(path)

path = '/Users/cudmore/box/data/sami/200108'
deleteMaskAndSkel(path)