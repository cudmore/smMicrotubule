import os, sys

#folderPath = '/Users/cudmore/box/data/sami/data/200414'
folderPath = '/Users/cudmore/box/data/sami/data/200421'
pathList = []
for dirpath, dirnames, files in os.walk(folderPath):
	for name in files:
		if name.endswith('_ch2.tif'):
			filePath = os.path.join(dirpath, name)
			#x,y,z = readVoxelSize(filePath, verbose=False)
			print(filePath)
			#pathList.append(filePath)
			#myConvert(path=filePath, imp=None)
