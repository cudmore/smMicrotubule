from ij import IJ, WindowManager

for imp in map(WindowManager.getImage, WindowManager.getIDList()):
	imp.close()