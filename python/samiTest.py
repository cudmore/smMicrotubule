from aicsimageio import AICSImage, omeTifWriter							

path = '/Users/cudmore/box/data/sami/191230/WT_Male/Cell_4_/4_1_5ADVMLEG1L1_ch2.tif'

img = AICSImage(path)

imgDims = img.dims  # returns string "STCZYX"
imgShape = img.shape  # returns tuple of dimension sizes in STCZYX order
#img.size("STC")  # returns tuple of dimensions sizes for just STC

print('imgDims:', imgDims, 'imgShape:', imgShape)

