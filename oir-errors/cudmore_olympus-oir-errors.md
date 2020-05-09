
## Subject:

Errors importing multi-channel deconvolved Olympus .oir files into ImageJ/Fiji using ImageJ/Fiji with Bioformats

## Description of problem

I am trying to use ImageJ/Fiji to open a 2-3 color channel Olympus .oir file that was saved on an Olympus FV300 scope after the Olympus software created a new deconvolved .oir file.

It is **working** in the 2017 lifeline version but **does not work** in newer ImageJ/Fiji

**Working** in ImageJ/Fiji lifeline:
	Version: 2.0.0-rc-59/1.51n
	Date: 2017-02-23T11:35:06-600
	Bioformats: Release 5.5.1, build date 25 May 2017
	
**Does not work** in a current version of ImageJ/Fiji:
	Version: 2.0.0-rc-69/1.52v
	Date: 2018-12-04T11:30:09+000
	Bioformats: Release 6.4.0, build darw 11 March 2020

## Example images

Here are example original and deconvolved .oir files:

 - original .oir: https://ucdavis.box.com/s/ppi2jncn2g0otg82br97uka6wj27c57d
 - deconvolved .oir: https://ucdavis.box.com/s/as8px5iwfgicn01eet81pkqfg4gba04o

## Symptoms:

In the current version of ImageJ/Fiji (direct download and after using update), dragging and dropping the deconvolved .oir causes channel 1 to be repeated into channel 2. We don't get the actual channel 2 image/stack. **For a three channel deconvolved .oir, the 3rd channel seems to be fine?** Same symptom if the deconvolved .oir file is imported using the menu "Plugin : Bio-Formats : Bio-Formats-Importer".

 - Drag and drop the **original** .oir file works in **both** the lifeline and current versions of ImageJ/Fiji.
 - Drag and drop the **deconvolved** .oir works in **lifeline ImageJ/Fiji**.
 - Drag and drop the **deconvolved** .oir **does not** work in the current version of ImageJ/Fiji, we get a copy of channel 1 in channel 2. Channel 3 (if it exists) seems fine.
 
I see the same problem/error in a ImageJ/Fiji Jython script. This script **works** in the lifeline version but **does not work** in the current version of ImageJ/Fiji. Specifically, in the newer version, channel 2 is a copy of channel 1 and for some reason Channel 3 seems correct ...

```
# Author: Robert Cudmore
# Date: 20200424

import os
from ij import IJ
from ij.plugin import ChannelSplitter
from loci.plugins import BF

#
# after you download these images, please set the paths to your local copies ...
original_OIR = '/Users/cudmore/box/share/deconvolded-oir-errors/8_5ADVMLEG1L1.oir'
deconvolved_OIR = '/Users/cudmore/box/share/deconvolded-oir-errors/8_5ADVMLEG1L1.oir'

#
# choose what to prcoess, either oringal or deconvolved
#filePath = original_OIR
filePath = deconvolved_OIR

# set up a base name to save channels 
# We will append one of ('_ch1.tif', '_ch2.tif', '_ch3.tif') to saveFileNoExtension
folderPath, fileName = os.path.split(filePath)
fileNameNoExtension = os.path.splitext(fileName)[0]
saveFileNoExtension = os.path.join(folderPath, fileNameNoExtension)

# open the .oir
imps = BF.openImagePlus(filePath) # list with one imp
if len(imps)==1:
	imp = imps[0]

# split the channel
channelArray = ChannelSplitter.split(imp) # returns an array of imp

for channelIdx, channelImp in enumerate(channelArray):
	saveFilePath = saveFileNoExtension + '_ch' + str(channelIdx+1) + '.tif'
	print('saving saveFilePath:', saveFilePath)
	IJ.save(channelImp, saveFilePath)
```

## Conclusion

Can anyone give some feedback?

I know it is annoying as this kind of question spans multiple groups including Olympus, ImageJ/Fiji, BioFormats, and me ... the end user.

This is not the first time I have posted about the inter-operability between Olymous .oir, ImageJ/Fiji, and Bioformats, please see:

[problems-opening-olympus-oir-files-using-bio-formats](https://forum.image.sc/t/problems-opening-olympus-oir-files-using-bio-formats/24747)

[problems-opening-olympus-oir-line-scan-files](https://forum.image.sc/t/problems-opening-olympus-oir-line-scan-files/24957/3)

## My Opinion

Proprietary/commercial file formats need to come to an end. The open-source community is now providing tools that greatly surpass these propietary/commercial formats. Proprietary/commercial formats are a waste of time for corporations like Olympus, open source developers like ImageJ/Figi, tool developers like BioFormats/Loci, and end users like me.

## Signature

Thanks,

Robert Cudmore

_____________________________________________
Assistant Professor
Department of Physiology and Membrane Biology
School of Medicine - University of California at Davis
https://robertcudmore.org

