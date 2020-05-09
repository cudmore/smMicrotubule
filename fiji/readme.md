These .py Jython files are to be run from within [Fiji](https://fiji.sc/).

Drag and drop a file onto Fiji and in the editor select menu run (or command+r)

From a command line on macOS, use the following to run `samiVolume.py`

```
/Applications/Fiji.app/Contents/MacOS/ImageJ-macosx samiVolume.py batch=../analysis/wt-female.txt
```

## Notes

menu 'help - update imagej ...'

v1.52v


```
# path to fiji.app we are using
# from java.lang.System import getProperty
#print(getProperty('user.dir'))


# __file__ is not available in Fiji/Jython
#myPath = os.path.dirname(os.path.abspath(__file__))

# __file__ is not available in Jython, we would normally just do
# .   myFilePath = os.path.dirname(os.path.abspath(__file__))
# this is required to do normal Python import <module> and/or from <module> import <function/class>
# see: https://forum.image.sc/t/can-the-script-engine-tell-a-script-in-which-folder-the-script-is-located/9805/4
# see: https://imagej.net/Jython_Scripting_Examples#Importing_other_.py_scripts_.28modules.29
from ij.io import OpenDialog, DirectoryChooser
myFilePath = OpenDialog.getLastDirectory() # would be __file__ in Python
sys.path.append(myFilePath)
```

### Importing custom .py code

see: https://forum.image.sc/t/can-the-script-engine-tell-a-script-in-which-folder-the-script-is-located/9805/5

```
import inspect
print(inspect.getsourcefile(lambda:0))
myFilePath = os.path.dirname(os.path.abspath(inspect.getsourcefile(lambda:0)))
print('    myFilePath:', myFilePath)
sys.path.append(myFilePath)

import simpleFileInfo
```

### ImagePlus

		'''
		myLut = channelImp.getLuts()
		print('    myLut:', myLut)
		'''
		
		'''
		# from (ImagePlus.GRAY8, ImagePlus.GRAY16, ImagePlus.GRAY32, ImagePlus.COLOR_256 or ImagePlus.COLOR_RGB)
		myType = channelImp.getType()
		print('    myType:', myType)
		'''
		
		'''
		ip = channelImp.getProcessor()
		print('    ip before reset:', ip)
		ip.resetMinAndMax()
		print('    ip after  reset:', ip)
		ip = channelImp.setProcessor(ip)
		'''
