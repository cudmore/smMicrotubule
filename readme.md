## Todo

- make a mask for nuclear (ch1) and remove those pixels from both (full, ring, eroded)
- try and use gaussian followed by erosion for maksing (not for skeleton)
- make sure samiViewer2.py just takes an output folder
- implement selection of n random (e.g. 100) measurements from skeleton, this requires to go through and figure out smallest n

- my masks are always 1-2 slices too big in z. Presumably due to z-spread (point spread function). Maybe remove 1-4 (probably 3-4 slices) at the top and bottom of the stack

- this points to a bigger problem with the shared santana confocal scope. Is it in alignment???

## Overview

1) convert raw data

Analysis

1) specify path to save/view everything in file gAnalysisPath.py
2) ./batch-samiAnalysisParallel.sh 
3) python samiVolume2.py # slowest
4) python samiDensity.py
5) run jupyter notebook densityAnalysis.ipynb
6) view with `python samiViewer2.py ../analysis/wt-female.txt

samiVolume3/
	samiVolume2.py used 
		erodedIterations0 = 3
		erodedIterations = 3
		
## Install aics and bImPy togethere

Need to use conda

```
# see: #https://github.com/AllenInstitute/aics-segmentation/blob/master/docs/installation_mac.md
cd ~/Sites
git clone https://github.com/AllenInstitute/aics-segmentation.git

cd ~/Sites/smMicrotubule

conda create -n sami_env python=3.6

conda activate sami_env

conda install nb_conda


pip install numpy
pip install itkwidgets==0.14.0

pip install -e ../aics-segmentation.

# install bimpy
pip install -e /Users/cudmore/Sites/bImPy/.

# downgrade tifffile (for aics)
pip install tifffile==0.15.1
```

remember

```
pip install seaborn statannot
```

## A) Workflow

### 1) Convert oir file to .tif and split channels

Drag and drop [fiji/samiConvert.py](fiji/samiConvert.py) onto Fiji. Manually specify path by editing the folderPath ...

```
# specify path here
# .endswith('_3ADVMLEG1L1.oir')
folderPath = '/Users/cudmore/box/data/sami/data/200421'
```

### 2) Add each file path to corresponding  [analysis/](analysis/) .txt file

Add **full path** for each new _ch2.tif file to corresponding file in the [analysis/](analysis/) folder. There are 4 [analysis/](analysis/) batch .txt files.

```
analysis/wt-female.txt
analysis/wt-male.txt
analysis/ko-female.txt
analysis/ko-male.txt
```

To simplify things, use [python/getFileList.py](python/getFileList.py)to quickly generate a list of all files for a given folder (and all subfolders). Just copy/paste the output into the corresponding analysis/ batch .txt file.


### 3) Run [python/samiAnalysis.py](python/samiAnalysis.py) for each analysis/ batch `.txt` file

This will use aics segmentation to make a filament mask and then use Skan to make a skeleton

```
cd python
python samiAnalysis.py batch=../analysis/wt-female.txt
python samiAnalysis.py batch=../analysis/wt-male.txt
python samiAnalysis.py batch=../analysis/ko-female.txt
python samiAnalysis.py batch=../analysis/ko-male.txt
```

Or run all 4 of those commands with [python/batch-samiAnalysis.sh](python/batch-samiAnalysis.sh)

```
cd python
./batch-samiAnalysis.sh
```

[python/samiAnalysis.py](python/samiAnalysis.py) will save a file with all the analysis results. For the `analysis/wt-female.txt` batch file, it will be [analysis/wt-female_results.csv](analysis/wt-female_results.csv).

```
	myCellNumber	filename	genotype	sex	branchType	len3d	euclideanDist	tortuosity	path
0	0	1_5ADVMLEG1L1_ch2.tif	wt	female	1	1.442065756	1.259733016	1.144739193	/Users/cudmore/box/data/sami/data/200108/WT_Female/Cell_1/1_5ADVMLEG1L1_ch2.tif
1	0	1_5ADVMLEG1L1_ch2.tif	wt	female	1	0.690498805	0.690498805	1	/Users/cudmore/box/data/sami/data/200108/WT_Female/Cell_1/1_5ADVMLEG1L1_ch2.tif
2	0	1_5ADVMLEG1L1_ch2.tif	wt	female	1	1.503575026	1.187815239	1.26583241	/Users/cudmore/box/data/sami/data/200108/WT_Female/Cell_1/1_5ADVMLEG1L1_ch2.tif
3	0	1_5ADVMLEG1L1_ch2.tif	wt	female	2	1.230391642	1.117784217	1.100741649	/Users/cudmore/box/data/sami/data/200108/WT_Female/Cell_1/1_5ADVMLEG1L1_ch2.tif
4	0	1_5ADVMLEG1L1_ch2.tif	wt	female	2	0.821482703	0.67197502	1.222489941	/Users/cudmore/box/data/sami/data/200108/WT_Female/Cell_1/1_5ADVMLEG1L1_ch2.tif
5	0	1_5ADVMLEG1L1_ch2.tif	wt	female	1	0.658364432	0.318208678	2.068970705	/Users/cudmore/box/data/sami/data/200108/WT_Female/Cell_1/1_5ADVMLEG1L1_ch2.tif
6	0	1_5ADVMLEG1L1_ch2.tif	wt	female	2	0.842093739	0.67197502	1.253162265	/Users/cudmore/box/data/sami/data/200108/WT_Female/Cell_1/1_5ADVMLEG1L1_ch2.tif
7	0	1_5ADVMLEG1L1_ch2.tif	wt	female	1	0.418505475	0.418505475	1	/Users/cudmore/box/data/sami/data/200108/WT_Female/Cell_1/1_5ADVMLEG1L1_ch2.tif
8	0	1_5ADVMLEG1L1_ch2.tif	wt	female	2	2.291680893	1.695334635	1.351757255	/Users/cudmore/box/data/sami/data/200108/WT_Female/Cell_1/1_5ADVMLEG1L1_ch2.tif
```

In this .csv file, each row is a branch segment. The branchType is from [0,1,2] with 0 meaning disconnected, 1 meaning connected on one end, and 2 meaning connected on both ends.

**todo** Go back to Skan documentation to verify this.


**todo**: I need to have `python/samiAnalysis.py` check if its output, `_dvMask.tif` already exists and not reprocess (it is getting slow). I don't think I can easily get unique cell numbers? All files need to be saved with unique names, e.g. `<yyyymmdd_cell<n>.oir`. Where `<yyyymmdd>` is the date it was acquired and `<n>` is a file index appended by Olympus software.

### 3.1) Description of algorithm to generate a Skeleton

This code uses The Allen Institute For Cell Science Segmentation (AICS-Segmentation). The Python source is at [https://github.com/AllenInstitute/aics-segmentation](https://github.com/AllenInstitute/aics-segmentation). We are using the 'Classic Image Segmentation Workflow' as compared to the Deep Learning part.

Algorithm is as follows. Both the **`f3_param`** and **`minArea`** are parameters provided by the user. Basically, start with raw image, contrast adjust and smooth, generate mask and then use [Skan](https://jni.github.io/skan/) to create a skeleton from the mask.

a) Automatically choose a contrast adjustment

```
low_ratio, high_ratio = my_suggest_normalization_param(struct_img0)
intensity_scaling_param = [low_ratio, high_ratio]
struct_img = intensity_normalization(struct_img0, scaling_param=intensity_scaling_param)
```

b) Edge preserving smoothing (similar to median filter)

```
structure_img_smooth = edge_preserving_smoothing_3d(struct_img)
```

c) Convert the contrast adjusted and smoothed image to a binary mask

```
	f3_param=[[0.5, 0.005]]
	bw = filament_3d_wrapper(structure_img_smooth, f3_param)

	scale_x is set based on the estimated thickness of your target filaments.
		For example, if visually the thickness of the filaments is usually 3~4 pixels,
		then you may want to set scale_x as 1 or something near 1 (like 1.25).
		Multiple scales can be used, if you have filaments of very different thickness.
	cutoff_x is a threshold applied on the actual filter reponse to get the binary result.
		Smaller cutoff_x may yielf more filaments, especially detecting more dim ones and thicker segmentation,
		while larger cutoff_x could be less permisive and yield less filaments and slimmer segmentation.
```

```
	minArea=20
	seg = skimage.morphology.remove_small_objects(bw>0, min_size=minArea, connectivity=1, in_place=False)
```

d) Generate a 1-pixel wide skeleton from the mask using [Skan](https://jni.github.io/skan/)

```
retDict0, mySkeleton = myAnalyzeSkeleton(out=out, imagePath=path) 
```

### 4) Volume analysis [fiji/samiVolume.py](fiji/samiVolume.py)

Code is now in

```
# HERE
/Users/cudmore/Sites/smMicrotubules/python/samiVolume2.py 

# NOT HERE
/Users/cudmore/Sites/bImPy/examples/edt/samiVolume2.py 
```

**Do not use Fiji version of this code !!!!**

This is entirly done in Fiji. Using the '3D ImageJ Suite' plugin from [https://imagejdocu.tudor.lu/plugin/stacks/3d_ij_suite/start](https://imagejdocu.tudor.lu/plugin/stacks/3d_ij_suite/start). The '3D ImageJ Suite' is not included in Fiji by default, please follow their install instructions.

Results of [fiji/samiVolume.py](fiji/samiVolume.py) look like this, where each row corresponds to an original .tif file and `pixelCount_1/_2/_3` correspond to the number of **pixels** in the largest three segmented volumes. Use the x/y/z voxel in um/pixel columns (xVoxel, yVoxel, zVoxel) to calculate the volume (in um^3) of the largest segmented object (e.g. `pixelCount_1`). The reason we end up with the 3 largest volumes is so we can manually decide if the segmentation was spotty, this happens when the number of pixels in _2 and _3 are large or close to pixelCount_1.

```
pixelCount_1,	pixelCount_2,	pixelCount_3
451316,	None,	None
375330,	117,	46
743800,	32,	None
307354,	34,	29
646537,	21,	8
```

**remember:** because z-spread is so absurd, maybe ignore it and assume it is 1 in all volume calculations.

**todo:** Add code to calculate the largest volume (pixelCount_1) in um^3 using xVoxel/yVoxel/zVoxel.

**todo:** These volume masks are not saved? 

#### Option 1: Run from a command-line (terminal). Where 'wt-female.txt' is a text file with a list of full .tif file paths

```
cd fiji
/Applications/Fiji.app/Contents/MacOS/ImageJ-macosx samiVolume.py batch=../analysis/wt-female.txt

```

Results are saved in `analysis/wt-female-volume.csv`. Be sure to move this file to a better location. It will be over-written each time a batch analysis is run (e.g. batch=../analysis/wt-female.txt).

####  Option 2: Run all conditions as a volume analysis batch

```
cd fiji
./volumeAnalysisBatch.sh
```

#### Option 3: [NOT WORKING] Run on one file from within Fiji

Drag and drop [fiji/samiVolume.py](fiji/samiVolume.py) into Fiji, open an image and run the code from the editor.

Once done, the results table need to be copy/pasted into a .csv text file. Try using `right-click` to `Save As ... csv`.


## B) Post Analysis using [python/samiPostAnalysis.py](python/samiPostAnalysis.py)

Once [python/samiAnalysis.py](python/samiAnalysis.py) is run on each analysis/ batch .txt file (e.g. analysis/wt-female.txt) we end up with .csv files containing the output. These are saved in the `analysis` folder as follows:

```
analysis/wt-female_results.csv
analysis/wt-male_results.csv
analysis/ko-female_results.csv
analysis/ko-male_results.csv
```

Now we can use the [python/samiPostAnalysis.py](python/samiPostAnalysis.py) as an engine to derive results quickly. This is easiest done using the provided Jupyter notebook in [python/exampleAnalysis.ipynb](python/exampleAnalysis.ipynb). To run this notebook locally, use...

```
cd python
jupyter notebook
```

## C) Calculate density of filaments in (full mask, eroded mask, ring mask)

Use the following, be sure to tweek where it grabs analysis from


```
# edit samiDensity.py
analysisRoot = '/Users/cudmore/Desktop/samiVolume' # remember, _ch1.tif DOES NOT have scale here !!!!!
analysisRoot = '/Users/cudmore/Desktop/samiVolume2' # remember, _ch1.tif DOES NOT have scale here !!!!!

python /Users/cudmore/Sites/bImPy/examples/edt/samiDensity.py
```

This pulls the **3 masks** (full, eroded, ring) from results of samiVolume2.py and the **skeleton** from smMicrotubules/python/samiAnalysis.py

It then calculates the density of filaments in each mask

Save to /Users/cudmore/Desktop/density-results.csv

This is then used by samiPostAanalysis.py

## Utilities

#### [python/readVoxelSize.py](python/readVoxelSize.py)

given a path, recursively generate a list of x/y/z voxel sizes for all .tif files in path and all subfolder.

#### [python/getFileList.py](python/getFileList.py)

given a path, recursively generate a list of all _ch2.tif files in path and all subfolder

## To Do

 - When I analyze skeleton, I need to break analysis down by branch type. final alanalysis might just be node-to-node branches?

https://jni.github.io/skan/complete_analysis.html

 - Look into this 'shape index', is it included in the Skan summary? It is used in Skan complete_analysis

https://scikit-image.org/docs/dev/api/skimage.feature.html#skimage.feature.shape_index

 - [Done] Look into using sesaborn, a wrapper around matplotlib? https://seaborn.pydata.org/

## Open Science Framework (osf) as a data repository

todo: Figure out if I can use their [API](https://developer.osf.io/) to pull data from within Jupyter notebooks.

https://osf.io/dashboard

my user profile

https://api.osf.io/v2/users/z6s3q/

https://api.osf.io/v2/users/z6s3q/nodes/

another way to get my nodes
https://api.osf.io/v2/nodes/48mq7/

link to my human readable page: https://osf.io/48mq7/

public link to human readable: https://osf.io/48mq7/

## Docker

**todo:** Look into Docker compose again. In particular, I want to mount volumes from Docker container onto client desktop so user/client can edit files.

https://hub.docker.com/u/cudmore

or ?

https://hub.docker.com/repository/docker/cudmore/mydocker

### Build and push from my end

```
docker build -t cudmore/mydocker .

# run locally just to make sure it works
docker run -p 8080:8080 cudmore/mydocker

# push to docker cloud
docker push cudmore/mydocker
```

#### For end users and clients

1) Download and install docker [https://www.docker.com/products/docker-desktop](https://www.docker.com/products/docker-desktop).

2) Check the install in a Terminal

```
docker --version
```

What do you see? Should be like 'Docker version 19.03.8, build afacb8b'. If you do not see that then **STOP**. Maybe Quit and rerun Terminal and try again?

3) Pull the `cudmore/smmicrotubule` image

```
docker pull cudmore/smmicrotubule
```

run `cudmore/smmicrotubule` in a web browser

```
docker run -p 8080:8080 cudmore/smmicrotubule
```

The `docker run` command will return a bunch of `//` addresses like this:

```
    To access the notebook, open this file in a browser:
        file:///root/.local/share/jupyter/runtime/nbserver-1-open.html
    Or copy and paste one of these URLs:
        http://826227aca5ec:8080/?token=6d46989842979160a6ba6c3b8cdaab4c93f11245918968dd
     or http://127.0.0.1:8080/?token=6d46989842979160a6ba6c3b8cdaab4c93f11245918968dd
```

Copy/paste one of these http:// addresses into the address bar of a Chrome browser window. If one of them does not work, try another one.

If the Chrome bowser opens correctly, navigate into '/notebooks/exampleAnalysis.ipynb' and use the Jupyter interface to run each cell (e.g. the Run button). 

If any of the cells return an error, send it to me. Some of the cell will produce `RuntimeWarning` in red because of stat tests that have too low n.

When you get to the cell after 'Plot all raw data', in my browser it is labelled 'In [6]', you should get some nice violin plots!!!
 
#### Docker notes

Remember, I am using `.dockerignore`. If I don't pay attention to details here, the docker image gets artificially big (like many GB).

```
__pycache__
data
fiji
old
pdf
tochristophe
```

## Allen Institute for cell science

#### see jupyter notebook here:

https://github.com/AllenInstitute/aics-segmentation/blob/master/lookup_table_demo/playground_filament3d.ipynb

#### Recipe

```
16. Pseudocode of classic image segmentation workflow for alpha tubulin

Input: I (original single channel 3D image stack)
Output: Final_segmentation (binary image of segmentation result)

Constant Parameters:
           normalization_param = [1.5, 8.0]
           F3_param = [[1,0.01]]
           min_size = 20
# pre-processing
I_norm = Auto_Contrast(I, normalization_param)
I_smooth = Edge_Preserving_Smoothing(I_norm)

# apply F3 filter
Seg = Filament3D(I_smooth, F3_param)

# size filtering
Final_segmentation = Size_Filter(Seg, min_size)
```

#### Install

```
ERROR: aicsimageprocessing 0.7.3 has requirement aicsimageio>=3.1.2, but you'll have aicsimageio 0.6.4 which is incompatible.

```
## Misc

 - Convert jupyter notebooks to readthedocs

https://matthew-brett.github.io/nb2plots/worked_example.html#worked-example

## Sami's Data

https://ucdavis.app.box.com/folder/107822482369

