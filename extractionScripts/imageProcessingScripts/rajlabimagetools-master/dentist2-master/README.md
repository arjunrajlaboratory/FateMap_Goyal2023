Dentist2
========
This repository contains tools to help analyze RNA FISH data from large image scans. Example input and output data are available on Dropbox [here](https://www.dropbox.com/sh/djy9mock6vp5j2d/AADqZDtkAEROWK7Jh85oCE68a?dl=0). 

For analysis of immunofluorescence data, consider trying the dentist2IF GUI. Instructions are available [here](https://github.com/arjunrajlaboratory/dentist2/wiki/Dentist2-Immunofluorescence-GUI)

Table of contents
=================
* [Dentist2](#dentist2)
* [Dependencies](#dependencies)
* [Expected input](#expected-input)
* [Quick start](#quick-start)
* [Stitching](#stitching)
* [D2ThresholdGUI](#d2ThresholdGUI)
  * [Main axes](#main-axes)
  * [Thumbnail axes](#thumbnail-axes)
  * [Threshold axes](#threshold-axes)
  * [Adding and deleting masks and nuclei](#adding-and-deleting-masks-and-nuclei)
* [Importing nuclei masks](#importing-nuclei-masks)
* [Loading single-channel .tif scan](#loading-single-channel-tif-scans)
* [Troubleshooting](#troubleshooting)

Dependencies
=============
Dentist2 uses the Bio-Formats MATLAB toolbox (bfmatlab) to read .nd2 files. This software can be downloaded from the [Open Microscopy Consortium](https://www.openmicroscopy.org/) [here](https://docs.openmicroscopy.org/bio-formats/6.3.1/users/matlab/). Dentist2 also uses a few functions from the [rajlabimagetools repository](https://github.com/arjunrajlaboratory/rajlabimagetools). 

Before running Dentist2, please download bfmatlab and rajlabimagetools and add them to your path by typing the following in your MATLAB console. 

```matlab
>>addpath('~/path/to/dentist2/')
>>addpath('~/path/to/bfmatlab/')
>>addpath('~/path/to/rajlabimagetools/')
>>savepath
```
Of note, Dentist2 was written and tested using MATLAB version 2020b. We have not tested Dentist2 using versions prior to 2018b. In addition, Dentist2 uses functions from a variety of built-in MATLAB toolboxes. We recommend installing all available toolboxes when downloading/updating MATLAB. 
 
Expected input
==============
Dentist2 expects images from a **tiled or pre-stitched rectangular scan** in a **single z-plane**. There should be one DAPI channel and one or more FISH channels. The FISH channels may be named 'YFP', 'GFP', 'CY3', 'A594', 'CY5', 'A647' '700' or 'CY7'. If you'd like to process other channels, you may need to modify the map in +d2utils/readND2Channels.m or load single-channel .tif scans as described [here](#loading-single-channel-tif-scans).  

At the moment, Dentist2 identifies nuclei based on DAPI signal and assigned RNA FISH spots to the nearest nucleus. In the future, we may add the option to input whole-cell segmentations and assign spots to the surrounding cell.

Quick start
============
Change your working directory to the folder containing your scan file (`>>cd('~/path/to/scan/'`). If you have multiple scans you wish to analyze with Dentist2, we recommend moving each scan file into separate folders. If your scan needs to be stitched, launch the stitching GUI by typing the following into your console. 
```matlab
>>h = d2stitchingGUI(scanDimensions, 'scanFileName');
```
The scanDimensions should be formatted as \[number of rows, numbers of columns\]. Use the stitching GUI to select the scan layout and control points as described [below](#stitching). When you close the GUI window, two files will be written to your working directory: 'scanSummary.txt' and 'tilesTable.csv'.

**Annotated screenshot of stitching GUI**

<img src="https://github.com/arjunrajlaboratory/dentist2/blob/master/diagrams/d2StitchingGUI.png" width="1100">

**Example scan patterns**

<img src="https://github.com/arjunrajlaboratory/dentist2/blob/master/diagrams/stitchGUIeg.png" width="500">

Next, launch the thresholding GUI by typing the following:
```matlab
>>h = launchD2ThresholdGUI();
```
or if you're loading a pre-stitched scan:
```matlab
>>h = launchD2ThresholdGUI('preStitchedScan', 'path/to/scanFile.nd2');
```
When first launching the threshold GUI, it may take several minutes for the software to stitch the scan, segment and identify nuclei, and identify RNA FISH spots.  The processed data are automatically saved to your working directory and can be loaded more quickly if you need to close MATLAB and restart the GUI. In addition, if the GUI handle (h) is still in your workspace, you can relaunch the GUI by typing `>>h.relaunchGUI;` without having to reload the data.

As described [below](#d2ThresholdGUI), use the threshold GUI to adjust the spot intensity threshold and mask, add, or delete erroneous spots and cells. When you are finished, you can export the data into a summarized table of spots per cell ('spotsSummary.csv') by clicking on the export button. When you close the GUI window, data for all spot calls, nuclei and masks will be saved to 'spots.csv', 'nuclei.csv' and 'masks.csv', respectively.

**Annotated screenshot of threshold GUI**

![thresholdGUI](https://github.com/arjunrajlaboratory/dentist2/blob/master/diagrams/D2threshGUI_2.png)

Note these hotkeys:
- 'up arrow': move up centroid list
- 'down arrow': move down centroid list
- 'z': zoom-mode
- 'p': pan-view mode
- 's': overlay spots checkbox
- 'n': overlay nuclei checkbox
- 'c': scatter checkbox
- 'a': mask all spots checkbox
- 'w': add spot mask
- 'x': add cell mask
- 'd': delete mask
- 'shift-s': save tables
- 'shift-e': export spots summary

Also note, if you click on one of the add mask, add cells, or delete cells buttons, but don't want to complete the operation, try clicking the Esc key. Alternatively, click on an empty area in the main axes then press Enter.

Stitching
==========
Use the scan pattern checkboxes &#9313; to select the pattern corresponding to how your scan images were acquired. If unsure, use the show new positions button &#9314; to display a random set of tiles consistent with the selected pattern. Note that some patterns may place the same tiles at the same positions in a scan, and it can be helpful to check several scan positions. 

Once you've chosen your scan pattern, use the show new positions button &#9314; or textbox underneath to select a position with plenty of nuclei on the border of the images. Click select row control points &#9315; to launch MATLAB's [cpselect](https://www.mathworks.com/help/images/ref/cpselect.html) tool. This should create a new figure window. Select matching pixels on the top and bottom borders of the images like so: 

![cpSelect Row](https://github.com/arjunrajlaboratory/dentist2/blob/master/diagrams/cpSelectRow.png)

Close the window when finished. Click select column control points &#9315; and select matching pixels on the right and left borders of the images like so: 

![cpSelect Col](https://github.com/arjunrajlaboratory/dentist2/blob/master/diagrams/cpSelectColumn.png)

Close the window when finished.

You may click on preview stitch &#9316; to view a stitched image based on the stored control points. If you are satisfied with the stitch, you may exit out of the GUI. Otherwise, you can change scan positions and select additional control points. 

When you exit out of the GUI, the files 'tilesTable.csv' and 'scanSummary.txt' will be saved to your working directory. TilesTable.csv contains the relative coordinates of all tiles in your scan.

D2ThresholdGUI
==============

Main axes
---------
### Scatter
When starting the threshold GUI, the main axes will display a scatter plot of the spatial coordinates of all identified nuclei. Each nucleus (i.e. point) will be colored according to the number of spots assigned to it. You can change the color scheme by selecting a different colormap in the colormap drop-down menu &#9318;. The nuclei colors will adjust automatically as you change the FISH channel, modify the spot intensity threshold or create/delete masks. 

If your scans is larger than 20,000 pixels in any dimension, the main axes will open with a zoomed-in view of your scan. The position of this zoomed-in view will be indicated in the thumbnail axes. You may zoom-out by 2X by right-clicking (or control-clicking on a Mac) the main axes. Note that the more zoomed-out view may slow down the GUI.

You can toggle between the scatterplot and image overlay by selecting the scatter checkbox &#9316;. However, if your current view is >64,000,000 pixels (e.g. >8,000 x 8,000), Dentist2 will not show the image overlay.  
### Image overlay
When toggling off the scatterplot, the main axes will show an overlay of DAPI and the current FISH channel cropped to your current view (note the exception above for very large views). You can use the sliders &#9315;.   


Thumbnail axes
---------------

For the thumbnail axis &#9313; the color of each square reflects the density of nuclei in the corresponding region of the scan. The black rectangle marks the zoomed-in region plotted in the main axis &#9312;. You can click and drag the black rectangle to move around the scan. Alternatively, you can double-click to center the view around the selected position. 

Threshold axes
--------------

The threshold axis &#9323; plots spot fluorescence intensity (x-axis) vs the rank of the intensity values (y-axis; log-scale). The blue vertical line indicates the brightness threshold for valid spots. The brightness threshold is also written in the text box in the upper right of the plot. Only spots whose brightness is ≥ this threshold are plotted as open circles in the main axis and exported to spotsSummary.csv. 

You can manually adjust the spot brightness threshold by dragging the blue vertical line or by typing an integer value into the text box in the upper right of the plot. To zoom into the plot, click the "zoom" button in the upper right then click on the plot and drag horizontally. To zoom-out, double click on the plot. 

By default, the plot includes all identified spots, even those which you mask. To exclude masked spots, click the "filter masked spots" button in the top right. Click the button again to return to the default plotting function.

Adding and deleting masks and nuclei
---------------------------------
### Spot masks
Use spot masks to filter erroneously labeled FISH spots. To add a spot mask, click the "add spot mask" button, draw a mask around the spots you wish to mask, then press Enter. To delete the mask, click the "delete mask" button, click inside each mask you wish to delete, then press Enter.

### Cell masks
Use cells masks to filter erroneously labeled nuclei. To add a cell mask, click the "add cell mask" button, draw a mask around the cells you wish to mask, then press Enter. All spots assigned to masked nuclei will be re-assigned to the nearest valid nuclei (within the maximum cell radius). The nuclei.csv table will store masked nuclei however, masked nuclei will not appear in the main axes and will not be exported to spotsSummary.csv. 

To delete cell masks, click the "delete mask" button, click inside each mask you wish to delete, then press Enter. Spots in your current view will be re-assigned to valid nuclei. 

### Delete cells
Use delete cells to permenantly delete labeled nuclei. Click the "delete cells" button, click on the nuclei wish to delete, then press Enter. As with masking, when you delete a nucleus, all spots assigned to the nucleus will be re-assigned to the nearest valid nucleus. Unlike masking, you cannot undo a nucleus deletion. Instead you must add a new nucleus. 

### Add cells
You can use the combination of "delete cells" and "add cells" to divide 2+ cells that were labelled as one. To add nuclei, click the "add cells" button, click on the locations you wish to label nuclei, then press Enter. All spots in your current view will be re-assigned to the nearest valid nuclei. Note that unlike nuclei identified by dentist2 initially, manually labeled nuclei will have an area of 0 in the nuclei.csv file. 

**Please note, if you click on one of the add mask, add cells, or delete cells buttons, but don't want to complete the operation, try clicking the Esc key. Alternatively, click on an empty area in the main axes then press Enter.**

Importing nuclei masks
===============================
If you use another algorithm to segment nuclei, you can input the nuclei masks into dentist2 using the 'nucleiMasks' parameter. The expected input is either a label matrix (2D image with the position of each mask labeled with a unique integer; see example [here](https://www.dropbox.com/s/f053aeddvd40vbd/RGFullScanDapi_cp_masks.tif?dl=0)) or a text file containing, per row, the boundary coordinates of each mask (i.e. mask outlines; see example [here](https://www.dropbox.com/s/cjcqwycls59qwul/RGFullScanDapi_cp_outlines.txt?dl=0)). Input the mask file when running launchD2ThresholdGUI() as below:
```matlab
>>h = launchD2ThresholdGUI('nucleiMasks', 'path/to/nucleiLabelMatrix.tif');
```
Dentist2 also includes an option to "stitch" masks from separate, overlapping tiles of a larger scan. To use this feature, use the 'cellPose' parameter and specify a folder containing the mask outlines for each tile (this feature currently does not support tiled label matrices currently). You will also need to use the 'cellPoseTileTable' parameter and specify the name of a file that lists the position of each tile. For example: 
```matlab
>>h = launchD2ThresholdGUI('cellPose', 'path/to/cp_directory/', 'cellPoseTileTable', 'cellPoseTilePositions.csv');
```
An example of the cellPoseTilePositions.csv file can be found [here](https://www.dropbox.com/s/28udv1pdlmy3wis/cellPoseTilePositions.csv?dl=0). 

You may use +d2utils/splitStitchedScan.m to generate both the overlapping image tiles from a pre-stitched scan and the cellPoseTilePositions.csv file. For each outline in tile i, dentist2 will determine if the outline overlaps any outlines from previous neighboring tiles, and if so, the overlapping outlines will be merged. This is to avoid duplicating nuclei that fall on tile boundaries. If your scan has a lot of nuclei, this step may take a very long time. You may be better off just running the segmenting algorithm on the larger (stitched) image.  

Note that if you ran your segmenting algorithm on resized images (to save on memory), you can specify the resize factor when running launchD2ThresholdGUI() as below:
```matlab
>>h = launchD2ThresholdGUI('cellPose', 'path/to/cp_outlines.txt', 'maskResizeFactor', 4); %Indicates a resize factor of 4. 
```
or 
```matlab
>>h = launchD2ThresholdGUI('cellPose', 'path/to/cp_directory/', 'cellPoseTileTable', 'cellPoseTilePositions.csv', 'maskResizeFactor', [2 4]); %Indicates a vertical resize factor of 2 and horizontal resize factor of 4. 
```
Loading single-channel .tif scans
=================================
To load pre-stitched, single-channel .tif files into dentist2, first create a cell array containing the file names (e.g. `preStitchedScanFiles={'alexa.tif', 'cy.tif', 'dapi.tif', 'trans.tif'};`). Next, create a cell array that specifies the channel type (either 'FISH', 'dapi', or 'other') for each of the .tif files (e.g. `channelTypes = {'FISH', 'FISH', 'dapi', 'other'};`). Dentist2 will find spots for all image files with channelType == 'FISH' and will find nuclei for the image file with channeltype == 'dapi'. Currently, dentist2 expect only one image file with channeltype == 'dapi'. Image files with channelType == 'other' will be loaded into dentist2 for viewing but will not be used for finding spots. Input the image files list and channelTypes list when running launchD2ThresholdGUI() as below:
```matlab
>>h = launchD2ThresholdGUI('preStitchedScanFilelist', preStitchedScanFiles, 'channelTypes', channelTypes); 
```  
An example dataset and script for loading the data into dentist2 can be found [here](https://www.dropbox.com/sh/yjzw0776ttac1qj/AAD3dmmnF8k8G9FWb2sEdmJ_a?dl=0).

Troubleshooting
===============
* Problem: "Error using X (line #). Undefined function 'A' for input arguments of type 'B'." 
  * Possible Reason: Some MATLAB toolkit is missing. Your MATLAB installation is trying to use a different function named 'A' than what was used for writing Dentist2. 
  * Possible Solution: If you can figure out which toolkit is missing, you can install that toolkit. Alternatively, you can try reinstalling MATLAB with all available toolkits.
* Problem: "Java exception occured: java.lang.OutOfMemoryError: Java heap space". 
  * Possible Reason: You are trying to read a very large image or scan into MATLAB using bfmatlab functions and your Java heap size is set too low. 
  * Possible Solution: Increase your Java heap size by selecting Preferences > MATLAB > General > Java Heap Memory. Instructions available [here](https://www.mathworks.com/help/matlab/matlab_external/java-heap-memory-preferences.html). 
* Problem: Dentist2 is missing FISH spots that are clearly present in the image, mostly at the cell periphery. 
  * Possible Reason: Dentist2 uses a "maximum distance from nucleus" parameter for filtering spots. The default is 200 pixels. If your cells are very large or your images were acquired at very high magnification/resolution, you may need to increase the max distance parameter. 
  * Possible Solution: You can adjust the maximum distance parameter using the method changeMaxDistance in the d2MainAxesController. For example, by typing `h.mainAxesCntrlr.changeMaxDistance(300);` into your MATLAB console, you can set the maximum distance to 300 pixels.
* Problem: Dentist2 is missing FISH spots that are clearly present in the image. Oddly, this is happening in blocks like so: 
![missingSpotsExample1](https://github.com/arjunrajlaboratory/dentist2/blob/master/diagrams/missingSpotsExample1.png)
  * Possible Reason: When finding spots, dentist2 automatically sets a minimum intensity threshold in order to limit the size of the spotsTable (see the function d2utils.findSpotsaTrous.m). The default is max(spotIntensity)/4. If your FISH signal intensity is low and there is bright schmutz in your images, then your FISH signal may fall below the minimum intensity threshold. Also note that denists2 finds spots in blocks in parallel to speed up processing. Thus different blocks may have different max(spotIntensity) values. 
  * Possible Solution: You can decrease the minimum intensity threshold for a specific channel using the method changeSpotFilterThreshold in the d2MainAxesController. For example, by typing `h.mainAxesCntrlr.changeSpotFilterThreshold('cy', 10);`  into your MATLAB console, you can set the minimum intensity threshold for 'cy' to max(spotIntensity)/10 (assuming h is your GUI object). Alternatively, you can change the minimum intensity threshold using the 'aTrousMinThreshFactor' parameter in the launchD2ThresholdGUI.m function, for example: `h = launchD2ThresholdGUI('aTrousMinThreshFactor', 10)`. 
    * Note that if you close the GUI and rerun `h = launchD2ThresholdGUI()` it will set the minimum intensity threshold to it's default max(spotIntensity)/4. If available, `h = launchD2ThresholdGUI()` will load the spots.csv table even if those spots were found using a different intensity threshold. However, if you try to refind spots (because spots.csv is missing or because you ran subtractBackground as described below), dentist2 will use the default threshold (or whatever value is stored in h.spotTable.threshFactor).
* Problem: My scan was acquired at 20X magnification and Dentist2 is doing a terrible job at identifying nuclei. 
  * Possible Reason: Dentist2 effectively uses 2 parameters to identify nuclei: 1. the minimum nucleus area in pixels and 2. the sensitivity parameter for ![MATLAB's adaptive thresholding function] (https://www.mathworks.com/help/images/ref/adaptthresh.html). The default values for these parameters were selected based on scans acquired at 60X magnification. 
  * Possible Solution: You can adjust the minimum nucleus area and  the adaptthresh sensitivity parameter using the adjustNucleiParameters method in the d2MainAxesController. For example, by typing `h.mainAxesCntrlr.adjustNucleiParameters(200, 'sensitivity', 0.2);` into your MATLAB console, you can set the minimum nucleus area to 200 pixels and the adaptthresh sensitivity to 0.2 (these values worked OK for 20X images of Calu-3 cells). Note that the DAPI signal intensity and contrast also affect nucleus masking but currently Dentist2 does not include methods for adjusting these features. 
  * Possible Solution 2: I've found that ![Cellpose](https://github.com/MouseLand/cellpose) does a much better job at masking nuclei in lower mag images (e.g. 20X) than the adaptive thresholding algorithm uses in Dentist2. You can find detailed instructions for running Cellpose ![here](https://cellpose.readthedocs.io/en/latest/) and a tutorial based on my own experience in the ![d2IF wiki](https://github.com/arjunrajlaboratory/dentist2/wiki/Dentist2-Immunofluorescence-GUI). 
* Problem: 
  * Possible Reason:
  * Possible Solution:

