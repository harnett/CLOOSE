# Welcome to CLOOSE!

To stream data to and from Matlab, the user needs to have the Instrument Control Toolbox installed.

## CLOOSE <img src="Logo/Logo.png" width="500" title="CLOOSE" alt="CLOOSE" align="right" vspace = "100">

Pipeline for running BCI experiments and online analysis of 1 and 2 photon imaging data

Copyright (C) 2025  Massachusetts Institute of Technology

## Paper Abstract

Brain-Computer Interfaces (BCI) have catalyzed advancements in both clinical applications and basic neuroscience research. However, technical barriers such as steep learning curves and complex synchronization requirements often impede their widespread adoption. In response to the increasing demand for optical closed-loop experiments and to the technical challenges associated with their implementation, we introduce CLOOSE (Closed-loop Optical Open-Source Experiments), a versatile platform designed to seamlessly perform closed-loop experiments on optical data. CLOOSE can easily interfaced with any frame-based imaging acquisition systems via a TCP/IP connection for data stream. Through benchmark tests, we validate CLOOSE's capability in ensuring real-time accurate image registration on both dense and sparse data, signal processing, and analysis at imaging frequencies of 1 kHz and above, making CLOOSE the first optical BCI compatible with voltage imaging. Throughout the paper we showcase CLOOSE's versatility in supporting several neurofeedback paradigms (from single neuron to population dynamics), multiplane imaging, and online z-tracking. CLOOSE functionality can be easily exploited by users with minimal coding ability thanks to an intuitive graphical user interface used to parametrize experiments and visualize tracking of CLOOSE performance. CLOOSE aims to streamline, standardize, and democratize online image analysis and neurofeedback BCI experiments by significantly lowering the entry level for performing a wide range of technically challenging in vivo imaging experiments. 
 

This code was written by Valerio Francioni.
For support, please open an [issue](https://github.com/CLOOSE_MS/issues).

### CITATION

If you use this package in your research, please cite the [paper](BIORXIV):

Valerio Francioni, Anna Beltramini, Linlin Z Fan and Mark T. Harnett (2025). CLOOSE  : An open-source platform to perform optical closed-loop experiments.

The online motion correction algorythm is based on [this paper](https://opg.optica.org/ol/abstract.cfm?uri=ol-33-2-156)

## Before you start
Make sure you have the following toolboxes installed with your Matlab 
1. DSP System Toolbox
2. Signal Processing toolbox
3. Image Processing Toolbox
4. If you are planning on running your BCI experiments using visual stimuli as feedback, make sure you have Psychtoolbox installed. Details can be found [here](http://psychtoolbox.org/download).

## Getting started
## Step 1 -  Make sure you can transfer imaging data via a TCP/IP connection internally
1. Open two instances of MATLB on the same machine
2. In the TCP_IP folder open stream_send_data.m in one instance of MATLAB
3. In the TCP_IP folder open stream_receive_data.m in the other instance of MATLAB
4. In the stream_send_data.m script, set the path to Generate_images folder (e.g. YourPath2GitHub\GitHub\CLOOSE_MS\Generate_images)
5. For this step, leave the IPv4 variable unchanged to 127.0.0.1 (Localhost) in both the stream_send_data.m and the stream_receive_data.m.
6. Run stream_send_data.m first, then run stream_receive_data.m
7. CONGRATULATIONS! YOU ARE STREAMING DATA BETWEEN TWO INSTANCES OF MATLAB 

## Step 2 -  Make sure you can transfer imaging data via a TCP/IP connection internally
1. Open two instances of MATLB on two different computers, connected to the same Network (ethernet will be faster but Wifi will work too)
2. In the TCP_IP folder open stream_send_data.m on PC
3. In the TCP_IP folder open stream_receive_data.m in the other PC
4. In the stream_send_data.m script, set the path to Generate_images folder (e.g. YourPath2GitHub\GitHub\CLOOSE_MS\Generate_images)
5. On the PC _sending_ data, set the IPv4 as the IPv4 of the receiving PC (e.g., 10.93.6.184). NB: You can find the IPv4 adress by typing ipconfig/all in your command window.
6. 5. On the PC _receiving_ data, set the IPv4 as the IPv4 of the sending PC (e.g., 10.18.1.121). 
8. Run stream_send_data.m first, then run stream_receive_data.m
9. CONGRATULATIONS! YOU ARE STREAMING DATA BETWEEN TWO DIFFERENT PCs.



## If you managed to run these two steps successfully you are 90% there. 
CLOOSE is setup to act at the receiving end of this process. All you have to do, is to stream data to it. The rest of the analysis can easily be setup using the GUI (next steps). All you have to do now, is to setup your acquisition machine (the PC where you are acquiring imaging data), in the same way as the stream_send_data.m is setup. 
## Critical steps for this are: 
1. Open a TCP/IP connection on your image acquisition device in the same way it's done in the stream_send_data.m script 
```	
 tcpipServer = tcpip(IPv4, port, 'NetworkRole', 'client', ...
    		'InputBufferSize', (xpixels*ypixels*buffsz), 'Terminator', 'CR/LF');
fopen(tcpipServer);
```

2. Identify in your code where the variable encoding for your image is stored. Normally acquisition devices plot your imaging data so you can visualize them online. Look for plotting functions in your code and you should be able to find your image. 
3. Vectorize your image in the same way it's done in the stream_send_data.m script. In this case the vectorized image is the variable 'stack'. NB The code below vectorizes a single frame 
```	
 vect_img = reshape(stack(:, :, iframe), [], 1);
```
4. Open CLOOSE and run your experiments.

### Using the GUI

<img src="https://www.suite2p.org/static/images/multiselect.gif" width="800" alt="selecting multiple ROIs in suite2p with Ctrl"/>


The suite2p output goes to a folder called "suite2p" inside your save_path, which by default is the same as the data_path. If you ran suite2p in the GUI, it loads the results automatically. Otherwise, you can load the results with File -> Load results or by dragging and dropping the stat.npy file into the GUI.

The GUI serves two main functions:

1. Checking the quality of the data and results.
	* there are currently several views such as the enhanced mean image, the ROI masks, the correlation map, the correlation among cells, and the ROI+neuropil traces
	* by selecting multiple cells (with "Draw selection" or ctrl+left-click), you can view the activity of multiple ROIs simultaneously in the lower plot
	* there are also population-level visualizations, such as [rastermap](https://github.com/MouseLand/rastermap)
2. Classify ROIs into cell / not cell (left and right views respectively)
	* the default classifier included should work well in a variety of scenarios.
	* a user-classifier can be learnt from manual curation, thus adapting to the statistics of your own data.
	* the GUI automatically saves which ROIs are good in "iscell.npy". The second column contains the probability that the ROI is a cell based on the currently loaded classifier.

Main GUI controls (works in all views):

1. Pan  = Left-Click  + drag
2. Zoom = (Scroll wheel) OR (Right-Click + drag)
3. Full view = Double left-click OR escape key
4. Swap cell = Right-click on the cell
5. Select multiple cells = (Ctrl + left-click) OR (SHIFT + left-click) AND/OR ("Draw selection" button)

You can add your manual curation to a pre-built classifier by clicking "Add current data to classifier". Or you can make a brand-new classifier from a list of "iscell.npy" files that you've manually curated. The default classifier in the GUI is initialized as the suite2p classifier, but you can overwrite it by adding to it, or loading a different classifier and saving it as the default. The default classifier is used in the pipeline to produce the initial "iscell.npy" file.


## Outputs

~~~~
F.npy: array of fluorescence traces (ROIs by timepoints)
Fneu.npy: array of neuropil fluorescence traces (ROIs by timepoints)
spks.npy: array of deconvolved traces (ROIs by timepoints)
stat.npy: array of statistics computed for each cell (ROIs by 1)
ops.npy: options and intermediate outputs
iscell.npy: specifies whether an ROI is a cell, first column is 0/1, and second column is probability that the ROI is a cell based on the default classifier
~~~~

# License

Copyright (C) 2023 Howard Hughes Medical Institute Janelia Research Campus, the labs of Carsen Stringer and Marius Pachitariu.

**This code is licensed under GPL v3 (no redistribution without credit, and no redistribution in private repos, see the [license](LICENSE) for more details).**

# Coming soon
See this **youtube [thread](SOME YOUTUBE LINK)** for GUI demonstrations.

### Logo
Logo was designed by ChatGPT :(:
