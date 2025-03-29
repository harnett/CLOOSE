# Welcome to CLOOSE! <img src="Logo/Logo.png" width="500" title="CLOOSE" alt="CLOOSE" align="right" vspace = "100">

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
1. Open a TCP/IP connection on your image acquisition device in the same way it's done in the stream_send_data.m script. xpixels is the number of lines in your image. ypixels is the number of columns in your image. Set buffsz to 20.  
```	
 tcpipServer = tcpip(IPv4, port, 'NetworkRole', 'client', ...
    		'InputBufferSize', (xpixels*ypixels*buffsz), 'Terminator', 'CR/LF');
fopen(tcpipServer);
```

2. Identify in your code where the variable encoding for your image is stored. Normally acquisition devices plot your imaging data so you can visualize them online. Look for plotting functions in your code and you should be able to find your image. 
3. Vectorize your image in the same way it's done in the stream_send_data.m script. In this case the vectorized image is the variable 'stack', which is tranformed into a column vector. NB The code below vectorizes a single frame (iframe)
```	
 vect_img = reshape(stack(:, :, iframe), [], 1);
```
4. Open CLOOSE and run your experiments.


## Using the GUI: 1. Parameters

The CLOOSE GUI is equipped with various settings to accommodate different experimental designs and uses. The user will set their path within the General Panel, and ensure all settings are correct for their specific experiment. The following parameters exist within the CLOOSE GUI: 
### General Panel
_Path_ : Directory where data will be saved 
_Trial Number_ : Readout of current trial number during the session
_% correct_ : Readout of percentage of trials that are correct
_% avers_ : Readout of percentage of trials that have been incorrect/aversive (Might go in future releases)
_Day_ : Current day of training 
### BCI Panel 
_Motion Corr_ checkbox : Opt to use CLOOSE online motion correction
_Stim Baseline_ checkbox : Opt to have stimuli or darkness presented during baseline
Display of Activity Levels : Readout of activity difference / grating angle during session
Trials : Number of trials during BCI session / determines length of BCI session 
Spatial freq. : 
Angle : 
Pixels(x) : Number of pixels being acquired by image acquisition software 
Lines (y) : Number of lines being acquired by image acquisition software
Frames (z) : 
Target : Required activity ratio (?) between the populations to achieve correct trial
I.T.I : Interstimulus Interval after completion of a trial (?)
Fold : Number of fold lines when using Scanbox subframe folding for high frame rates
“Baseline” : Push button for starting run of baseline recording session
“Test” : Push button for 
“Reward” : Push button for manually delivering reward through the GUI
“Run BCI” : Push button for starting run of BCI recording session 
Optogenetics Panel
PMT Gating Time : 
Led ON time :
Frames light ON : 
Selection (Opto_Only, Opto_and_Rew, Reward_Only) :
Drifting Gratings Panel
	Selection and seconds (Grey, Black, Randomize drift) : 
Screening (Initial or Final) :
Angles (8 or 12) : Selection of number of angles to be presented 
Trial Type (Standard or Opto) :  
Spatial freq. : 
Temporal freq. :  
Frames (z) : Number of frames that will be collected for the baseline recording
Wait for Trigger checkbox : 
“Rotating bars” : Push button for  
“Retinotopy” : Push button for 
“Drifting Gratings” : Push button for 
“Plot retinotopy” : Push button for 
ROI Panel 
1st frame : 
# frames : Number of frames streamed for motion correction / ROI drawing`
Plane # : Plane to stream frames from if doing dual-plane imaging 
Plane tot. :  
Selection (Mean, Max, Std, All Frames) : How to generate image for motion correction / ROI
“Load Image” : Push button for streaming image for motion correction / ROI drawing
“Draw ROI” : Push button for manually drawing single ROIs
“Load ROI” : Push button for loading saved ROIs from previous sessions
“Save ROI” : Push button for saving drawn or adjusted ROIs
Alignment Only checkbox : Opt to use CLOOSE GUI for FoV/ROI alignment 
Rois :


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
