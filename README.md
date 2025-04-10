# Welcome to CLOOSE! <img src="Logo/Logo.png" width="500" title="CLOOSE" alt="CLOOSE" align="right" vspace = "100">
# 1. Overview
Pipeline for running BCI experiments and online analysis of 1 and 2 photon imaging data

Copyright (C) 2025  Massachusetts Institute of Technology

CLOOSE is currently under development and not all the features might yet be compatible with one another. If in the meantime you are interested in trying it out for your experiments, get in touch with us at harnett at mit dot edu (Mark Harnett) or valeriof at mit dot edu (Valerio Francioni).  

## 1.1 Paper Abstract

Brain-Computer Interfaces (BCI) have catalyzed advancements in both clinical applications and basic neuroscience research. However, technical barriers such as steep learning curves and complex synchronization requirements often impede their widespread adoption. In response to the increasing demand for optical closed-loop experiments and to the technical challenges associated with their implementation, we introduce CLOOSE (Closed-loop Optical Open-Source Experiments), a versatile platform designed to standardize, simplify, and accelerate closed-loop experiments with functional imaging data. CLOOSE interfaces easily with any frame-based imaging system via TCP/IP, allowing for real-time data streaming and distributed computation. Benchmark tests validate CLOOSE's real-time accuracy in image registration, signal processing, and analysis at imaging frequencies ≥1 kHz, making it the first optical BCI compatible with voltage indicators. Throughout the paper we showcase CLOOSE's versatility in supporting several neurofeedback paradigms: from single neuron to population dynamics, multiplane imaging, and online (and offline) z-tracking. CLOOSE’s functionality is easily accessible to users with minimal coding experience through an intuitive graphical interface for experiment setup and real-time performance monitoring. By significantly lowering the barrier to performing technically demanding in vivo imaging experiments, CLOOSE advances the streamlining, standardization, and democratization of online image analysis and neurofeedback BCI paradigms, opening new avenues for precise manipulation and real-time readout of neuronal activity in health and disease.
 

This code was written by Valerio Francioni.
For support, please open an [issue](https://github.com/harnett/CLOOSE/issues).

## 1.2 Citations

If you use this package in your research, please cite the [paper](BIORXIV):

Valerio Francioni, Anna Beltramini, Linlin Z Fan and Mark T. Harnett (2025). CLOOSE  : An open-source platform to perform optical closed-loop experiments.

The online motion correction algorythm is based on [this paper](https://opg.optica.org/ol/abstract.cfm?uri=ol-33-2-156)

## 1.3 Before you start
Make sure you have the following toolboxes installed with your Matlab 
1. DSP System Toolbox
2. Signal Processing toolbox
3. Image Processing Toolbox
4. If you are planning on running your BCI experiments using visual stimuli as feedback, make sure you have Psychtoolbox installed. Details can be found [here](http://psychtoolbox.org/download).

# 2. Getting started
## Step 1 -  Make sure you can transfer imaging data via a TCP/IP connection internally
1. Open two instances of MATLB on the same machine
2. In the TCP_IP folder open stream_send_data.m in one instance of MATLAB
3. In the TCP_IP folder open stream_receive_data.m in the other instance of MATLAB
4. In the stream_send_data.m script, set the path to Generate_images folder (e.g. YourPath2GitHub\GitHub\CLOOSE\Generate_images)
5. For this step, leave the IPv4 variable unchanged to 127.0.0.1 (Localhost) in both the stream_send_data.m and the stream_receive_data.m.
6. Run stream_send_data.m first, then run stream_receive_data.m
7. CONGRATULATIONS! YOU ARE STREAMING DATA BETWEEN TWO INSTANCES OF MATLAB ON THE SAME PC

## Step 2 -  Make sure you can transfer imaging data via a TCP/IP connection internally
1. Open two instances of MATLB on two different computers, connected to the same Network (ethernet will be faster but Wifi will work too)
2. In the TCP_IP folder open stream_send_data.m on one PC
3. In the TCP_IP folder open stream_receive_data.m in the other PC
4. In the stream_send_data.m script, set the path to Generate_images folder (e.g. YourPath2GitHub\GitHub\CLOOSE\Generate_images)
5. On the PC _sending_ data, set the IPv4 as the IPv4 of the receiving PC (e.g., 10.93.6.184). NB: You can find the IPv4 adress by typing ipconfig/all in your command window.
6. On the PC _receiving_ data, set the IPv4 as the IPv4 of the sending PC (e.g., 10.18.1.121). 
7. Run stream_send_data.m first, then run stream_receive_data.m
8. CONGRATULATIONS! YOU ARE STREAMING DATA BETWEEN TWO DIFFERENT PCs.

## If you managed to run these two steps successfully you are 90% there. 
CLOOSE is setup to act at the receiving end of this process. All you have to do, is to stream data to it. The rest of the analysis can easily be setup using the GUI (next steps). All you have to do now, is to setup your acquisition machine (the PC where you are acquiring imaging data), in the same way as the stream_send_data.m is setup. 
## Critical steps for this are: 
1. Open a TCP/IP connection on your image acquisition device in the same way it's done in the stream_send_data.m script. xpixels is the number of lines in your image. ypixels is the number of columns in your image. Set buffsz to 20.  
```	
 tcpipServer = tcpip(IPv4, port, 'NetworkRole', 'client', ...
    		'InputBufferSize', (xpixels*ypixels*buffsz), 'Terminator', 'CR/LF');
fopen(tcpipServer);
```

2. Identify in your microscope code where the variable encoding for your image is stored. Normally acquisition devices plot your imaging data so you can visualize them online. Look for plotting functions in your code and you should be able to find your image. 

3. Vectorize your image in the same way it's done in the stream_send_data.m script. In this case the vectorized image is the variable 'stack', which is tranformed into a column vector. NB The code below vectorizes a single frame (iframe)
```	
 vect_img = reshape(stack(:, :, iframe), [], 1);
```
4. Open CLOOSE and run your experiments.


# 3.1. Using the GUI: Parameters (Inputs)

The CLOOSE GUI is equipped with various settings to accommodate different experimental designs and uses. The user will set their path within the General Panel, and ensure all settings are correct for their specific experiment. The following parameters exist within the CLOOSE GUI: 

### General Panel
_Trial Number_ : Readout of current trial number during the session

_Path_ : Directory where data will be saved 

_save_ : Will save your session if checked 

_% correct_ : Readout of percentage of trials in which the experimental model reached target

_% avers_ : Readout of percentage of trials that have been incorrect/aversive (Will be removed in future versions)

_Day_ : Current day of training (Will be removed in future versions)

### BCI Panel 
_Motion Corr_ checkbox : Opt to use CLOOSE online motion correction

_Stim Baseline_ checkbox : Opt to have stimuli (visual or auditory) presented during baseline recording

_Activity Levels_ Display: Readout of feedback during session. Reflects how close the model is to reaching target activation

_Trials_ : Number of trials to be acquired during the BCI session

_Frames (z)_ : Lenght of a single trial in terms of number of frames

_target_ : Fraction of trials in which the animal will reach target based on baseline activity. Effectively determines the difficulty of the task 

_I.T.I_ : Inter trial Interval 

_Baseline (z)_ : Lenght of the baseline session, in terms of number of frames

_Fold_: For Scanbox users only. If you know, you know

_Test_ Push button: Currently for debugging and offline plotting in GUI. Might be removed in future releases.

_Baseline_ Push button: Starts running the baseline session

_Run BCI_ Push button: Starts running the closed-loop bci session

### MC Panel 

_Nq_ : Number of quadrants that will be used for motion correction

_NqX_ : x pixels (columns) of each motion correction quadrant

_NqY_ : y pixels (rows) of each motion correction quadrant

_Dsmp_ : Downsampling factor for motion correction. If set to 5 for example, it will only motion correct every 5th frame

### Exp Design Panel

_1 Pop_ : Will setup the baseline recording and the closed-loop part of your BCI experiment, to perfom a 1 population experiment. Uncertain about what that is? Chck out our paper (Link above) 

_2 Pop_ : Will setup the baseline recording and the closed-loop part of your BCI experiment, to perfom a 2 populations experiment. Uncertain about what that is? Chck out our paper (Link above) 

_Pop dynamics_ : Will setup the baseline recording and the closed-loop part of your BCI experiment, to perfom a population dynamics experiment. Uncertain about what that is? Chck out our paper (Link above) 

### Image Panel

_Pixels(x)_ : Number of columns being acquired by image acquisition software 

_Lines (y)_ : Number of lines being acquired by image acquisition software

_Quality_ : Quality of the data being streamed into CLOOSE. Can either be uint8, uint16 or uint32. Do you need a different format? Get in touch

_Dsmp_ : CLOOSE will average and process your data every nth image. 

### Feedback Panel

_Visual_ : CLOOSE will map neuronal activity to a visual feedback stimulus (rotating gabor grating).

_Audio_ : CLOOSE will map neuronal activity to an auditory feedback stimulus (rotating gabor grating).

_Spatial freq_ : Spatial Frequency of the feedback stimulus. 

_Angle_ : Target (Rewarded?) Angle during the closed-loop session.

### Signal Panel

_DF/F0_ : Will use DF/F0 as the main signal to drive the BCI. Reccomended for calcium inidcators and iGluSnfr experiments

_Spikes_ : Will use firing rates as the main signal to drive the BCI. Reccomended for experiments using voltage indicators

_Freq. Cutoff_ : The high-pass filter frequency cutoff for spike detection using voltage indicators

_Nstd_ : Number of standard deviations above noise for spike detection using voltage indicators


### Pop Dynamics Panel
  
_PCA_ : Principal Component Analysis for dimensionality reduction. Only for users selecting Pop Dynamics in the Exp Design panel.

_tSNE_ : t-distributed Stochastic Neighbor Embedding for dimensionality reduction. Only for users selecting Pop Dynamics in the Exp Design panel.

_Ndim_ : In how many dimensions you want your data collapsed. Can be anything between 1 and the number of ROIs for PCA. It will be equal 2 for tSNE.

### Online z-tracking Panel

_z-tracking_ checkbox : Opt for online z-tracking. Plotting will pop up in a separate figure panel

_Load Reference Planes_ : Load a z-stack (1 image per plane) of  .tiffs to be used as a reference volume to compare your incoming images against

_N Planes_ : The number of planes in your z-stack 

_Sq Size_: Determnines the size of the square (Sq x Sq) used for z-tracking

_Load image_ : Instead of loading a previously acquired z-stack, you can use this button to acquire a new z-stack online.

### ROI Activity Panel

Plots the activity traces (either DF/F0 or spikes) online

### ROI Panel

_Alignment Only_ checkbox : opt if you're using CLOOSE to align your FOV in x, y and z with previous days

_Add Filename here_ : Load a FOV from a different recording to use for online motion correction

_# frames_ : Number of frames streamed to CLOOSE for motion correction and ROI drawing

_Plane #_ : Plane(s) to use for BCI if doing multi-plane imaging 

_Plane tot_ : Total number of planes acquired during imaging 

_Selection (Mean, Max, Std, All Frames)_ : How to generate image for motion correction/ROI - We reccomend using mean (other options might be removed in future releases)

_Load Image_ Push button: for streaming image for motion correction / ROI drawing

_Draw ROI_ Push button: For manually drawing ROIs

_Load ROI_ Push button: For loading saved ROIs from previous sessions

_Save ROI_ Push button: For saving drawn or adjusted ROIs

# 3.2. Outputs

### ROIs Folder 

_RoiInfo.ROIpos_ : A 1 x number of ROI Polygon Object, containing the x and y coordinates of your ROIs

_RoiInfo.ROIplane_ : A 1 x number of ROI vector, containing the plane ID (which plane your ROI was drawn on) for your ROIs

_RoiInfo.subpop_ : For each ROI, whether it belongs to P1 or P2 or Pdynamics

_RoiInfo.dir_ : The directory used for loading your ROIs

### Baseline Folder 

In the baselineMask file:

_F0_ : Starting F0 for your ROIs (1 x N rois)

_timeStamp_ : For each loop during the baseline recording, a timestamp (1 x N frames x plane_tot)

_angle_ : For each loop, the angle presented, if you opted in for Stim Baseline (1 x N frames)

_ImTranslateXY_ : estimated x and y motion (2 x N frames)

_subpop1_ : IDs of ROIs belonging to P1 

_subpop2_ : IDs of ROIs belonging to P2

_mean_ : For each ROI, mean activity during baseline recordings (1 x N rois)

_std_ : For each ROI, the activity standard deviation during baseline recordings (1 x N rois)

_roiData_ : For each roi, frame by frame activity during baseline recording (N rois x N frames)

_roiMask_ : Binary mask for each roi

_nroi_ : number of rois

_stimBaseline_ : 1 if Stim Baseline checkbox opted during baseline recording. 0 otherwise

_ROIplane_ : A 1 x number of ROI vector, containing the plane ID (which plane your ROI was drawn on) for your ROIs

_roiFile_ : Directory used to load the ROIs

_zzz_ : 200 resampled trials (size estimated by Frames (z) in the BCI Panel). Used to estimate target

_propTrialsTresh_ : For each activity value, the proportion of trials that would be successful based on baseline recording

_selected_tresh_ : Threshold for target activity

_muFit_ : mean of a fitted Gaussian to activity levels (will be removed)

_sigmaFit_ : std of a fitted Gaussian to activity levels (will be removed)

_selected range_ : minimum activity level reached during baseline. Estimate the range of feedback: From it to target feedback gets split into a user-defined number of bins

_F_all_ : For each roi, frame by frame activity during the 200 trails resampling of baseline recording (N rois x N frames)

_RefImage_ : The image used as reference for motion correction

_plane num_ : Plane IDs with ROIs used for running the BCI

_plane tot_ : Total number of planes recorded (includes plane num + any other plane recorded but not used for BCI)

The Threshold_Values file is redundant and will be removed in future releases

### BCI Folder 

_subpop1_ : ROI IDs for P1

_subpop2_ : ROI IDs for P2

_TStampGlobal_ : For each loop during the BCI recording, a timestamp (1 x N frames x plane_tot)

_Angle_ : Angle presented during visual feedback (1 x N frames). -1 is dark. -5 is grey

_targetAngle_ : target angle

_zscore_ : The 1D activity used to drive the feedback (has nothing to do with a classical z-score. Sorry!)

_bline_loaded_ : Directory to the baseline utilized for the mapping of neural activity and stimulus feedback during the BCI

_ImTranslateXY_ : estimated x and y motion (2 x N frames)

_NewTrial_ : Frame marking the beginning of a new trial

_F0_ : For each roi, for each frame, the F0 value used for estimating DF

_zshift_ : For each frame, the estimated z-shift

_Rotation_ : For each activity level, the rotation of the visual stimulus, if it wasn't binned

_ROIPerFrame_ : For each roi, frame by frame activity during baseline recording (N rois x N frames)

_RefImage_ : The image used as reference for motion correction

_plane num_ : Plane IDs with ROIs used for running the BCI

_plane tot_ : Total number of planes recorded (includes plane num + any other plane recorded but not used for BCI)


# 3.3. Using the GUI: Procedures

### DEMOs and videos on how to (procedurally) use the GUI will be posted here soon. Stay tuned!


# 4. Additional (Misc) Info
## License

**This code is licensed under GPL v3 (no redistribution without credit, and no redistribution in private repos, see the [license](https://docs.github.com/en/repositories/managing-your-repositorys-settings-and-features/customizing-your-repository/licensing-a-repository) for more details).**

## Coming soon
### DEMOs and videos on how to (procedurally) use the GUI will be posted here soon. Stay tuned!

## Logo
Logo was designed by ChatGPT :(:
