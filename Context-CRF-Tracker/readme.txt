Spatio-temporal Cell tracking MATLAB code (v1.0)

=======================================================================================
Instructions for running the Context-Aware Spatio-temporal Cell tracking MATLAB code
=======================================================================================

While running the code for the first time, from within the 'Context-CRF-Tracker' folder, run

>>setup


The tracker provided generates spatio-temporal correspondences between cells in a 4D spatio-temporal stack of confocal images. The code provided can be used for stacks with images across two observational time points and there is NO limit on the number of z-slices at each of the two observational time points. Note that, although the provided code is written for two time points, it can be trivially extended to generate results on 4D stacks of any size. The code is also capable of handling 'spatial-only'/'temporal-only' collection of confocal images as input.

To generate correspondence/cell tracking results, use the MATLAB function

>>correspondence_mat = cell_tracker(cord, mode, properties_spat, properties_temp);


--------
Inputs:
--------

mode:
------

There are 3 operating modes for this code. 

1. 'spatial' : When the target is to obtain tracking results between cells across multiple z-slices collected at the same observational time point.

2. 'temporal' : When the target is to obtain temporal lineages between confocal image slices collected across multiple observational time points but approximately at the same 'z depth' in the tissue.

3. 'spatio-temporal' : When the target is to obtain correspondences between 2D cell image segments both across space and time.

cord: 
------

For 'spatial' and 'temporal' modes: A MATLAb cell data structure containing 'K' elements where 'k'th element is a structure containing segmented point clouds and other strutural parametes for 2D cell segments at z=k (for spatial) or at t=k (for temporal). A sample cord{k} data structure is provided ('sample_cord_t_z_n.mat').

For 'spatio-temporal' mode: A MATLAB cell data-strcuture containing 2 elements. Each of these 2 members, i.e., cord{1} and cord{2} are the segmented collection of image slices at time points t_1 and t_2 respectively. cord{1} and cord{2} are cell data structures and cord{t}{z}(c_n) is a structure containing segmented point clouds and other strutural parametes for 2D cell segment 'c_n' on the image slice at depth 'z' (e.g. 1, 2, 3 etc.) from the top at 't'th observational time point. A sample of the cord{t}{z} is provided as 'sample_cord_t_z_n.mat'. 

Please note that we have used the modified Watershed algorithm given by Mkrtchyan et al. in ICIP 2011 (“Efficient cell segmentation and tracking of developing plant meristem”, Intl. Conf. on Image Processing, 2011) for generating these segmented cellular data structures. However, any other segmentation algorithm could be used for this purpose. Codes for generating such data structures using the 'level-set' segmentation method can be found at http://www.ee.ucr.edu/~amitrc/celltracking.php

properties_spat / properties_temp:
------------------------------------

A collection of variables directly related to the CRF parameters used for generating the pairwise similarity scores. The structures for the temporal tracking properties (properties_temp) additionally contain the cell division parameters. Explanation on the various parameters can be found in the MedIA 2014 paper. Sample parameter sets are provided. Please note that these parameters are highly data-specific and the values provided in the property structures may need modifications before applying to a different dataset.

--------
Output:
--------

A MATLAB cell data structure containing all the association/correspondence matrices between pairs of image slices. For example, for spatio-temporal cell tracking, correspondence_mat{t_1,t_2}{z,z} contains the association matrix between cells at image slices collected at the depth 'z' in the tissue and at observatioanl time points 't_1' and 't_2'.

-----------------------------------------------------------------------------------------------------------------------
                                                   Updates
-----------------------------------------------------------------------------------------------------------------------

We are currently working on a version 1.1 - a much faster implementation which will be updated at the same place. Keep checking this website for periodic updates on codes.

-----------------------------------------------------------------------------------------------------------------------
                                                   Acknowledgment
-----------------------------------------------------------------------------------------------------------------------

We have used the loopy belief propagation MATLAB code written by Dr. Talya Meltzer and the toolbox is downloaded from 
http://www.cs.huji.ac.il/~talyam/inference.html. We would like to thank Dr. Min Liu and a number of Matlab Central file exchange authors for their useful tools and codes.

-----------------------------------------------------------------------------------------------------------------------
                                                   Copyright
-----------------------------------------------------------------------------------------------------------------------

Copyright(2014) Anirban Chakraborty and Amit K. Roy-Chowdhury, University of California, Riverside.
Email: <anirban.chakraborty@email.ucr.edu> or <amitrc@ece.ucr.edu>.

------------------------------------------------------------------------------------------------------------------------

