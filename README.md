# Tracking unknown number of multiple stationary and moving sources in 2D using particle filter

This repository presents how to use the Rao-Blackwellized particle filtering for tracking unknown number of 2D targets [proposed by Simo Särkkä et. al.](http://becs.aalto.fi/en/research/bayes/rbmcda/mt_demo.html). The original docmumentation for the method can be read [here](http://becs.aalto.fi/en/research/bayes/rbmcda/). Specifically, this script is a modified version of the [original script by Särkkä et. al.](http://becs.aalto.fi/en/research/bayes/rbmcda/html_doc_demos/src/demos/mt_demo/kf_mt_demo_dp.html), adapted for the real example of tracking unknown/multiple number of sound sources in complete 2D space represented using azimuth and elevation angles, also referred as direction of arrival (DOA) estimation.

This work was used as a baseline to compare the performance of a [deep neural network (DNN) for tracking multiple moving sources.](https://github.com/sharathadavanne/seld-net) Check citation information below. To read more about the [general approaches for sound event localization and tracking refer here](https://www.aane.in/research/sound-event-localization-and-tracking).

This script reads the frame-wise 2D target location from a CSV file. Each row of the CSV file consists of the time in seconds, and estimated
location. As the estimated location, we use the DOA of a sound source represented using the azimuth (in 0-360 degree range) and elevation angles (in 0-180 degree range). If more than one source occurs in the same frame, the consecutive rows will contain the spatial location of each sources with an identical time stamp.

The output of the script is the list of DOA tracks in `*.mat` format. Additionally, this script also visualizes the estimated tracks, and its respective reference. 

For custom dataset, you will need to tune the variables in the `PARAMETER` section of the script. I have made brief notes on how these parameters can be tuned based on my understanding through tuning, and from the original documentation. 

The repository provides two example frame-wise sequences - stationary and moving sources - to test the code. The visualization of the tracking result given by the script for the stationary example is shown below.

<p align="center">
   <img src="https://github.com/sharathadavanne/multiple-target-tracking/blob/master/images/stationary_azi.png" width="400" height= "300" title="Visualizing the azimuth part">
   <img src="https://github.com/sharathadavanne/multiple-target-tracking/blob/master/images/stationary_ele.png" width="400" height= "300" title="Visualizing the elevation part">
</p>

Similarly the visualization for the moving sources is shown below.

<p align="center">
   <img src="https://github.com/sharathadavanne/multiple-target-tracking/blob/master/images/moving_azi.png" width="400" height= "300" title="Visualizing the azimuth part">
   <img src="https://github.com/sharathadavanne/multiple-target-tracking/blob/master/images/moving_ele.png" width="400" height= "300" title="Visualizing the elevation part">
</p>

## Keywords
Multiple target tracking, multiple object tracking, acoustic tracking, particle filter, bayesian filter, DOA tracking

## Getting Started
* The Matlab script `multisignal_tracking_2d.m` uses two repositories saved in the `matlab_packages` folder. If needed these can be downloaded directly from the original repository of [RBMCDA Toolbox for Matlab V1.0](http://becs.aalto.fi/en/research/bayes/rbmcda/install.html) and [EKF/UKF Toolbox for Matlab V1.3](http://becs.aalto.fi/en/research/bayes/ekfukf/install.html).
* Choose the example to track (stationary/moving sources) in [line 31](https://github.com/sharathadavanne/multiple-target-tracking/blob/596e1fc962505117649fe62856513eedebaed647/multisignal_tracking_2d.m#L31) of the script and run the script. This should produce the visualization similar to the one seen above. The obtained visualization may not be identical to the one above, since the particle filter uses a random process.
* For custom sequences/dataset, tune the parameters according to the description in the script. Preferably, have a bunch of sequences each for training and a validation split. Tune the parameters on the training split, use the best parameters obtained on the validation split. In a similar fashion to the work in [Sound event localization, detection, and tracking of multiple overlaping and moving sources using convolutional recurrent neural network.](https://github.com/sharathadavanne/seld-net)


## Citation
If you are using this code or the test sequence in any format, then please consider citing the following papers

> Sharath Adavanne, Archontis Politis, Joonas Nikunen, and Tuomas Virtanen, "Sound event localization and detection of overlapping sources using convolutional recurrent neural network" in IEEE Journal of Selected Topics in Signal Processing (JSTSP 2018)

> Sharath Adavanne, Archontis Politis and Tuomas Virtanen, "Localization, detection, and tracking of multiple moving sources using convolutional recurrent neural network" submitted in IEEE Workshop on Applications of Signal Processing to Audio and Acoustics (WASPAA 2019)

## License
Similar to the original repository of Simo Särkkä et. al., this repository is distributed under the [GNU General Public License (version 2 or later)](http://www.gnu.org/copyleft/gpl.html).
