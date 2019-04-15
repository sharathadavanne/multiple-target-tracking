% Rao-Blackwellized Monte Carlo Data Association
% Toolbox for Matlab 7.x,
% Version 1.0, February 19. 2008
%
% Copyright (C) 2003-2008 Simo Särkkä, <simo.sarkka@hut.fi>
%               2008      Jouni Hartikainen <jmjharti@lce.hut.fi>
%
% This software is distributed under the GNU General Public
% Licence (version 2 or later); please refer to the file
% Licence.txt, included with the software, for details.
%
%
% RBMCDA for known and constant number of targets:
%   MCDA_INIT           Particle structure initialization 
%   KF_MCDA_PREDICT     KF RBMCDA prediction
%   KF_MCDA_UPDATE      KF RBMCDA update 
%   KF_MCDA_SMOOTH      RTS smooothing for RBMCDA particles
%   EKF_MCDA_PREDICT    EKF RBMCDA prediction
%   EKF_MCDA_UPDATE     EKF RBMCDA update
%   EKF_MCDA_SMOOTH     ERTS smooothing for RBMCDA particles
%   UKF_MCDA_PREDICT    UKF RBMCDA prediction
%   UKF_MCDA_UPDATE     UKF RBMCDAupdate
%   UKF_MCDA_SMOOTH     URTS smooothing for RBMCDA particles
%
%
% RBMCDA for unknown and time-varying number of targets:
%   NMCDA_INIT          Particle structure initialization
%   KF_NMCDA_PREDICT    KF RBMCDA prediction
%   KF_NMCDA_UPDATE     KF RBMCDA update 
%   KF_NMCDA_COLLECT    Collects and smooths the target trajectories
%   EKF_NMCDA_PREDICT   EKF RBMCDA prediction
%   EKF_NMCDA_UPDATE    EKF RBMCDA update
%   EKF_NMCDA_COLLECT   Collects and smooths the target trajectories
%   UKF_NMCDA_PREDICT   UKF RBMCDA prediction
%   UKF_NMCDA_UPDATE    UKF RBMCDAupdate
%   UKF_NMCDA_COLLECT   Collects and smooths the target trajectories
%
%
% Functions for handling particles
%   CALCULATE_MEAN      Calculation the mean estimates of particles
%   EFF_PARTICLES       Calculates the effective number of particles
%   GET_WEIGHTS         Extraction of weights from particle structures 
%   NORMALIZE_WEIGHTS   Weights normalization for particle structures    
%   RESAMPLE            Resampling for particle filtering
%   SET_WEIGHTS         Setting of particle weights
%
%
% Misc:
%   BIN_PDF             PDF of Binomial distribution
%   CATEG_RND           Draw random category
%   GAMMA_CDF           CDF of Gamma distribution
%   GAMMA_PDF           PDF of Gamma distribution
%   GMM_PDF             PDF of Gaussian mixture
%   GMM_RND             Random samples of Gaussian mixture distribution
%   POISSON_RND         Random samples of Poisson distribution
%   
%
% /DEMOS/
%   /CLUTTER_DEMO/             
%      CLUTTER_DEMO       Single target tracking with clutter measurements
%      CLUTTER_PLOT1      Plots the filtered estimates using particles
%      CLUTTER_PLOT2      Plots the smoothed estimates using particles
%
%   /BOT_TT_DEMO/          
%      AZ_H               Measurement model function
%      AZ_DH_DX           1st order derivative of the measurement model 
%      EKF_BOT_TT_DEMO    Bearings only tracking of two targets demo with EKF/RBMCDA
%      UKF_BOT_TT_DEMO    Bearings only tracking of two targets demo with UKF/RBMCDA
%
%   /MULTISIGNAL_DEMO/
%      MULTISIGNAL_DEMO     Estimation of multiple 1D Gaussian signals demo 
%      MULTISIGNAL_DEMO_DP  Same as above with the death processing in the prediction step
%      MULTISIGNAL_PLOT     Plots filtered and smoothed signals
% 
%   /MT_DEMO/
%      KF_MT_DEMO         Estimation of multiple 2D Gaussian targets with linear measurements
%      KF_MT_DEMO_DP      Same as above with the death processing in the prediction step
%      KF_MT_DEMO2        Same as KF_MT_DEMO with the restriction of data associations
%      KF_MT_DEMO2_DP     Same as above with the death processing in the prediction step     
%     
%  Currently not documented demos:
%
%   /EKF_MT_DEMO
%      AZ_H               Measurement model function
%      AZ_H_2A            Measurement model function with 2D attributes
%      AZ_DH_DX           1st order derivative of the measurement model 
%      AZ_DH_DX_2A        1st order derivative of the measurement model with 2D attributes
%      EKF_MT_DEMO1       Same as KF_MT_DEMO where the measument model is the bearings only model
%      EKF_MT_DEMO1_DP    Same as above with the death processing in the prediction step
%      EKF_MT_DEMO2       Same as EKF_MT_DEMO1 with additional 2D attribute measurements
%      EKF_MT_DEMO2_DP    Same as above with the death processing in the prediction step
%
%
