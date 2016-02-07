%BRANCH2_CLEANED This script contains settings for GPFA training and procesing
%
%        DESCRIPTION: This script carried is loaded by gpfa4run
%
% Marcell Stippinger, 2016

%dataset
settings.animal = 3;
%section in the maze to analyze
settings.section.in      = 'mid_arm';
settings.section.out     = 'lat_arm';
settings.debug           = false;
settings.namevar         = 'run';
%segmentation and filtering of silent neurons
settings.bin_size        = 0.04; %duration (s)
settings.min_firing      = 1.0; %minimium firing rate (Hz)
settings.filterTrails    = false; % filter trails with irregular speed/spike count?
% GPFA training
settings.n_folds         = 3; %CV folds
settings.zDim            = 5; %latent dimension
settings.maxTime         = 0; %maximum segmentation time 0 if use all
