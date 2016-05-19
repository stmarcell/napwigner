%GPFA4RUNTEMPLATE This script contains settings for GPFA training and procesing
%
%        DESCRIPTION: This script is loaded by gpfa4run
%
% Marcell Stippinger, 2016

%dataset
settings.animal          = 29;
%section in the maze to analyze
settings.section.in      = 'all';
settings.section.out     = 'all';
settings.debug           = false;
settings.namevar         = 'run';
%segmentation and filtering of silent neurons
settings.bin_size        = 0.05; %duration (s)
settings.min_firing      = 0.1; %minimium firing rate (Hz)
settings.filterTrails    = false; % filter trails with irregular speed/spike count?
% GPFA training
settings.train           = 0; %redo training
settings.n_folds         = 3; %CV folds
settings.zDim            = 3; %latent dimension
settings.maxTime         = 0; %maximum segmentation time, 0 if use all
settings.maxLength       = 100; %maximum unsplitted time series length, 0 if no limit
