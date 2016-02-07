%BRANCH2_CLEANED This script contains a modularized version of the analysis
%        included in the script branch2d.m, that process the HC-5 database.
%
%        DESCRIPTION: This script carried out most of the analysis in the files
%        branch2.m using functions. See branch2.m for further details.
%Version 1.0 Ruben Pinzon@2015


clc, close all; %clear all;
run('gpfa4runSettings.m');

basepath        = '~/hc-5/';
[files, animals, roots] = get_matFiles(basepath);


%========================Paramteres and variables==========================
data            = load(files{settings.animal});
mazesect        = data.Laps.MazeSection;
events          = data.Par.MazeSectEnterLeft;
Fs              = data.Par.SamplingFrequency;
X               = data.Track.X;
Y               = data.Track.Y;
eeg             = data.Track.eeg;
time            = linspace(0, length(eeg)/Fs,length(eeg));
speed           = data.Track.speed;
wh_speed        = data.Laps.WhlSpeedCW;
isIntern        = data.Clu.isIntern;
numLaps         = length(events);
spk_clust       = get_spikes(data.Spike.totclu, data.Spike.res);
n_cells         = length(isIntern);
n_pyrs          = sum(isIntern==0);
TrialType       = data.Par.TrialType;
BehavType       = data.Par.BehavType+1;
clear data
% String descriptions
Typetrial_tx    = {'left', 'right', 'errorLeft', 'errorRight'};
Typebehav_tx    = {'first', 'regular', 'uncertain'};
Typeside_tx     = {'left', 'right', 'right', 'left'};
% GPFA training
showpred        = false; %show predicted firing rate
train_split     = true; %train GPFA on left/right separately?
name_save_file  = 'trainedGPFA';
test_lap        = 10;

%%
% ========================================================================%
%==============   (1) Extract trials              ========================%
%=========================================================================%

D = extract_laps(Fs, spk_clust, ~isIntern, X, Y, speed, wh_speed, events, ...
                 struct('in','all','out','all'), TrialType, BehavType);

%show one lap for debug purposes 
if settings.debug
    figure(test_lap)
    raster(D(test_lap).spikes), hold on
    plot(90.*D(test_lap).speed./max(D(test_lap).speed),'k')
    plot(90.*D(test_lap).wh_speed./max(D(test_lap).wh_speed),'r')
end

% ========================================================================%
%==============  (2)  Extract Running Sections    ========================%
%=========================================================================%

%S = get_section(D, in, out, debug, namevar); %lap#1: sensor errors
S = extract_laps(Fs, spk_clust, ~isIntern, X, Y, speed, wh_speed, events, ...
                settings.section, TrialType, BehavType);

% ========================================================================%
%============== (3) Segment the spike vectors     ========================%
%=========================================================================%
%load run model and keep the same neurons
% run = load([roots{animal} '_branch2_results40ms.mat']);

[R,keep_neurons]       = segment(S, 'spike_train', settings.bin_size, Fs, ...
                                 1, settings.min_firing, settings.maxTime);
%%
% ========================================================================%
%============== (4)         Train GPFA            ========================%
%=========================================================================%
laps_all               = select_laps(BehavType, TrialType, 'run');
M                      = trainGPFA(R, laps_all, settings.zDim, showpred, ...
                                   settings.n_folds);

if train_split
    %[R_left, R_right]  = split_trails(R);
    laps_left          = select_laps(BehavType, TrialType, 'run', 1);
    laps_right         = select_laps(BehavType, TrialType, 'run', 2);
    
    if settings.filterTrails
        R_left         = filter_laps(R_left);
        R_right        = filter_laps(R_right,'bins');
    end

    M_left             = trainGPFA(R, laps_left, settings.zDim, showpred, ...
                                   settings.n_folds);
    M_right            = trainGPFA(R, laps_right, settings.zDim, showpred, ...
                                   settings.n_folds);
end

%%
% ========================================================================%
%============== (5)    Save data                  ========================%
%=========================================================================%
fprintf('Will save at %s\n',[roots{settings.animal} name_save_file])
save([roots{settings.animal} name_save_file '_' settings.namevar '.mat'], ...
     'M', 'M_left', 'M_right', 'R', 'keep_neurons')

%%
% ========================================================================%
%============== (6)    Show Neural Trajectories   ========================%
%=========================================================================%

colors = cgergo.cExpon([2 3 1], :);
labels = [R.type];
Xorth = show_latent({M},R,colors, labels);

%%
%=========================================================================%
%=========(7) Compare mean spike counts              =====================%
%=========================================================================%
figure(7)
set(gcf,'position',[100 100 500*1.62 500],'color','w')
plot(mean([R_left.y],2),'r','displayname','wheel after left')
hold on
plot(mean([R_right.y],2),'b','displayname','wheel after right')
ylabel('Average firing rate')
xlabel('Cell No.')
set(gca,'fontsize',14)
savefig()
%=========================================================================%
%=========(8) Compute loglike P(run|model_run)       =====================%
%=========================================================================%

load([roots{animal} name_save_file])
R           = shufftime(R);
%Classification stats of P(run events|model) 
models      = {M_left, M_right};
Xtats       = classGPFA(R, models);
cm          = [Xtats.conf_matrix];
fprintf('hitA: %2.2f%%, hitB: %2.2f%%\n', 100*cm(1,1),100*cm(2,2))

%show likelihood given the models
% plot show likelihood given the models
label.title = 'P(run_j | Models_{left run, right run})';
label.modelA = 'Left alt.';
label.modelB = 'Right alt.';
label.xaxis = 'j';
label.yaxis = 'P(run_j| Models_{left run, right run})';
compareLogLike(R, Xtats, label)

%XY plot
cgergo = load('colors');

label.title = 'LDA classifier';
label.xaxis = 'P(run_j|Model_{left run})';
label.yaxis = 'P(run_j|Model_{right run})';
LDAclass(Xtats, label, cgergo.cExpon([2 3], :))

%=========================================================================%
%=========(9) Compute loglike P(wheel|run_model)     =====================%
%=========================================================================%
%#TODO: Separate this part v in a different script

in              = 'wheel';
out             = 'wheel';
maxTime         = 6;
allTrials       = true; %use all trials of running to test since they are 
                        %all unseen to the wheel model

S = get_section(D, in, out, debug, namevar); %lap#1: sensor errors 
W = segment(S, bin_size, Fs, keep_neurons,...
                [namevar '_spike_train'], maxTime);
W = filter_laps(W);
W = W(randperm(length(W))); 

models      = {M_left, M_right};
Xtats       = classGPFA(W, models,[],allTrials);
cm          = [Xtats.conf_matrix];
fprintf('hitA: %2.2f%%, hitB: %2.2f%%\n', 100*cm(1,1),100*cm(2,2))

% plot show likelihood given the models
label.title = 'P(wheel_j | run model)';
label.modelA = 'Run rigth alt.';
label.modelB = 'Run left alt.';
label.xaxis = 'j';
label.yaxis = 'P(wheel_j|run model)';
compareLogLike(R, Xtats, label)

%XY plot
label.title = 'Class. with Fisher Disc.';
label.xaxis = 'P(wheel_j|run right)';
label.yaxis = 'P(wheel_j|run left)';
LDAclass(Xtats, label)     