%BRANCH2_CLEANED This script contains a modularized version of the analysis
%        included in the script branch2d.m, that process the HC-5 database.
%
%        DESCRIPTION: This script carried out most of the analysis in the files
%        branch2.m using functions. See branch2.m for further details.
%Version 1.0 Ruben Pinzon@2015


clc, close all; %clear all;
run('gpfa4runSectionSettings.m');

basepath        = '~/marcell/_Data_hc-5/';
workpath        = '~/marcell/napwigner/work/';
[files, roots, animals] = get_matFiles(basepath,'.*_BehavElectrData\.mat');


%========================Paramteres and variables==========================
fprintf('\nSelecting %d: %s\n\n',settings.animal,files{settings.animal});
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
trained         = false;


project         = strrep(roots{settings.animal},'/','');
savepath        = [workpath roots{settings.animal}];
fn              = [project '_' name_save_file '_' settings.namevar '.mat'];
if ~exist(savepath,'dir')
    mkdir(savepath);
end

    
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

%lap#1: sensor errors
S = extract_laps(Fs, spk_clust, ~isIntern, X, Y, speed, wh_speed, events, ...
                settings.section, TrialType, BehavType);

% ========================================================================%
%============== (3) Segment the spike vectors     ========================%
%=========================================================================%
%load run model and keep the same neurons
if isCommandWindowOpen() && exist([savepath fn],'file')
    fprintf('Will load from %s\n', [savepath fn]);
    info = load([savepath fn], 'M', 'laps', 'R', 'keep_neurons', 'settings');
    keep_neurons = info.keep_neurons;
    fprintf('Successfully loaded file, you may skip Section (4).\n');
else
    keep_neurons = 1;
end
    
[R,keep_neurons]       = segment(S, 'spike_train', settings.bin_size, Fs, ...
                                 keep_neurons, settings.min_firing, settings.maxTime);

laps.all               = select_laps(BehavType, TrialType, settings.namevar);
laps.left              = select_laps(BehavType, TrialType, settings.namevar, 1);
laps.right             = select_laps(BehavType, TrialType, settings.namevar, 2);


%%
% ========================================================================%
%============== (4)         Train GPFA            ========================%
%=========================================================================%
try
    M.all                  = trainGPFA(R, laps.all, settings.zDim, showpred, ...
                                       settings.n_folds,'max_length',settings.maxLength);

    if train_split
        if settings.filterTrails
            keep_laps      = filter_laps(R(laps.left));
            laps.left      = laps.left(keep_laps);
            keep_laps      = filter_laps(R(laps.right),'bins');
            laps.right     = laps.right(keep_laps);
        end

        M.left             = trainGPFA(R, laps.left, settings.zDim, showpred, ...
                                       settings.n_folds,'max_length',settings.maxLength);
        M.right            = trainGPFA(R, laps.right, settings.zDim, showpred, ...
                                       settings.n_folds,'max_length',settings.maxLength);
    end
    trained = true;
catch ME
    fprintf('Error training GPFA: %s\n', ME.identifier);
    rethrow(ME);
end

%%
% ========================================================================%
%============== (5)    Save / Use saved data      ========================%
%=========================================================================%

if trained
    fprintf('Will save at %s\n', [savepath fn]);
    save([savepath fn], 'M', 'laps', 'R', 'keep_neurons', 'settings');
    trained = false; %#ok<NASGU>
    exit;
else
    M = info.M; %#ok<UNRCH>
end


%%
% ========================================================================%
%============== (6)    Show Neural Trajectories   ========================%
%=========================================================================%

%colors = cgergo.cExpon([2 3 1], :);
colors = hsv(4);
labels = [R.type];
Xorth = show_latent({M.all}, R, colors, labels, Typetrial_tx);

%%
%=========================================================================%
%=========(7) Compare mean spike counts              =====================%
%=========================================================================%
figure(7)
set(gcf,'position',[100 100 500*1.62 500],'color','w')
plot(mean([R(laps.left).y],2),'r','displayname','wheel after left')
hold on
plot(mean([R(laps.right).y],2),'b','displayname','wheel after right')
ylabel('Average firing rate')
xlabel('Cell No.')
set(gca,'fontsize',14)
savefig()

figure(71)
plot_timescales({M.left, M.right, M.all}, colors, {'trained_{left}', 'trained_{right}', 'trained_{all}'})

%%
%=========================================================================%
%=========(8) Compute loglike P(run|model_run)       =====================%
%=========================================================================%

%load([roots{settings.animal} name_save_file])
Rs           = R;
%Rs           = shufftime(R);
%Rs           = shuffspike(R);
[Rs, Malt.left, Malt.right, Malt.all] = normalizedatanmodel(Rs, M.left, M.right, M.all);
%sufficient: 1:4 or 12:16 or 21:22 or 23:24 or 25:30? OR 31:35
%inconclusive: 5:10 and 17:18 and
%spy = [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 20 21 22 23 24 25 26 27 28 29 30 36 37 38 39 40 41 42 43 44 45 46];
%for l = 1:length(Rs)
%    Rs(l).y(spy,:)=0;
%end
%Classification stats of P(run events|model) 
models      = {M.left, M.right};
models      = {Malt.left, Malt.right};
Xtats       = classGPFA(Rs, models);
cm          = [Xtats.conf_matrix];
fprintf('hitA: %2.2f%%, hitB: %2.2f%%\n', 100*cm(1,1),100*cm(2,2))

% %show likelihood given the models
% % plot show likelihood given the models
label.title = 'P(run_j | Models_{left run, right run})';
label.modelA = 'Left alt.';
label.modelB = 'Right alt.';
%label.modelB = 'Global model';
label.xaxis = 'j';
label.yaxis = 'P(run_j| Models_{left run, right run})';
compareLogLike(Rs, Xtats, label)

%XY plot
cgergo = load('colors');

label.title = 'LDA classifier';
label.xaxis = 'P(run_j|Model_{left run})';
label.yaxis = 'P(run_j|Model_{right run})';
LDAclass(Xtats, label, cgergo.cExpon([2 3], :))

%%
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
