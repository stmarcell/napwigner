%% Analysis of Buzsaki database i01_maze06_MS.002


clc, close all; clear all;

basepath        = '/media/bigdata/i01_maze05.005/';
animal          = 'i01_maze05_MS.005';
basepath        = '/media/bigdata/i01_maze06.002/';
animal          = 'i01_maze06_MS.002';
basepath        = '/media/bigdata/i01_maze06.005/';
animal          = 'i01_maze06_MS.005';
basepath        = '/media/bigdata/i01_maze08.001/';
animal          = 'i01_maze08_MS.001';
%basepath        = '/media/bigdata/i01_maze08.004/';
%animal          = 'i01_maze08_MS.004';
% basepath        = '/media/bigdata/i01_maze13.003/';
% animal          = 'i01_maze13_MS.003';
% basepath        = '/media/bigdata/i01_maze15.002/';
% animal          = 'i01_maze15_MS.002';

obj             = load([basepath animal '_BehavElectrData.mat']);
clusters        = obj.Spike.totclu;
laps            = obj.Laps.StartLaps(obj.Laps.StartLaps~=0); %@1250 Hz
%Adding the end of the last lap because Laps only contains the start
laps(end+1)     = obj.Par.SyncOff;
mazesect        = obj.Laps.MazeSection;
wheelspeed      = obj.Laps.WhlSpeedCW;
XT              = obj.Track.X;
YT              = obj.Track.Y;
X               = obj.Spike.X;
Y               = obj.Spike.Y;
events          = obj.Par.MazeSectEnterLeft;
Fs              = obj.Par.SamplingFrequency;
eeg             = obj.Track.eeg;
speed           = obj.Track.speed;
isIntern        = obj.Clu.isIntern;
numLaps         = numel(laps)-1;
time_seg        = linspace(0, length(eeg)/Fs, length(eeg));
[spk, spk_lap, X, Y, X_lap, Y_lap]  = get_spikes(clusters, obj.Spike.res, laps, X, Y);

N               = size(spk_lap,2);
typetrial       = {'left', 'right', 'errorLeft', 'errorRight'};
trialcolor      = hsv(5);
% Extract spks when the mouse if running and in the wheel to calculate
% neural trajectories. 

for lap = 1:numLaps  
    %(a) Runing in the wheel. Detected based on the speed of the wheel that
    %is a better indicator than the EnterSection time stamp
    idx_lap = laps(lap):laps(lap+1);

    wheelNonZero = find(wheelspeed(idx_lap)~=0) +  laps(lap);
    if ~isempty(wheelNonZero)
        length_wheel(lap) = numel(wheelNonZero)/Fs;
        %last lap does not have wheel
        for neu=1:N
            idx = spk_lap{lap,neu}>=wheelNonZero(1)...
                          & spk_lap{lap,neu}<=wheelNonZero(end);
            SpkWheel_lap{lap,neu} = spk_lap{lap,neu}(idx);
        end
    end
    
    %(b) Runing in the maze. Extracted based on the
    %EnterSection time stamp without considering left-right
    idx_run = [events{lap}(2,1), sum(events{lap}(7:8,1))];
    int_at_maze(lap, :) = idx_run;    
    length_run(lap) = (idx_run(2)-idx_run(1))/Fs;
    %sect 1:enter, 6:exit
    for neu=1:N
        idx = spk_lap{lap,neu}>=idx_run(1) & spk_lap{lap,neu}<=idx_run(end);
        SpkRun_lap{lap,neu} = spk_lap{lap,neu}(idx) - idx_run(1);
        X_Run_Lap{lap,neu}  = X_lap{lap, neu}(idx);
        Y_Run_Lap{lap,neu}  = Y_lap{lap, neu}(idx);
    end
    %Type of trial
    trial{lap}          = typetrial{obj.Laps.TrialType(laps(lap))};
    color(lap,:)        = trialcolor(obj.Laps.TrialType(laps(lap)),:);
end
%%
%Only pyramidal cells
disp('Starting spatial segmentation')
%Segment base on spatial coordinates rather than time.
%interpolate the position to the longest time
[leftT, rightT, failed_trial] = protoLap(XT, YT, length_run, trial, X_Run_Lap, ...
                                         Y_Run_Lap, int_at_maze, Fs, animal, isIntern);

%script to extract the grids
segments     = 40;
roiDims      = [20 150]; %width and length of ROI
connectgrids = 1;
ctrNeuron    = 10; % neuron to plot just see things are going OK
show         = 0;
verbose      = false;

%count spikes inside grids and get normalized firing rate
[rate, time_per_bin, centers] = normFiringRate(XT, YT, X_Run_Lap, Y_Run_Lap, int_at_maze,...
                      segments, show, connectgrids, roiDims, leftT,...
                      rightT, ctrNeuron, trial, verbose);

X_pyr = Fs*rate(~isIntern,:,:);
%% Get data in the format
% Command based GPFA based on DataHigh Library
%
bin_width = 30 ; %30mm
%remove failed trails
laps_success = 1:numLaps;
laps_success(failed_trial) = [];
for ilap = 1 : numel(laps_success)
    real_lap = laps_success(ilap);
    D(ilap).data = X_pyr(:,:,real_lap);
    D(ilap).condition = trial{real_lap}; 
    D(ilap).epochColors = color(real_lap, :);
    D(ilap).trialId = real_lap;
    D(ilap).T = size(D(1).data,2);
    D(ilap).centers = centers(:,:,real_lap);
end

firing_thr      = 0.5; % Minimum norm. firing rate which 
                        % neurons should be kept
m               = mean([D.data],2);
keep_neurons    = m >= firing_thr;
fprintf('%d neurons remained with firing rate above %2.2f Hz\n',...
            sum(keep_neurons),firing_thr)
% Remove low firing rate neurons
for itrial = 1:length(D)
    D(itrial).y = D(itrial).data(keep_neurons,:);
end
%Prepare data for the datahigh library


%% DataHigh(D,'DimReduce')
cells           = sum(keep_neurons);
mask            = false(1,length(D));
yDim            = size(D(1).data, 1);
useSqrt         = 1; % square root tranform?    
show_cv         = true;
zDim            = 10;

%prellocating variables
test_trials     = 1:4; % one left and one right, 3 folds
folds           = size(test_trials,1);



test_mask              = mask;
test_mask(test_trials) = true;
train_mask             = ~test_mask;

train_data      = D(train_mask);
test_data       = D(test_mask);
%training of the GPFA with already binned data (1 ms)
[params, gpfa_traj, ll_tr] = gpfa_mod(train_data,zDim,...
                                         'bin_width', 1);

%Posterior of test data given the trained model
[traj, ll_te]   = exactInferenceWithLL(test_data, params,'getLL',1);
% orthogonalize the trajectories
[Xorth, Corth]  = orthogonalize([traj.xsm], params.C);
traj            = segmentByTrial(traj, Xorth, 'data');
traj            = rmfield(traj, {'Vsm', 'VsmGP', 'xsm'});

%Validation with LNO
cv_gpfa = cosmoother_gpfa_viaOrth_fast(test_data,params,zDim);
cv_gpfa_cell    = struct2cell(cv_gpfa);
y_est           = cell2mat(cv_gpfa_cell(8,:));
y_real          = [test_data.y];
xx              = (y_est - y_real).^2;
mse_fold        = sum(sum(xx));

if show_cv
   figure( 33 )
   image(xx)
   title('Square error (y_{true} - y_{est})^2')
   xlabel('Concatenated laps')
   ylabel('Neuron Num.')
   %show one good and one bad prediction
   [~, rse]        = sort(sum(xx,2));


   for j = 1 : floor(cells/12)
       figure(30+j)
       axis tight
       for i = 1 : 12
           subplot(4,3,i)
           plot(y_est(rse(i+(j-1)*12),:),'r'),hold on 
           plot(y_real(rse(i+(j-1)*12),:),'b')
           ylim([min(y_est(rse(i+(j-1)*12),:)) max(y_real(rse(i+(j-1)*12),:))])
           plot(repmat((1 : length(test_data) - 1) * traj(1).T, 2, 1), ylim, 'k')
           xlim([0, traj(1).T*length(test_data)])
           title(sprintf('Cell %d', rse(i+(j-1)*12)))
       end
   end
end

%% Plot the LV (gpfa_traj.data) and X-Y position to compute place fields in
%latent space
meantraj = zeros(1,2*D(1).T);
for idx_hdv = 1 : zDim
    for ilap = 1:length(D) 
        
       NoLap = D(ilap).trialId;   
       traj  = exactInferenceWithLL(D(ilap), params,'getLL',0); 
        
       x       = XT(int_at_maze(NoLap,1):int_at_maze(NoLap,2));
       y       = YT(int_at_maze(NoLap,1):int_at_maze(NoLap,2));
       spe     = speed(int_at_maze(NoLap,1):int_at_maze(NoLap,2));
       spe     = spe(1:50:end);
       pos     = [x(1:50:end)./1200 y(1:50:end)./1000];
        
       spe_int = interp1(spe, linspace(0, length(spe), 40), 'spline');
       
       
       xDim    = length(pos);
       pDim    = size(traj.xsm, 1);
       hidvar  = zeros(xDim, pDim);
       for i = 1 : pDim
            hidvar(:, i)  = interp1(traj.xsm(i,:), linspace(0,39,xDim), 'spline');
       end
       feat              = [pos, spe, hidvar];
       corX(:,:, ilap )  = corr(feat);
        
       figure(idx_hdv)
       title(sprintf('Hidden variable %d, animal %s',idx_hdv, animal))
       hold on, grid on, box on
       xlabel('Space (mm)')
       ylabel('Space (mm)')
       zlabel('Amplitude') 
       campos([4875.4   -5094.3    22.1])
       plot(D(ilap).centers(1, :), D(ilap).centers(2, :), 'color', [0.4 0.4 0.4])
       plot3(D(ilap).centers(1, :), D(ilap).centers(2, :), traj.xsm(idx_hdv,:),...
             'color',D(ilap).epochColors)
       plot3(D(ilap).centers(1, :), D(ilap).centers(2, :), spe_int./max(spe_int),...
             'color',[0.3 0.8 0.3])
       if strcmp(trial{D(ilap).trialId}, 'left')  
          meantraj(1:D(1).T) =  traj.xsm(3,:) + meantraj(1:D(1).T) ;
       else
          meantraj(D(1).T + 1: end) =  traj.xsm(idx_hdv,:) + meantraj(D(1).T + 1: end);
       end   
    end
end

figure
imagesc(mean(corX,3))
%% Check covarince kernel

Tdif         = repmat((1:D(1).T)', 1, D(1).T) - repmat(1:D(1).T, D(1).T, 1);
figure(idx_hdv)
title(sprintf('GPs covariances, animal %s',animal))

for idx_hdv = 1 : zDim
    subplot(2,5,idx_hdv)
    K            = (1 - params.eps(idx_hdv)) * ...
                    exp(-params.gamma(idx_hdv) / 2 * Tdif.^2) +...
                    params.eps(idx_hdv) * eye(D(1).T);
    imagesc(K)
end