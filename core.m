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
basepath        = '/media/bigdata/i01_maze08.004/';
animal          = 'i01_maze08_MS.004';
basepath        = '/media/bigdata/i01_maze13.003/';
animal          = 'i01_maze13_MS.003';
basepath        = '/media/bigdata/i01_maze15.002/';
animal          = 'i01_maze15_MS.002';

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
trialcolor      = hot(5);
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
%Only pyramidal cells

%Segment base on spatial coordinates rather than time.
%interpolate the position to the longest time
[leftT, rightT, failed_trial] = protoLap(XT, YT, length_run, trial, X_Run_Lap, ...
                                         Y_Run_Lap, int_at_maze, Fs, animal, isIntern);

%script to extract the grids
segments     = 40;
roiDims      = [20 200]; %width and length of ROI
connectgrids = 1;
ctrNeuron    = 0; % neuron to plot just see things are going OK
show         = 0;

%count spikes inside grids and get normalized firing rate
[rate, duration] = normFiringRate(XT, YT, X_Run_Lap, Y_Run_Lap, int_at_maze,...
                      segments, show, connectgrids, roiDims, leftT,...
                      rightT, ctrNeuron, trial);

X_pyr = Fs*rate(~isIntern,:,:);
%% Get data in the format
% Command based GPFA based on DataHigh Library
%
bin_width = 30 ; %30mm
%remove failed trails
laps_success = 1:numLaps;
laps_success(failed_trial) = [];
for ilap = 1 : numel(laps_success)
    D(ilap).data = X_pyr(:,:,laps_success(ilap));
    D(ilap).condition = trial{laps_success(ilap)}; 
    D(ilap).epochColors = color(laps_success(ilap), :);
    D(ilap).trialId = laps_success(ilap);
    D(ilap).T = size(D(1).data,2);
end

dims            = 3:20; % Target latent dimensions

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
removeCells     = randperm(cells);
mask            = false(1,length(D));
yDim            = size(D(1).data, 1);
useSqrt         = 1; % square root tranform?    
show_cv         = false;

%prellocating variables
test_trials     = [1:4; 5:8; 12:15]; % one left and one right, 3 folds
lat         = []; % Struct for the latent variables
ll_te       = 0;  % these store the likelihood
ll_tr       = 0;  % these store the likelihood
mse_fold    = 0;  % and mse for one fold
folds       = size(test_trials,1);
ll_train    = zeros(length(dims),folds); % keeps a running total of likelihood
ll_test     = zeros(length(dims),folds); % keeps a running total of likelihood
mse         = zeros(length(dims),folds); % and mse    
paramsGPFA  = cell(length(dims), folds);
orth_traje_tr  = cell(length(dims), folds); %training orth_traje
orth_traje_te  = cell(length(dims), folds); %training orth_traje

for idim = 1 : length(dims)
    for ifold = 1 : folds  % three-fold cross-validation        
        % prepare masks:
        % test_mask isolates a single fold, train_mask takes the rest
        test_mask = mask;
        test_mask(test_trials(ifold, :)) = true;

        train_mask = ~test_mask;

        train_data = D(train_mask);
        test_data = D(test_mask);
        %training of the GPFA
        [params, gpfa_traj, ll_tr] = gpfa_mod(train_data,dims(idim),...
                                                 'bin_width', bin_width);

        %Posterior of test data given the trained model
        [traj, ll_te] = exactInferenceWithLL(test_data, params,'getLL',1);
        % orthogonalize the trajectories
        [Xorth, Corth] = orthogonalize([traj.xsm], params.C);
        traj = segmentByTrial(traj, Xorth, 'data');
        traj = rmfield(traj, {'Vsm', 'VsmGP', 'xsm'});

        %Validation with LNO
        cv_gpfa = cosmoother_gpfa_viaOrth_fast...
                                  (test_data,params,1:idim);
        cv_gpfa_cell    = struct2cell(cv_gpfa);
        y_est           = cell2mat(cv_gpfa_cell(6 + idim,:));
        y_real          = [test_data.data];
        xx = (y_est - y_real).^2;
        mse_fold = sum(sum(xx));
        if show_cv
           figure( 33 )
           image(xx)
           title('Square error (y_{true} - y_{est})^2')
           xlabel('Concatenated laps')
           ylabel('Neuron Num.')
           %show one good and one bad prediction
           [~, rse]        = sort(sum(xx,2));
            
           
           for j = 1 : ceil(cells/12)
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
        mse(idim, ifold) = mse_fold;
        ll_test(idim, ifold) = ll_te;
        ll_train(idim, ifold) = ll_tr;
        paramsGPFA{idim, ifold} = params;
        orth_traje_tr{idim, ifold} = gpfa_traj;
        orth_traje_te{idim, ifold} = traj;
        fprintf('Dimension %d, fold %d', idim, ifold)
    end        
end
% best dimension and best across folds
[a, foldmax_te] = max(sum(ll_test));
[a, imax_te] = max(ll_test(:,foldmax_te));
[a, foldmax_tr] = max(sum(ll_train));
[a, imax_tr] = max(ll_test(:,foldmax_tr));

results.ll_test = ll_test;
results.ll_train = ll_train;
results.mse = mse;
results.dim = imax_te;
results.GPFA = paramsGPFA;
results.traj_tr = orth_traje_tr;
results.traj_te = orth_traje_te;

% Setup figure to sumarize results
close all
til = sprintf('Spatial binning');
annotation('textbox', [0 0.9 1 0.1], ...
    'String', til, ...
    'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center',...
    'Fontsize',18)
set(gcf,'color', 'w', 'position', [100 100 1400 700])

% projection over the three best dimensions (SVD)
for itrial = 1:length(train_data)
    p = orth_traje_tr{imax_tr, foldmax_tr}(itrial).data;
    c = orth_traje_tr{imax_tr, foldmax_tr}(itrial).epochColors;
    subplot(2,3,[1 2 4 5]), grid on
    plot(p(1,:), p(2,:), 'Color', c,...
          'linewidth',2); hold on
end
for itrial = 1:length(test_data)
    p = orth_traje_te{imax_te, foldmax_te}(itrial).data;
    c = orth_traje_te{imax_te, foldmax_te}(itrial).epochColors;
    subplot(2,3,[1 2 4 5]), grid on
    plot(p(1,:), p(2,:), 'Color', 0.5*c,...
          'linewidth',2); hold on
end
xlabel('Eigenvector 1')
ylabel('Eigenvector 2')
zlabel('Eigenvector 3')

% MSE
subplot(2,3,3)
mean_mse = mean(mse,2);
plot(mean_mse,'linewidth',2, 'linestyle','-.'), hold on
plot(mse)
plot(imax_te, mean_mse(imax_te),'r*', 'markersize',10)
xlabel('Latent Dimension')
ylabel('Mean Sqaure Error (LNO)')

% LogLike
subplot(2,3,6)
mean_ll_test  = mean(ll_test,2);
mean_ll_train = mean(ll_train,2);
offset_te = mean(mean_ll_test);
offset_tr = mean(mean_ll_train);
plot(dims,mean_ll_test,'linewidth',2, 'linestyle','-.', 'color','k'), hold on
plot(dims,mean_ll_train,'linewidth',2, 'linestyle','-.', 'color','b')
plot(dims,ll_test,'k')
plot(dims,ll_train,'b')
plot(imax_te, mean_ll_test(imax_te),'k*', 'markersize',10)
plot(imax_tr, mean_ll_train(imax_tr),'bs', 'markersize',10)
xlabel('Latent Dimension')
ylabel('Log Likelihood')
title('Train','color','b')
box off
namePNG = sprintf('Results/%s_cells',animal);
print(gcf,[basepath namePNG],'-dpng')

figure(iwin + 1)
numDim = 7; % 9 latent variables
til = sprintf('9 Latent Variables');
annotation('textbox', [0 0.9 1 0.1], ...
    'String', til, ...
    'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center',...
    'Fontsize',18)
set(gcf,'color', 'w', 'position', [100 100 1400 700])
duration = linspace(0, window, T); % from he extraction program
F    = orth_traje_tr{numDim, foldmax_tr};
for l= 1:length(F)
   for v = 1:dims(numDim)
       subplot(3, 3, v)
       plot(duration, F(l).data(v, :),'color',F(l).epochColors), hold on       
   end
end
namePNG = sprintf('Results/%s_LV_w%d',animal, window);
print(gcf,[basepath namePNG],'-dpng')    


save([basepath 'Results/' animal '_run_results.mat'],'results')
