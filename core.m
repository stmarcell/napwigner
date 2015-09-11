%% Analysis of Buzsaki database with GPFA

clc, close all; clear all;

basepath        = '/media/bigdata/';
[files, animals, roots]= get_matFiles(basepath);
% % basepath        = '/media/bigdata/i01_maze15.002/';
% % animal          = 'i01_maze15_MS.002';

for rat = 1: length(animals)
    %=====================================================================%
    %=============          LOAD DB MAT FILE      ========================%
    %=====================================================================%
    
    y_obj             = load(files{rat});
    y_clusters        = y_obj.Spike.totclu;
    y_laps            = y_obj.Laps.StartLaps(y_obj.Laps.StartLaps~=0); 
    y_laps(end+1)     = y_obj.Par.SyncOff;
    y_wheelspeed      = y_obj.Laps.WhlSpeedCW;
    y_XT              = y_obj.Track.X;
    y_YT              = y_obj.Track.Y;
    y_X               = y_obj.Spike.X;
    y_Y               = y_obj.Spike.Y;
    y_events          = y_obj.Par.MazeSectEnterLeft;
    Fs                = y_obj.Par.SamplingFrequency;
    y_speed           = y_obj.Track.speed;
    y_isIntern        = y_obj.Clu.isIntern;
    num_laps           = numel(y_laps)-1;
    [y_spk, y_spk_lap, y_X, y_Y, X_lap, Y_lap]  = get_spikes(y_clusters, y_obj.Spike.res, y_laps, y_X, y_Y);

    num_cells         = size(y_spk_lap,2);
    typetrial         = {'left', 'right', 'errorLeft', 'errorRight'};
    trialcolor        = hsv(5);
    %=====================================================================%
    %====       Extract spks when the mouse is running   =================%
    %=====================================================================%
    y_int_maze    = zeros(num_laps, 2);
    y_len_run     = zeros(num_laps, 1);
    for lap = 1:num_laps  
        %Runing in the maze. Extracted based on the
        %EnterSection time stamp without considering left-right
        idx_run          = [y_events{lap}(2,1), sum(y_events{lap}(7:8,1))];        
        y_int_maze(lap, :) = idx_run;    
        y_len_run(lap) = (idx_run(2)-idx_run(1))/Fs;
        %sect 1:enter, 6:exit
        for neu=1:num_cells
            idx = y_spk_lap{lap,neu}>=idx_run(1) & y_spk_lap{lap,neu}<=idx_run(end);
            y_x_Lap{lap,neu}  = X_lap{lap, neu}(idx);
            y_y_Lap{lap,neu}  = Y_lap{lap, neu}(idx);
        end
        %Type of trial
        y_trial{lap}          = typetrial{y_obj.Laps.TrialType(y_laps(lap))};
        y_color(lap,:)        = trialcolor(y_obj.Laps.TrialType(y_laps(lap)),:);
    end
    clear idx_run idx
    %=====================================================================%
    %===============       Spatial Segmentation          =================%
    %=====================================================================%
    disp('Starting spatial segmentation')
    show         = 1;
    segments     = 60;
    roiDims      = [segments 100]; %[Num rois, height of the roi]
    ctrNeuron    = [95 13]; % neuron and Lap to plot just see things are going OK
    verbose      = false;

    %Segment base on spatial coordinates.
    %Find the typical left and right trajectories
    figure(1)
    set(gcf, 'position', [100, 100, 900, 700])
    y_pr = protoLap(y_XT, y_YT, y_len_run,...
                    y_trial, y_int_maze, Fs, show, roiDims);
    title(sprintf('Spatial segmentation rat %s, bin %3.3f ms',animals{rat},y_pr.bin_size_le))
    saveas(gcf,sprintf('%s_segments_(%d).png',roots{rat},segments))
    close 
    
    %count spikes inside grids and get normalized firing rate
    %(XT, YT, X_Run_Lap, Y_Run_Lap, int_at_maze,gridsR, gridsL, ctrNeuron, trial, verbose)
    [y_rate, time_bin] = normFiringRate(y_XT, y_YT, y_x_Lap, y_y_Lap, y_int_maze,...
                           y_pr.rois_ri, y_pr.rois_le,ctrNeuron, y_trial, verbose);

    X_pyr = Fs*y_rate(~y_isIntern,:,:);
    
    %=====================================================================%
    %===============       Get data in the format DH     =================%
    %=====================================================================%
   
    % Get data in the format
    % Command based GPFA based on DataHigh Library
    %
    %remove failed trails
    num_laps_hit = 1:num_laps;
    num_laps_hit(y_pr.failures) = [];
    D = struct;
    for ilap = 1 : numel(num_laps_hit)
        real_lap = num_laps_hit(ilap);
        D(ilap).data = X_pyr(:,:,real_lap);
        D(ilap).condition = y_trial{real_lap}; 
        D(ilap).epochColors = y_color(real_lap, :);
        D(ilap).trialId = real_lap;
        D(ilap).T = size(D(1).data,2);
        D(ilap).centers_left = y_pr.centers_le;
        D(ilap).centers_right = y_pr.centers_ri;
    end

    firing_thr      = 0.5; % Minimum norm. firing rate which 
                            % neurons should be kept
    m               = mean([D.data],2);
    keep_neurons    = m >= firing_thr;
    fprintf('%d neurons remained with firing rate above %2.2f Hz\n',...
                sum(keep_neurons),firing_thr)
    num_cells_pyr = sum(keep_neurons);
    % Remove low firing rate neurons
    for itrial = 1:length(D)
        D(itrial).y = D(itrial).data(keep_neurons,:);
    end
    %Prepare data for the datahigh library
    time_bin_mu = mean(time_bin(:, num_laps_hit));
    fprintf('average bin time per lap %3.3f\n', time_bin_mu)
    
    figure(2)
    set(gcf,'position',[69,719,1820,311])
    imagesc([D.y]), hold on
    plot(repmat([1:numel(num_laps_hit)]*segments,2,1), ylim, 'w')
    xlabel(sprintf('Segmented %d laps each one in %d bins', numel(num_laps_hit), segments))
    ylabel('Pyramidal Cells (EC and CA1)')
    title(sprintf('Spike counts rat %s, bin %3.3f ms',animals{rat},y_pr.bin_size_le))
    saveas(gcf,sprintf('%s_spikecount_bin_(%d).png',roots{rat},segments))
    close
    
    clear y*
    
end 
%% 
cells           = sum(keep_neurons);
mask            = false(1,length(D));
yDim            = size(D(1).data, 1);
useSqrt         = 1; % square root tranform?    
show_cv         = true;
zDim            = 10;

%prellocating variables
test_trials     = 1:6; % one left and one right, 3 folds
folds           = size(test_trials,1);

test_mask              = mask;
test_mask(test_trials) = true;
train_mask             = ~test_mask;

train_data      = D(train_mask);
test_data       = D(test_mask);
%training of the GPFA with already binned data (1 ms)
[params, gpfa_traj] = gpfa_mod(train_data,zDim,...
                                         'bin_width', 1);

%Posterior of test data given the trained model
[traj, ll]   = exactInferenceWithLL(test_data, params,'getLL',1);
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
        
       x       = y_XT(y_int_maze(NoLap,1):y_int_maze(NoLap,2));
       y       = y_YT(y_int_maze(NoLap,1):y_int_maze(NoLap,2));
       spe     = y_speed(y_int_maze(NoLap,1):y_int_maze(NoLap,2));
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
       if strcmp(y_trial{D(ilap).trialId}, 'left')  
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

%% SPWs
close all, clear *spw
SPWs        = load([basepath, '__p111111-0101_SPWs.txt']);
spw_cnt     = 1;
for ispw = 1:2:length(SPWs)-1 
    figure, hold on
    spw_start    = SPWs( ispw );
    spw_end      = SPWs( ispw + 1 );
    len_spw(spw_cnt)      = (spw_end - spw_start)/1000;
    space        = 0;
    cnt_n        = 1;
    for neu=1:num_cells
        if ~y_isIntern(neu)
            idx = y_spk{neu}>=spw_start & y_spk{neu}<=spw_end;
            Spk_spw{spw_cnt,neu} = y_spk{neu}(idx) - spw_start;

            plot(Spk_spw{spw_cnt,neu},space*ones(1,length(Spk_spw{spw_cnt,neu})),'.')
            space     = space + 10;
            cnt_n     = cnt_n + 1;
        end
    end
    spw_cnt     = spw_cnt + 1;
    
end

%show the SPWs


