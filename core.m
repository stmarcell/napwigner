%% Analysis of Buzsaki database with GPFA

%=========================================================================%
%=================      Preprocessing             ========================%
%=========================================================================%
clc, close all; clear all;
cd /media/LENOVO/HAS/CODE/Wigner-Pattern

basepath        = '/media/bigdata/';
[files, animals, roots]= get_matFiles(basepath);
%needed to compensate for the difference in the measurements among mice
roi_height      = [100, 100, 100, 100, 100, 50, 100];

%create different segmenation sizes
bin_sizes       = 40:5:85;

for rat = 1: 1%length(animals)
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
    [y_spk, y_spk_lap, y_X, y_Y, y_X_lap, y_Y_lap]  = get_spikes(y_clusters, y_obj.Spike.res, y_laps, y_X, y_Y);

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
            y_x_Lap{lap,neu}  = y_X_lap{lap, neu}(idx);
            y_y_Lap{lap,neu}  = y_Y_lap{lap, neu}(idx);
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
    show         = 0;
    verbose      = false;
    ctrNeuron    = [0 5]; % neuron and Lap to plot just see things are going OK
    
    for ibin = 1 : length(bin_sizes)
        fprintf('bin size : %d', bin_sizes(ibin))
        
        segments     = bin_sizes(ibin);
        roiDims      = [segments roi_height(rat)]; %[Num rois, height of the roi]
        

        %Segment base on spatial coordinates.
        %Find the typical left and right trajectories
        if show
            figure(1)
            set(gcf, 'position', [100, 100, 900, 700])
            title(sprintf('Spatial segmentation rat %s, bin %3.3f ms',animals{rat},q_pr.bin_size_le))
            set(gcf, 'PaperUnits', 'centimeters');
            set(gcf, 'PaperPosition', [0 0 18 11]);
        end
        q_pr = protoLap(y_XT, y_YT, y_len_run,...
                        y_trial, y_int_maze, Fs, show, roiDims);
        
        %rat #6 has different spatial scales
        if rat ~= 6 && show
           axis([-100, 1000, 0, 1100])
        end
        if show
            saveas(gcf,sprintf('%s_segments_(%d).png',roots{rat},segments))
        end
        %count spikes inside grids and get normalized firing rate
        %(XT, YT, X_Run_Lap, Y_Run_Lap, int_at_maze,gridsR, gridsL, ctrNeuron, trial, verbose)
        [q_rate, q_time_bin] = normFiringRate(y_XT, y_YT, y_x_Lap, y_y_Lap, y_int_maze,...
                               q_pr.rois_ri, q_pr.rois_le,ctrNeuron, y_trial, verbose);

        X_pyr = Fs*q_rate(~y_isIntern,:,:);

        %=================================================================%
        %===========       Get data in the format DH     =================%
        %=================================================================%

        % Get data in the format
        % Command based GPFA based on DataHigh Library
        %
        %remove failed trails
        num_laps_hit = 1:num_laps;
        num_laps_hit(q_pr.failures) = [];
        D = struct;
        for ilap = 1 : numel(num_laps_hit)
            real_lap = num_laps_hit(ilap);
            D(ilap).data = X_pyr(:,:,real_lap);
            D(ilap).condition = y_trial{real_lap}; 
            D(ilap).epochColors = y_color(real_lap, :);
            D(ilap).trialId = real_lap;
            D(ilap).T = size(D(1).data,2);            
            D(ilap).bin_size = segments;
            D(ilap).timebin = q_time_bin(ilap, :);
            D(ilap).prepro = q_pr;
            
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
            q_spikes    = D(itrial).data(keep_neurons,:);
            D(itrial).y = q_spikes;
            D(itrial).sparsness = sum(sum(q_spikes == 0))/numel(q_spikes);
        end
        %Prepare data for the datahigh library
        %q_time_bin_mu = mean(q_time_bin(:, num_laps_hit));
        %fprintf('average bin time per lap %3.3f\n', q_time_bin_mu)
        if show
            figure(2)
            set(gcf,'position',[69,719,1820,311])
            imagesc([D.y]), hold on
            plot(repmat([1:numel(num_laps_hit)]*segments,2,1), ylim, 'w')
            xlabel(sprintf('Segmented %d laps each one in %d bins', numel(num_laps_hit), segments))
            ylabel('Pyramidal Cells (EC and CA1)')
            title(sprintf('Spike counts rat %s, bin %3.3f ms',animals{rat},q_pr.bin_size_le))
            set(gcf, 'PaperUnits', 'centimeters');
            set(gcf, 'PaperPosition', [0 0 25 10]);
            saveas(gcf,sprintf('%s_spikecount_bin_(%d).png',roots{rat},segments))
        end
        
        save(sprintf('%s_spikecount_bin_(%d).mat',roots{rat},segments), 'D');
         
        close all
        clear q*
    end    
    
    clear y*
    
end 


%% Run PPCA and FA and GPFA, EID

clc, clear all;
cd /media/LENOVO/HAS/CODE/Wigner-Pattern

basepath                = '/media/bigdata/';
[files, animals, roots] = get_matFiles(basepath, '/*count*.mat');
useSqrt                 = 1; % square root tranform?    
zDim                    = [2:2:30]; %hidden space dimensions
rat_curr                = '';
num_dims                = numel(zDim);
colors                  = jet(11);
iter                    = 1;
for irat = 11 : 20 %:length(files)
   
   y_obj            = load(files{irat});
   D                = y_obj.D; 
   num_laps         = length(D);
   num_cells        = size(D(1).y, 1);
   fprintf('Processing animal %s\n', animals{irat})
   
   %get alternations to construct balanced folds (Righ/left)
   cv_alternation = zeros(1, num_laps);
   for ilap = 1 : num_laps
       if strcmp(D(ilap).condition, 'left') || strcmp(D(ilap).condition, 'errorRight')
            cv_alternation(ilap) = 1;      
       end       
   end
   %Paramters of cross validation 3 folds
   %save the folds equal for the same animal, for comparisons
   if ~strcmp(rat_curr, animals{irat})
       cv_mask         = false(1,num_laps);
       cv_trials       = randperm(num_laps);
       cv_fold_idx     = floor(linspace(1,num_laps, 4));  %three chunks
       rat_curr        = animals{irat};
   end
   num_folds = length(cv_fold_idx)- 1;
   
   %======================================================================%
   %================       Training Validation           =================%
   %======================================================================%
   for s = {'fa', 'ppca'}
      eval(sprintf('mse_%s=zeros(1, num_dims);',s{1}))
      eval(sprintf('like_%s= zeros(1, num_dims);',s{1})); 
   end
   parfor_progress(num_dims*num_folds);

   for idim = 1 : num_dims
       tic
        for ifold = 1 : num_folds
            parfor_progress(-1,sprintf('Bin %d,q=%d  ',D(1).T, zDim(idim)));
            
            gp_test_mask = cv_mask;
            gp_test_mask(cv_trials(cv_fold_idx(ifold):...
                             cv_fold_idx(ifold+1)-1)) = true;
            gp_train_mask = ~gp_test_mask; 
            
            for s = {'fa', 'ppca'}
                gp_test_data  = [D(gp_test_mask).y];
                if ~strcmp(s{1},'gpfa')
                    gp_train_data = [D(gp_train_mask).y];
                    
   
                    params   = fastfa(gp_train_data,zDim(idim),'typ',s{1});

                    % compute likelihood on test data FA
                    [chugs, like_fold] = fastfa_estep(gp_test_data,params); 
                    cvdata        = cosmoother_fa(gp_test_data,params);
                    mse_fold      = sum(sum((cvdata-gp_test_data).^2));
                    
                else
                    %training of the GPFA with already binned data (1)
                    [params, gpfa_traj] = gpfa_mod(D(gp_train_mask),zDim(idim),...
                                                         'bin_width', 1);
                    %Posterior of test data given the trained model
                    [gp_traj, gp_ll]   = exactInferenceWithLL(D(gp_test_mask), params,'getLL',1);

                    like_fold = gp_ll;
                    % orthogonalize the trajectories

                    [gp_Xorth, gp_Corth]     = orthogonalize([gp_traj.xsm], params.C);
                    gp_traj                  = segmentByTrial(gp_traj, gp_Xorth, 'data');

                    %Validation with LNO
                    cv_gpfa = cosmoother_gpfa_viaOrth_fast(D(gp_test_mask),params,zDim(idim));
                    cv_gpfa_cell       = struct2cell(cv_gpfa);
                    cvdata             = cell2mat(cv_gpfa_cell(end,:));
                    
                    mse_fold           = sum(sum((cvdata - gp_test_data).^2));
                    
                    
                end
                if idim == 1 && ifold == 1
                    figure(200+iter)
                    stairs(gp_test_data(10,:),'color',[0.4 0.4 0.4]), hold on
                    plot(repmat([1:sum(gp_test_mask)]*D(1).T,2,1),ylim,'color',[0.8 0.8 0.8])
                    axis([0 sum(gp_test_mask)*D(1).T -5 30])
                end
                if ifold == 1 && zDim(idim) == 10 
                    figure(200+iter)
                    
                    if strcmp(s{1},'fa')
                        plot(cvdata(10,:), 'r')
                    else
                        plot(cvdata(10,:), 'b')
                    end
                    axis([0 sum(gp_test_mask)*D(1).T -5 30])
                    drawnow
                end
                
                % add up the likelihood and LNO errors across folds
                eval(sprintf('mse_%s(idim) = mse_%s(idim) + mse_fold;',s{1},s{1}))
                eval(sprintf('like_%s(idim) = like_%s(idim) + like_fold;',s{1},s{1}));
            end
            
        end
   end
   parfor_progress(0);
   typ = {'fa', 'ppca'};
   for j = 1 : 2
       mse_s = eval(sprintf('mse_%s;', typ{j}));
       like_s = eval(sprintf('like_%s;', typ{j}));
       figure(j)
       set(gcf, 'position', [100, 100, 900, 500],'color', 'w')
       subplot(1,2,1)
       plot(zDim,mse_s - mean(mse_s),'-x', 'color', colors(iter,:)),hold on
       
   
       xlabel('Latent dimension')
       ylabel('LNO MSE (centered)')
       set(gca, 'fontsize', 14)
       subplot(1, 2, 2)
       plot(zDim,like_s - mean(like_s),'-x', 'displayname','fa','color', colors(iter,:)),hold on 
       text(10, -iter*100, sprintf('%d bins %2.2e',D(1).T, mean(mse_s)),...
       'color', colors(iter,:),'fontsize', 14)
       xlabel('Latent dimension')
       ylabel('Posterior loglike (centered)')
       set(gca, 'fontsize', 14)
       
       drawnow
   end
   
   iter = iter + 1;
end

for j = 1 : 2
    figure(j)
    set(gcf,'NextPlot','add'); axes;
    h = title(sprintf('%s',animals{irat}));
    set(gca,'Visible','off');
    set(h,'Visible','on');
    saveas(gcf,sprintf('%s_val%s.png',roots{irat},typ{j}))
end
bin_sizes       = 40:5:85;

for j = 201 : 210
    figure(j)
    title(sprintf('%s',animals{irat}));
    set(gcf, 'color', 'w')
    ylim([0 60])
    saveas(gcf,sprintf('%s_predictionFA_PPCA_bin%d.png',roots{irat},bin_sizes(j-200))) 
end
%% Run the GPFA with 10 dims

%=========================================================================%
%============= GPFA based on HIGHDATA  lib (Byron yu et al.) =============%
%=========================================================================%
clc, close all; clear all;
cd /media/LENOVO/HAS/CODE/Wigner-Pattern

basepath                = '/media/bigdata/';
[files, animals, roots] = get_matFiles(basepath, '/*count*.mat');
useSqrt                 = 1; % square root tranform?    
zDim                    = 10; %hidden space dimensions
rat_curr                = '';
num_dims                = numel(zDim);

parfor_progress(length(files));
disp('Starting...')
for irat = 1 : length(files)
    parfor_progress(-1, animals{irat});
    y_obj            = load(files{irat});
    D                = y_obj.D; 
    num_laps         = length(D);
    num_cells        = size(D(1).y, 1);

    %get alternations to construct balanced folds (Righ/left)
    cv_alternation = zeros(1, num_laps);
    for ilap = 1 : num_laps
       if strcmp(D(ilap).condition, 'left') || strcmp(D(ilap).condition, 'errorRight')
            cv_alternation(ilap) = 1;      
       end       
    end
    %Paramters of cross validation 3 folds
    %save the folds equal for the same animal, for comparisons
    if ~strcmp(rat_curr, animals{irat})
       cv_mask         = false(1,num_laps);
       cv_trials       = randperm(num_laps);
       cv_fold_idx     = floor(linspace(1,num_laps, 4));  %three chunks
       rat_curr        = animals{irat};
    end
    num_folds = length(cv_fold_idx)- 1;
    %======================================================================%
    %================       Training Validation           =================%
    %======================================================================%
    for s = {'fa','gpfa'}
       eval(sprintf('cv_data_%s = [];',s{1}))
       eval(sprintf('model_%s = [];',s{1})) 
    end
    
    ifold = 1;
    gp_test_mask = cv_mask;
    gp_test_mask(cv_trials(cv_fold_idx(ifold):...
                     cv_fold_idx(ifold+1)-1)) = true;
    gp_train_mask = ~gp_test_mask; 
    
    
    for s = {'fa','gpfa'}
        gp_test_data  = [D(gp_test_mask).y];
        
        if ~strcmp(s{1},'gpfa')
            gp_train_data = [D(gp_train_mask).y];


            params   = fastfa(gp_train_data,zDim,'typ',s);

            % compute likelihood on test data FA
            [chugs, like_fold] = fastfa_estep(gp_test_data,params); 
            cvdata        = cosmoother_fa(gp_test_data,params);
            mse_fold      = sum(sum((cvdata-gp_test_data).^2));

        else
            %training of the GPFA with already binned data (1)
            [params, gpfa_traj] = gpfa_mod(D(gp_train_mask),zDim,...
                                                 'bin_width', 1);
            %Posterior of test data given the trained model
            [gp_traj, gp_ll]   = exactInferenceWithLL(D(gp_test_mask), params,'getLL',1);

            like_fold = gp_ll;
            % orthogonalize the trajectories

            [gp_Xorth, gp_Corth]     = orthogonalize([gp_traj.xsm], params.C);
            gp_traj                  = segmentByTrial(gp_traj, gp_Xorth, 'data');

            %Validation with LNO
            cv_gpfa = cosmoother_gpfa_viaOrth_fast(D(gp_test_mask),params,zDim);
            cv_gpfa_cell       = struct2cell(cv_gpfa);
            cvdata             = cell2mat(cv_gpfa_cell(end,:));

            mse_fold           = sum(sum((cvdata - gp_test_data).^2));
        end        
        if mod(irat, 10) == 0
            figure
            stairs(gp_test_data(1,:)); hold on
            plot(cvdata(1,:))
            drawnow
        end
        % add up the likelihood and LNO errors across folds
        eval(sprintf('cv_data_%s = cvdata;',s{1}))
        eval(sprintf('model_%s = params;',s{1}))
    end
   
    
    R.modelFA       = model_fa;   
    R.modelGPFA     = model_gpfa;   
    R.cv_dataFA     = cv_data_fa;
    R.cv_dataGPFA   = cv_data_gpfa;
    R.testdata      = gp_test_data;
    %======================================================================%
    %================       Save Data                     =================%
    %======================================================================%
    
    %update data
    save(sprintf('%s_PredictionDim10_bin{%d}.mat',roots{irat},D(1).T), 'R');

    clear post_ll mse_fold tmp* R D  

end
parfor_progress(0);

%======================================================================%
%======================================================================%
%%
clc, close all; clear all;
cd /media/LENOVO/HAS/CODE/Wigner-Pattern

basepath                = '/media/bigdata/';
[cvfiles, animals]      = get_matFiles(basepath, '/*PredictionDim10*.mat');


for irat = 1 : 10 : length(cvfiles)
    load(cvfiles{irat});
    bins        = strfind(cvfiles{irat}, '{');
    num_bins    = str2num(cvfiles{irat}(bins+1:bins+2));
    
    
    cvdata      = R.cv_dataGPFA;
    test_data   = R.testdata;
    
    num_laps    = length(test_data)/num_bins;

    
    %find the error between real and prediction 
    mse_fold           = sum((cvdata - test_data).^2,2);
    [tmp_a, tmp_b]     = sort(mse_fold);
    
    figure(irat)
    cnt         = 1;
    for ifig = 1 : length(tmp_b)        
        
        subplot(2,5, cnt)
        stairs(test_data(tmp_b(ifig),:)), hold on
        plot(repmat([1 : num_laps - 1]*num_bins, 2, 1),ylim,'color',[0.6 0.6 0.6])
        plot(cvdata(tmp_b(ifig),:),'r')
        axis([0 length(test_data) -10 50])
        cnt = cnt + 1;
        if cnt == 11;
           set(gcf,'NextPlot','add'); axes;
           h = title(sprintf('%s bin(%d)',animals{irat},num_bins));
           set(gca,'Visible','off');
           set(h,'Visible','on');
           cnt = 1;
           figure()
        end
    end
    
end


%% Plot the LV (gpfa_traj.data) and X-Y position to compute place fields in
%latent space

clc, close all; clear all;
cd /media/LENOVO/HAS/CODE/Wigner-Pattern

basepath        = '/media/bigdata/';
[files, animals, roots]= get_matFiles(basepath);
zDim            = 10;
%create different segmenation sizes

animal          = 4;

cvresults       = get_matFiles(basepath, '/*PredictionDim10*.mat');
data            = get_matFiles(basepath, '/*count*.mat');

load(data{animal*10-9})
load(cvresults{animal*10-9}) %trained model paramters

params          = R.modelGPFA;

y_obj           = load(files{animal});
y_clusters      = y_obj.Spike.totclu;
y_laps          = y_obj.Laps.StartLaps(y_obj.Laps.StartLaps~=0); 
y_laps(end+1)   = y_obj.Par.SyncOff;
y_XT            = y_obj.Track.X;
y_YT            = y_obj.Track.Y;
y_X             = y_obj.Spike.X;
y_Y             = y_obj.Spike.Y;
y_events        = y_obj.Par.MazeSectEnterLeft;
y_speed         = y_obj.Track.speed;

meantraj = zeros(1,2*D(1).T);
for idx_hdv = 1 : zDim
    for ilap = 1:length(D) 
        
       NoLap = D(ilap).trialId;   
       gp_traj  = exactInferenceWithLL(D(ilap), params,'getLL',0); 
        
      
       idx_run            = [y_events{NoLap}(2,1), sum(y_events{NoLap}(7:8,1))];        
       
       
       x       = y_XT(idx_run(1):idx_run(2));
       y       = y_YT(idx_run(1):idx_run(2));
       spe     = y_speed(idx_run(1):idx_run(2));
       spe     = spe(1:50:end);
       pos     = [x(1:50:end)./1200 y(1:50:end)./1000];
        
       spe_int = interp1(spe, linspace(0, length(spe), 40), 'spline');
       
       
       xDim    = length(pos);
       pDim    = size(gp_traj.xsm, 1);
       hidvar  = zeros(xDim, pDim);
       for i = 1 : pDim
            hidvar(:, i)  = interp1(gp_traj.xsm(i,:), linspace(0,39,xDim), 'spline');
       end
       feat              = [pos, spe, hidvar];
       corX(:,:, ilap )  = corr(feat);
        
       figure(idx_hdv)
       title(sprintf('X_{%d}, animal %s',idx_hdv, animals{animal}))
       set(gcf, 'position', [100, 100, 700, 500], 'color', 'w')
       hold on, grid on, box on
       xlabel('X (mm)')
       ylabel('Y (mm)')
       zlabel('Magnitude') 
       campos([2985   -1206    32.1])
       plot(x,y, 'color', [0.8 0.8 0.8])
       plot(D(ilap).prepro.centers_le(1, :), D(ilap).prepro.centers_le(2, :), 'color', [0.8 0.8 0.8])
       plot(D(ilap).prepro.centers_ri(1, :), D(ilap).prepro.centers_ri(2, :), 'color', [0.8 0.8 0.8])
       if strcmp(D(ilap).condition, 'left') 
           plot3(D(ilap).prepro.centers_le(1, :), D(ilap).prepro.centers_le(2, :), gp_traj.xsm(idx_hdv,:),...
                 'color',D(ilap).epochColors)       
       else
           plot3(D(ilap).prepro.centers_ri(1, :), D(ilap).prepro.centers_ri(2, :), gp_traj.xsm(idx_hdv,:),...
                 'color',D(ilap).epochColors)  
       end   
    end
    drawnow
    saveas(gcf,sprintf('%s_projXpaceHidden(%d).png',roots{animal},idx_hdv))
end

figure
set(gcf, 'position', [800, 100, 700, 500], 'color', 'w')
title(sprintf('Correlation %s', animals{animal}))
imagesc(mean(corX,3))
xlabel('x, y, speed, latent vars.')
ylabel('x, y, speed, latent vars.')
set(gca, 'fontsize', 14)
colorbar
saveas(gcf,sprintf('%s_corrX.png',roots{animal}))

% Check covarince kernel

Tdif         = repmat((1:D(1).T)', 1, D(1).T) - repmat(1:D(1).T, D(1).T, 1);

figure(30)
set(gcf, 'position', [100, 100, 700, 500], 'color', 'w')
figure(20)
set(gcf, 'position', [800, 100, 700, 500], 'color', 'w')
colors = jet(10);
for idx_hdv = 1 : zDim
    subplot(2,5,idx_hdv)
    K            = (1 - params.eps(idx_hdv)) * ...
                    exp(-params.gamma(idx_hdv) / 2 * Tdif.^2) +...
                    params.eps(idx_hdv) * eye(D(1).T);
    for j = 1:10
        X_sampled = mvnrnd(zeros(1, length(K)), K);      
        plot(X_sampled, 'color', colors(j,:)), hold on
    end
    figure(30)
    subplot(2,5,idx_hdv)
    imagesc(K)
    figure(20)
end


set(gcf,'NextPlot','add'); axes;
h = title(sprintf('GPs Covariance %s',animals{animal}));
set(gca,'Visible','off');
set(h,'Visible','on');
saveas(gcf,sprintf('%s_GPsK.png',roots{animal}))

figure(30)
set(gcf,'NextPlot','add'); axes;
h = title(sprintf('Sampled variables %s',animals{animal}));
set(gca,'Visible','off');
set(h,'Visible','on');
saveas(gcf,sprintf('%s_SampledX.png',roots{animal}))
