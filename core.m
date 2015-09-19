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
%% Run the GPFA

%=========================================================================%
%============= GPFA based on HIGHDATA  lib (Byron yu et al.) =============%
%=========================================================================%
clc, close all; clear all;
cd /media/LENOVO/HAS/CODE/Wigner-Pattern

basepath                = '/media/bigdata/';
[files, animals, roots] = get_matFiles(basepath, '/*count*');
useSqrt                 = 1; % square root tranform?    
zDim                    = [2:2:30]; %hidden space dimensions
rat_curr                = '';
num_dims                = numel(zDim);

for irat = 1 : 20 %length(files)
    
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
   like         = zeros(1, length(cv_fold_idx)-1 );
   mse          = zeros(1, length(cv_fold_idx)-1 );
   model        = {};
   parfor_progress(num_dims*num_folds)
   for idim = 1 : num_dims
       tic
        for ifold = 1 : num_folds         
            parfor_progress(-1,sprintf('Bin %d,q=%d ',D(1).T, zDim(idim)));
            
            gp_test_mask = cv_mask;
            gp_test_mask(cv_trials(cv_fold_idx(ifold):...
                             cv_fold_idx(ifold+1)-1)) = true;
            gp_train_mask = ~gp_test_mask; 
            gp_train_data = D(gp_train_mask);
            gp_test_data = D(gp_test_mask);

            %training of the GPFA with already binned data (1)
            [gp_params, gpfa_traj] = gpfa_mod(gp_train_data,zDim(idim),...
                                                 'bin_width', 1);
            %Posterior of test data given the trained model
            [gp_traj, gp_ll]   = exactInferenceWithLL(gp_test_data, gp_params,'getLL',1);
            
            like_fold = gp_ll;
            % orthogonalize the trajectories

            [gp_Xorth, gp_Corth]     = orthogonalize([gp_traj.xsm], gp_params.C);
            gp_traj                  = segmentByTrial(gp_traj, gp_Xorth, 'data');

            %Validation with LNO
            cv_gpfa = cosmoother_gpfa_viaOrth_fast(gp_test_data,gp_params,zDim(idim));
            cv_gpfa_cell       = struct2cell(cv_gpfa);
            cvdata           = cell2mat(cv_gpfa_cell(end,:));
            gp_y_real          = [gp_test_data.y];
            gp_xx              = (cvdata - gp_y_real).^2;
            
            mse_fold           = sum(sum(gp_xx));
           % add up the likelihood and LNO errors across folds
            mse(idim)          = mse(idim) + mse_fold;
            like(idim)         = like(idim) + like_fold; 
            
            model{idim, ifold}.GPparams = gp_params;         
            
            clear gp* 
        end
        
   end
   parfor_progress(0)
   
   R.gpfa       = model;
   R.resultLL   = like;
   R.resultLNO  = mse;
   
   figure(1)
   set(gcf, 'position', [100, 100, 900, 500])
   subplot(121)
   set(gcf, 'PaperUnits', 'centimeters');
   set(gcf, 'PaperPosition', [0 0 25 11]);   
   plot(zDim, mse - mean(mse), '-x'), hold on, grid on

   ylabel('MSE (centered)')
   xlabel('Latent Dimension GPFA')
   text(10, 1, 'left/right train(test)')
   text(1,1.01,sprintf('%1.2e',mse_fold(1,1)))
   colors = lines(3);
   for ifold = 1 : num_folds
      tmp_test  = model{1,ifold}.test;
      tmp_train = model{1, ifold}.training;
      %number fo left alternations in test set
      tmp_left_test = sum(cv_alternation(tmp_test)==1);
      tmp_righ_test = numel(cv_alternation(tmp_test))-tmp_left_test;      
      tmp_left      = sum(cv_alternation(tmp_train)==1);
      tmp_right     = numel(cv_alternation(tmp_train))-tmp_left;
      text(15, 1-ifold*0.02, sprintf('%d/%d (%d/%d)',tmp_left, tmp_right, tmp_left_test, tmp_righ_test), 'color', colors(ifold,:))
      text(2,0.7-ifold*0.02,sprintf('Lap %d\t ',find(tmp_test==1)), 'color', colors(ifold,:))
   end
   %ylim([0.5, 1.1]) 
   
   subplot(122)  
   plot(zDim, like - mean(like), '-x'), hold on, grid on
   text(1,-1.001,sprintf('%1.2e',abs(post_ll(1,1))))
   ylabel('Posterir LogLike (centered)')
   xlabel('Latent Dimension')   
   
   
   set(gcf,'NextPlot','add'); axes;
   h = title(sprintf('%s bin(%d)',animals{irat},D(1).T));
   set(gca,'Visible','off');
   set(h,'Visible','on'); 
   
   
   saveas(gcf,sprintf('%s_resultGPFA_bin(%d).png',roots{irat},D(1).T))
   close all
   %======================================================================%
   %================       Save Data                     =================%
   %======================================================================%
   
   %update data
   save(sprintf('%s_CV_bin(%d).mat',roots{irat},D(1).T), 'R');
   
   clear post_ll mse_fold tmp* R D  
   
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

for irat = 11 : 10 :length(files)
   
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
   for s = {'gpfa'}
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
            
            for s = {'gpfa'}
                gp_test_data  = [D(gp_test_mask).y];
                if ~strcmp(s{1},'gpfa')
                    gp_train_data = [D(gp_train_mask).y];
                    
   
                    params   = fastfa(gp_train_data,zDim(idim),'typ',s);

                    % compute likelihood on test data FA
                    [chugs, like_fold] = fastfa_estep(gp_test_data,params); 
                    cvdata        = cosmoother_fa(gp_test_data,params);
                    mse_fold      = sum(sum((cvdata-gp_test_data).^2));
                    
                else
                    %training of the GPFA with already binned data (1)
                    [gp_params, gpfa_traj] = gpfa_mod(D(gp_train_mask),zDim(idim),...
                                                         'bin_width', 1);
                    %Posterior of test data given the trained model
                    [gp_traj, gp_ll]   = exactInferenceWithLL(D(gp_test_mask), gp_params,'getLL',1);

                    like_fold = gp_ll;
                    % orthogonalize the trajectories

                    [gp_Xorth, gp_Corth]     = orthogonalize([gp_traj.xsm], gp_params.C);
                    gp_traj                  = segmentByTrial(gp_traj, gp_Xorth, 'data');

                    %Validation with LNO
                    cv_gpfa = cosmoother_gpfa_viaOrth_fast(D(gp_test_mask),gp_params,zDim(idim));
                    cv_gpfa_cell       = struct2cell(cv_gpfa);
                    cvdata             = cell2mat(cv_gpfa_cell(end,:));
                    
                    mse_fold           = sum(sum((cvdata - gp_test_data).^2));
                    
                    
                end
                if idim == 1 && ifold == 1
                    figure(200)
                    stairs(gp_test_data(1,:)), hold on
                end
                if ifold == 1 && zDim(idim) == 10 
                    figure(200)
                    plot(cvdata(1,:))
                    drawnow
                end
                
                % add up the likelihood and LNO errors across folds
                eval(sprintf('mse_%s(idim) = mse_%s(idim) + mse_fold;',s{1},s{1}))
                eval(sprintf('like_%s(idim) = like_%s(idim) + like_fold;',s{1},s{1}));
            end
            
        end
   end
   parfor_progress(0);
   
   figure(irat)
   set(gcf, 'position', [100, 100, 900, 500])
   subplot(1,2,1)
   plot(zDim,mse_gpfa - mean(mse_gpfa),'-bx'),hold on
   xlabel('Latent dimension with GPFA')
   ylabel('LNO MSE (centered)')
   set(gca, 'fontsize', 14)
   subplot(1, 2, 2)
   plot(zDim,like_gpfa - mean(like_gpfa),'-rx'),hold on 
   xlabel('Latent dimension with GPFA')
   ylabel('Posterior loglike (centered)')
   set(gca, 'fontsize', 14)
   
   set(gcf,'NextPlot','add'); axes;
   h = title(sprintf('%s bin(%d)',animals{irat},D(1).T));
   set(gca,'Visible','off');
   set(h,'Visible','on');
   
   drawnow
   saveas(gcf,sprintf('%s_valGPFA_bin(%d).png',roots{irat},D(1).T))
   
end


%% Plot the LV (gpfa_traj.data) and X-Y position to compute place fields in
%latent space
meantraj = zeros(1,2*D(1).T);
for idx_hdv = 1 : zDim
    for ilap = 1:length(D) 
        
       NoLap = D(ilap).trialId;   
       gp_traj  = exactInferenceWithLL(D(ilap), gp_params,'getLL',0); 
        
       x       = y_XT(y_int_maze(NoLap,1):y_int_maze(NoLap,2));
       y       = y_YT(y_int_maze(NoLap,1):y_int_maze(NoLap,2));
       spe     = y_speed(y_int_maze(NoLap,1):y_int_maze(NoLap,2));
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
       title(sprintf('Hidden variable %d, animal %s',idx_hdv, animal))
       hold on, grid on, box on
       xlabel('Space (mm)')
       ylabel('Space (mm)')
       zlabel('Amplitude') 
       campos([4875.4   -5094.3    22.1])
       plot(D(ilap).centers(1, :), D(ilap).centers(2, :), 'color', [0.4 0.4 0.4])
       plot3(D(ilap).centers(1, :), D(ilap).centers(2, :), gp_traj.xsm(idx_hdv,:),...
             'color',D(ilap).epochColors)
       plot3(D(ilap).centers(1, :), D(ilap).centers(2, :), spe_int./max(spe_int),...
             'color',[0.3 0.8 0.3])
       if strcmp(y_trial{D(ilap).trialId}, 'left')  
          meantraj(1:D(1).T) =  gp_traj.xsm(3,:) + meantraj(1:D(1).T) ;
       else
          meantraj(D(1).T + 1: end) =  gp_traj.xsm(idx_hdv,:) + meantraj(D(1).T + 1: end);
       end   
    end
end

figure
imagesc(mean(corX,3))
%% Check covarince kernel

Tdif         = repmat((1:D(1).T)', 1, D(1).T) - repmat(1:D(1).T, D(1).T, 1);
figure(2)
title(sprintf('GPs covariances, animal %s',animal))
colors = jet(50);
for idx_hdv = 1 : zDim
    subplot(2,5,idx_hdv)
    K            = (1 - gp_params.eps(idx_hdv)) * ...
                    exp(-gp_params.gamma(idx_hdv) / 2 * Tdif.^2) +...
                    gp_params.eps(idx_hdv) * eye(D(1).T);
    for j = 1:10
        X_sampled = mvnrnd(zeros(1, length(K)), K);      
        plot(X_sampled, 'color', colors(j,:)), hold on
    end
    %imagesc(K)
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


