% Script for meeting 22 Sept 2015 that show the general pre and processing
% of running section in the HC-5 database with GPFA
%

clc, close all; clear all;
cd /media/LENOVO/HAS/CODE/Wigner-Pattern

basepath        = '/media/bigdata/';
[files, animals, roots]= get_matFiles(basepath);

%======================================================================%
%================       Figure 1: the task            =================%
%======================================================================%
animal = 2; %choose one animal to show the task

y_obj             = load(files{animal});
y_clusters        = y_obj.Spike.totclu;
y_laps            = y_obj.Laps.StartLaps(y_obj.Laps.StartLaps~=0); 
y_laps(end+1)     = y_obj.Par.SyncOff;
y_wheelspeed      = y_obj.Laps.WhlSpeedCW;
y_eeg              = y_obj.Track.eeg;
y_XT              = y_obj.Track.X;
y_YT              = y_obj.Track.Y;
y_X               = y_obj.Spike.X;
y_Y               = y_obj.Spike.Y;
y_events          = y_obj.Par.MazeSectEnterLeft;
Fs                = y_obj.Par.SamplingFrequency;
y_speed           = y_obj.Track.speed;
y_isIntern        = y_obj.Clu.isIntern;
y_spks            = y_obj.Spike.res;
num_laps           = numel(y_laps)-1;
%
close all
pause(3)

figure(1)
subplot(2,5,[1 2 6 7])
%show all the trajectories
set(gcf, 'position', [100 185 1600 702], 'color', 'w')
for ilap = 1 : num_laps
   lap_st   = y_laps(ilap);
   lap_en   = y_laps(ilap+1);
   
   x        = y_XT(lap_st:lap_en);
   y        = y_YT(lap_st:lap_en);
   plot(x, y, 'color', [0.6 0.6 0.6]), hold on     
end
set(gca, 'fontsize', 14)
axis([-50 1300 0 1050])
xlabel('x (mm)', 'fontsize', 14)
ylabel('y (mm)','fontsize', 14), grid on
subplot(2,5,3)
axis([0 1 -100 1000])
set(gca, 'fontsize', 14, 'xtick', [])
ylabel('Wheel speed (m/s)')

subplot(2,5,8)
axis([0 1 -100 1000])
set(gca, 'fontsize', 14, 'xtick', [0 1])
xlabel('time (s)'), box on
ylabel('Animal speed (m/s)')

subplot(2,5,[4 5 9 10])
set(gca, 'fontsize', 14, 'yaxislocation', 'right')
xlabel('time (s)'), box on
ylabel('Cells')

pause(1)
%shown the evolution of one left and ne right trajectories
colors  = jet(110);
colors(2,:) = []; %delete green
for ilap = 1:2
    lap_st   = y_laps(ilap);
    dt       = 250; %200 ms
    h        = [];
       
    for t = 1:dt/2:y_laps(ilap+1) - y_laps(ilap)
        if ~isempty(h)
            delete(h)
        end
        dt_x    = y_XT(lap_st + t:lap_st + t + dt);
        dt_y    = y_YT(lap_st + t:lap_st + t + dt);
        dt_wh   = y_wheelspeed(lap_st + t:lap_st + t + dt);
        dt_sp   = y_speed(lap_st + t:lap_st + t + dt);
        figure(1)
        subplot(2,5,3)
        time    = linspace(lap_st +t,lap_st +t +dt, length(dt_wh))./1250;
        plot(time, dt_wh,'r', 'linewidth', 2), 
        axis([time(1) time(end) -100 1000])        
        ylabel('Wheel speed (m/s)')
        set(gca, 'fontsize', 14, 'xtick', [])
        
        subplot(2,5,8)
        time    = linspace(lap_st +t,lap_st +t +dt, length(dt_wh))./1250;
        plot(time, dt_sp,'r', 'linewidth', 2), 
        axis([time(1) time(end) -100 1000])        
        ylabel('Animal speed (m/s)')
        xlabel('time(s)')
        set(gca, 'fontsize', 14, 'xtick', [time(1)  mean(time([end, 1])) time(end)])
        set(gca,'XTickLabel',sprintf('%2.1f\n',[time(1)  mean(time([end, 1])) time(end)]))
       
        subplot(2,5,[1 2 6 7])
        h = plot(dt_x, dt_y, 'r', 'linewidth', 3);
        
        
       subplot(2,5,[4 5 9 10])
       for j = lap_st + t : lap_st + t + dt
           spikes   = find(y_spks==j);
           if ~isempty(spikes)
              cells    = y_clusters(spikes);
              cells(y_isIntern(cells)==1) = [];
              if ~isempty(cells)
                  for c = 1 : numel (cells)
                      plot(j/1250*[1 1],[cells(c) cells(c)+1.5],'color', colors(cells(c),:), 'linewidth',3), hold on
                  end
              end
              
           end
       end
       axis([time(1) time(end) 0 110])
       set(gca, 'fontsize', 14, 'xtick', [time(1)  mean(time([end, 1])) time(end)])
       set(gca,'XTickLabel',sprintf('%2.2f\n',[time(1)  mean(time([end, 1])) time(end)]))

       xlabel('time (s)'), box on
       ylabel('Cells')
       hold off 
        drawnow
        %pause(0.05)
        
        
    end
    delete(h)
    pause(1)
end

%% Show one runing sequence and a SPWs

idx_run            = [y_events{1}(2,1), y_events{1}(2,1)+1250]; 

figure(1)
subplot(3,1,[2 3])
set(gcf, 'position', [100 185 700 900], 'color', 'w')
colors  = jet(110);
for j = idx_run(1) : idx_run(2)
   spikes   = find(y_spks==j);
   if ~isempty(spikes)
      cells    = y_clusters(spikes);
      cells(y_isIntern(cells)==1) = [];
      if ~isempty(cells)
          for c = 1 : numel (cells)
              plot(j/1250*[1 1],[cells(c) cells(c)+1],'color', colors(cells(c),:), 'linewidth',2), hold on
          end
      end
   end
end
set(gca,'Fontsize',14)
xlabel('Time (seconds)')
title('Spike rater (only pyr cells)')
xlim(idx_run./1250)
subplot(3,1,1)
plot((idx_run(1):idx_run(2))./1250, y_eeg(idx_run(1):idx_run(2)))
xlim(idx_run./1250)
title('EEG (Running sect.)')
set(gca,'Fontsize',14)
saveas(gcf,sprintf('%s_running.png',roots{animal}))


idx_spw             = floor([849279.5693	849501.9818]*1.250);
figure(2)
subplot(3,1,[2 3])
set(gca,'Fontsize',14)
xlabel('Time (seconds)')
title('Spike rater (only pyr cells)')
set(gcf, 'position', [100 185 700 900], 'color', 'w')

for j = idx_spw(1) : idx_spw(2)
   spikes   = find(y_spks==j);
   if ~isempty(spikes)
      cells    = y_clusters(spikes);
      cells(y_isIntern(cells)==1) = [];
      for c = 1 : numel (cells)
              plot(j/1250*[1 1],[cells(c) cells(c)+1],'color', colors(cells(c),:), 'linewidth',2), hold on
      end
   end
end
xlim(idx_spw./1250)
subplot(3,1,1)
plot((idx_spw(1):idx_spw(2))./1250, y_eeg(idx_spw(1):idx_spw(2)))
xlim(idx_spw./1250)
title('EEG (SPW)')
set(gca,'Fontsize',14)
saveas(gcf,sprintf('%s_spws.png',roots{animal}))

%% Toy construction for GPFA
%See python file covFuncToy.py
%% run a single GPFA
clc, close all; clear all;
cd /media/LENOVO/HAS/CODE/Wigner-Pattern

basepath                = '/media/bigdata/';
[files, animals, roots] = get_matFiles(basepath, '/*).mat');
useSqrt                 = 1; % square root tranform?    
zDim                    = 10; %hidden space dimensions
rat_curr                = '';
for irat = 1 : length(files)
    
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
   post_ll = zeros(length(zDim), length(cv_fold_idx)-1 );
   mse_fold = zeros(length(zDim), length(cv_fold_idx)-1 );
   model = {};
   for idim = 1 : length(zDim)
       tic
        for ifold = 1 : num_folds         
            fprintf('fold %d', ifold)
            gp_test_mask = cv_mask;
            gp_test_mask(cv_trials(cv_fold_idx(ifold):...
                             cv_fold_idx(ifold+1)-1)) = true;
            gp_train_mask = ~gp_test_mask; 
            gp_train_data = D(gp_train_mask);
            gp_test_data = D(gp_test_mask);

            %training of the GPFA with already binned data (1)
            figure(1)
            %show all the trajectories
            set(gcf, 'position', [100 185 1210 702], 'color', 'w')
            pause(5)
            [gp_params, gpfa_traj, ll, llt, diffC] = gpfa_mod(gp_train_data,zDim(idim),...
                                                 'bin_width', 1);
            %Posterior of test data given the trained model
            [gp_traj, gp_ll]   = exactInferenceWithLL(gp_test_data, gp_params,'getLL',1);
            
            post_ll(idim, ifold) = gp_ll;
            % orthogonalize the trajectories

            [gp_Xorth, gp_Corth]     = orthogonalize([gp_traj.xsm], gp_params.C);
            gp_traj            = segmentByTrial(gp_traj, gp_Xorth, 'data');
            gp_traj            = rmfield(gp_traj, {'Vsm', 'VsmGP', 'xsm'});

            %Validation with LNO
            cv_gpfa = cosmoother_gpfa_viaOrth_fast(gp_test_data,gp_params,zDim(idim));
            cv_gpfa_cell       = struct2cell(cv_gpfa);
            gp_y_est           = cell2mat(cv_gpfa_cell(end,:));
            gp_y_real          = [gp_test_data.y];
            gp_xx              = (gp_y_est - gp_y_real).^2;
            
            mse_fold(idim, ifold) = sum(sum(gp_xx));
            
            
            model{idim, ifold}.GPparams = gp_params;
            model{idim, ifold}.GPtrajec = gpfa_traj;
            model{idim, ifold}.y_est    = gp_y_est;
            model{idim, ifold}.error    = gp_xx;
            model{idim, ifold}.training = gp_train_mask;
            model{idim, ifold}.test     = gp_test_mask;
            
            clear gp* 
        end
        fprintf('Rat %s (%d/%d),zDim %d,loglike %5.2f, mse %5.2f, done in %3.2fs\n',...
                     animals{irat},irat, length(files), zDim(idim),mean(post_ll(idim, :)),mean(mse_fold(idim, :)), toc)
   end
end

%% Figure 2: sparseness

clc, close all; clear all;
cd /media/LENOVO/HAS/CODE/Wigner-Pattern

basepath                = '/media/bigdata/';
[files, animals, roots] = get_matFiles(basepath, '/*).mat');
useSqrt                 = 1; % square root tranform?    
rat_curr                = '';
bin_sizes               = 40:5:85;

for irat = 51 : 60
    
   y_obj            = load(files{irat});
   D                = y_obj.D; 
   num_laps         = length(D);
   num_cells        = size(D(1).y, 1);
   sparse_ratio (irat-50, :) = [D.sparsness];
   
end

figure(2)
set(gcf, 'position', [100 185 970 600], 'color', 'w')
plot(bin_sizes, sparse_ratio, 'color', [0.7 0.7 0.7]), hold on
plot(bin_sizes, mean(sparse_ratio,2), 'color', [1 0.3 0.3], 'linewidth', 2)
xlabel('number of bins')
ylabel('Ratio zero elements/all elements')
title(sprintf('Sparseness ratio animal %s', animals{irat}))
set(gca,'fontsize', 14 )
grid on
xlim([bin_sizes(1) bin_sizes(end)])

saveas(gcf,sprintf('%s_sparseness.png',roots{irat}))

%% Postetior likelihood for EID

clc, close all; clear all;
cd /media/LENOVO/HAS/CODE/Wigner-Pattern

basepath                = '/media/bigdata/';
[files, animals, roots] = get_matFiles(basepath, '/*CV*');
useSqrt                 = 1; % square root tranform?    
rat_curr                = '';
bin_sizes               = 40:5:85;
zDim                    = [2:2:20]; %hidden space dimensions
colors                  = jet(10);

figure(2)
set(gcf, 'position', [600 185 970 600], 'color', 'w')
figure(1)
set(gcf, 'position', [100 185 970 600], 'color', 'w')
ibin2 = 0;
for ibin = 1 : 10%length(files)
   
   ibin2 = ibin2 + 1;
   r = load(files{ibin});    
   LL(ibin2,:) = mean(r.R.resultLL,2); 
   
   LNO(ibin2,:) = mean(r.R.resultLNO,2);
   
   figure(1)
   plot(zDim, LL(ibin2,:)-mean(LL(ibin2,:)),'-s',...
       'displayname', sprintf('%d bins', bin_sizes(ibin2)),...
       'color', colors(ibin2,:)), hold on
   text(14, 0.006-ibin2*35, sprintf('%d bins %2.2e',bin_sizes(ibin2), LL(ibin2,1)),...
       'color', colors(ibin2,:),'fontsize', 14)
   
   figure(2)
   plot(zDim, LNO(ibin2,:)-mean(LNO(ibin2,:)),'-s',...
       'displayname', sprintf('%d bins', bin_sizes(ibin2)),...
       'color', colors(ibin2,:)), hold on
   text(14, 15e4-ibin2*1.5e4, sprintf('%d bins %2.2e',bin_sizes(ibin2), LNO(ibin2,1)),...
       'color', colors(ibin2,:),'fontsize', 14)
end

ylabel('MSE (centered)')
xlabel('Latent dimension')
set(gca,'fontsize', 14 )
xlim([0 zDim(end)])
grid on
title(sprintf('MSE LNO val. %s',animals{ibin}))
saveas(gcf,sprintf('%s_MSE_LNO.png',roots{ibin}))

figure(1)
ylabel('Posterior loglike (centered)')
xlabel('Latent dimension')
set(gca,'fontsize', 14 )
xlim([0 zDim(end)])
grid on
title(sprintf('Posterior LogLike %s',animals{ibin}))
saveas(gcf,sprintf('%s_loglike.png',roots{ibin}))


%% 
%D
load('/media/bigdata/i01_maze06.002/i01_maze06.002_spikecount_bin_(40).mat'); 
conditions = strsplit([D.condition],'t');
idx = 0;
for cell = 1:40
    y10 = R.gpfa{5,1}.y_est;
    y2 = R.gpfa{1,1}.y_est;
    y = R.gpfa{1,1}.y_real;

    num_test_laps = sum(R.gpfa{5,1}.test);
    totalLaps = numel(R.gpfa{5,1}.test);

    bins = length(R.gpfa{5,1}.Xorth);
    idx = idx + 1;
    subplot(2,5,idx)
    stairs(y(cell, :)), hold on
    plot(repmat([1 : num_test_laps - 1]*bins,2,1),ylim,'color',[0.7 0.7 0.7]);
    plot(y10(cell, :),'r')
    plot(y2(cell, :),'g')
    axis([0 length(y) 0 30])
    if idx == 10
        idx = 1;
        figure()
    end
end

test_mask = r.R.gpfa{1,1}.test;
cnt = 1;
for ilap = 1 : totalLaps-1
    
   if test_mask(ilap) == 1
      
      if strcmp(conditions{ilap},'righ') 
          figure(1)
          stairs(y(cell, cnt:cnt + bins-1), 'b'), hold on
          plot(y10(cell, cnt:cnt + bins-1),'r')
          plot(y2(cell, cnt:cnt + bins-1),'g')
      else
          figure(2)
          stairs(y(cell, cnt:cnt + bins-1),'b'), hold on
          plot(y10(cell, cnt:cnt + bins-1),'r')
          plot(y2(cell, cnt:cnt + bins-1),'g')
          
      end
      cnt = cnt + bins; 
   end
end
