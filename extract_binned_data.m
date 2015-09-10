%% Analysis of Buzsaki database
% wth SPWs file \01__hc5_maze06_002\__p111111-0101_SPWs.txt
% creates the matrices with different bin size for the model selection
clc, close all; clear all;
%ToDO: The median trajectory is noisy, should be smoothed. 
%ToDO: Animal 6 has different scale in x, y (divided by 100)
%ToDD: Plot spikes_xy one cell/lap to have a groundtrue of the spike count

basepath = {'/media/bigdata/i01_maze05.005/',...
           '/media/bigdata/i01_maze06.002/',...
           '/media/bigdata/i01_maze06.005/',...
           '/media/bigdata/i01_maze08.001/',...
           '/media/bigdata/i01_maze08.004/',...
           '/media/bigdata/i01_maze13.003/'              
            };
animal =  {'i01_maze05_MS.005',...
           'i01_maze06_MS.002',...
           'i01_maze06_MS.005',...
           'i01_maze08_MS.001',...
           'i01_maze08_MS.004',...
           'i01_maze13_MS.003'};
       
%height of the bin. depends on the variability of the animal's trajectories
height_bin = [200, 100, 150, 100, 100, 50];
       
       
num_segs  = 35:5:80; % number of segments to extract, is equivalent to bin_size

for rat = 1 : length(animal)
    fprintf('Processing animal %s\n', animal{rat});
    e_obj             = load([basepath{rat} animal{rat} '_BehavElectrData.mat']);
    e_clusters        = e_obj.Spike.totclu;
    e_laps            = e_obj.Laps.StartLaps(e_obj.Laps.StartLaps~=0); %@1250 Hz
    %Adding the end of the last lap because Laps only contains the start
    e_laps(end+1)     = e_obj.Par.SyncOff;
    e_XT              = e_obj.Track.X;
    e_YT              = e_obj.Track.Y;
    e_X               = e_obj.Spike.X;
    e_Y               = e_obj.Spike.Y;
    e_events          = e_obj.Par.MazeSectEnterLeft;
    e_Fs              = e_obj.Par.SamplingFrequency;
    e_isIntern        = e_obj.Clu.isIntern;
    e_numLaps         = numel(e_laps)-1;
    [e_spk, e_spk_lap, e_X, e_Y, e_X_lap, e_Y_lap]  = get_spikes(e_clusters, e_obj.Spike.res, e_laps, e_X, e_Y);

    e_N               = size(e_spk_lap,2);
    e_typetrial       = {'left', 'right', 'errorLeft', 'errorRight'};
    e_trialcolor      = hsv(5);
    
    %========== Extract spks when the mouse if running====================%
    %=====================================================================%

    for lap = 1:e_numLaps          

        %(b) Runing in the maze. Extracted based on the
        %EnterSection time stamp without considering left-right
        e_idx_run = [e_events{lap}(2,1), sum(e_events{lap}(7:8,1))];
        e_int_at_maze(lap, :) = e_idx_run;    
        e_length_run(lap) = (e_idx_run(2)-e_idx_run(1))/e_Fs;
        %sect 1:enter, 6:exit
        for neu=1:e_N
            e_idx = e_spk_lap{lap,neu}>=e_idx_run(1) & e_spk_lap{lap,neu}<=e_idx_run(end);
            e_X_Run_Lap{lap,neu}  = e_X_lap{lap, neu}(e_idx);
            e_Y_Run_Lap{lap,neu}  = e_Y_lap{lap, neu}(e_idx);
        end
        %Type of trial
        e_trial{lap}          = e_typetrial{e_obj.Laps.TrialType(e_laps(lap))};
        e_color(lap,:)        = e_trialcolor(e_obj.Laps.TrialType(e_laps(lap)),:);
    end
    
    %====================Extract median trajectory========================%
    %=====================================================================%
    disp('Starting spatial segmentation')
    %Segment base on spatial coordinates rather than time.
    %interpolate the position to the longest time
    [e_leftT, e_rightT, e_failed_trial] = protoLap(e_XT, e_YT, e_length_run, e_trial,...
                                    e_X_Run_Lap,e_Y_Run_Lap, e_int_at_maze,...
                                    e_Fs, animal{rat}, e_isIntern);
   

    %====================Spatial Segementation============================%
    %=====================================================================%                            
    for bi = 1 : length(num_segs) 
        %For validation purposes plot the trajectories and the rois
        figure(100+rat)
        subplot(2,5,bi), hold on         
        for ilap = 1 : e_numLaps
           %postion of the animal per lap 
           xt          = e_XT(e_int_at_maze(ilap,1):e_int_at_maze(ilap,2));
           yt          = e_YT(e_int_at_maze(ilap,1):e_int_at_maze(ilap,2));
           plot(xt, yt,'color',[0.6, 0.6, 1])
        end
        title(sprintf('bin %d',num_segs(bi)))
        axis([-50 1000 0 1100])
        
        
        %options to the function below
        roiDims      = [20 height_bin(rat)]; %width and length of ROI
        connectgrids = 1;
        ctrNeuron    = 0; % neuron to plot just see things are going OK
        show         = 0;
        verbose      = false;
        segments     = num_segs(bi);
        fprintf('Segmenting with %d bins\n', num_segs(bi));

        %count spikes inside grids and get normalized firing rate
        [e_rate, e_time_per_bin, e_centers] = normFiringRate(e_XT, e_YT, e_X_Run_Lap, e_Y_Run_Lap, e_int_at_maze,...
                              segments, show, connectgrids, roiDims, e_leftT,...
                              e_rightT, ctrNeuron, e_trial, verbose);
        
               
        e_X_pyr = e_Fs*e_rate(~e_isIntern,:,:);                  
        %remove failed trails
        e_laps_success = 1:e_numLaps;
        e_laps_success(e_failed_trial) = [];
        for ilap = 1 : numel(e_laps_success)
            real_lap = e_laps_success(ilap);
            e_D(ilap).data = e_X_pyr(:,:,real_lap);
            e_D(ilap).condition = e_trial{real_lap}; 
            e_D(ilap).epochColors = e_color(real_lap, :);
            e_D(ilap).trialId = real_lap;
            e_D(ilap).T = size(e_D(1).data,2);
            e_D(ilap).centers = e_centers(:,:,real_lap);
        end

        firing_thr      = 0.5; % Minimum norm. firing rate which 
                               % neurons should be kept
        m               = mean([e_D.data],2);
        keep_neurons    = m >= firing_thr;
        fprintf('%d neurons remained with firing rate above %2.2f Hz\n',...
                    sum(keep_neurons),firing_thr)
        % Remove low firing rate neurons
        for itrial = 1:length(e_D)
            e_D(itrial).y = e_D(itrial).data(keep_neurons,:);
        end
        %Prepare data for the datahigh library
        e_time_bin = mean(e_time_per_bin(:, e_laps_success)./1.25);
        fprintf('average bin time per lap %3.3fms (%3.3f - %3.3f)ms\n', mean(e_time_bin), min(e_time_bin), max(e_time_bin))
        txt = sprintf('%s_binned_at_%d.mat',[basepath{rat}  animal{rat}],num_segs(bi));
        save(txt,'e_D', 'e_leftT', 'e_rightT', 'e_time_bin')
        fprintf('Saved file %s\n\n', txt)
        
        drawnow 
        if show
            for i = 1 : e_numLaps
                figure(i)
                txt = sprintf('%s_binned_at_%d_cell_%d_lap_%d.png',[basepath{rat}  animal{rat}],num_segs(bi), ctrNeuron, i);
                saveas(gcf, txt)
                close
            end
        end
                         
    end
    clear e_*                            
end