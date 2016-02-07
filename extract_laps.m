function D = extract_laps(Fs, spk, keep_neurons, X, Y, speed, wh_speed, events, event_range, TrialType, BehavType)
%EXTRACT LAP Takes the HC-5 database and divides the vectors into laps. More
%            information about the database in
%            crcns.org/files/data/hc-5/crcns-hc5-data-description.pdf
%
%            INPUTS:
%            Fs: sampling frequency (Default: 1250)
%            spk: spikes times
%            keep_neuron: whether spike train is to be generated
%            X:  position of the animal in the x-axis for the whole experiment
%            Y:  position of the animal in the y-axis for the whole experiment
%            speed: speed of the animal in during the whole experiment
%            wh_speed: speed of the runing wheel
%            events: vector containing the time start/end for each lap and section
%                    in the maze
%            event_range: struct with fields in and out containing section IDs
%                         that describe entering and leaving the area of interest
%            TrialType: vector indicating the type of each lap, options are
%                   1: right alternation, 2: left alt., 3: wrong right alt.
%                   4: wrong left alt.
%            BehavType: vector indicating the behavior in each lap, options are
%                   1: first, 2: regular, 3: uncertain
%
%
%see also branch2, branch2_cleaned.m
%Revision:  Feb05: simplified things and variables
%Marcell Stippinger, 2016


numLaps         = length(events);
n_cells         = length(spk);
if length(keep_neurons) == 1
    n_kept      = true(1,n_cells);
else
    n_kept      = sum(keep_neurons);
end
%n_pyrs          = sum(isIntern==0);
%kernel          = gausswin(0.1*Fs);
color           = hsv(4);
side_select     = [1 2 2 1];

% Extract spks when the mouse is running 
for lap = 1:numLaps  
    
    % in some rare cases both left and right sections are visited
    % but non-visited section have entering and leaving time "0"
    [sect_in,  time_in]  = get_section_id(event_range.in, events{lap}(:,1));
    [sect_out, time_out] = get_section_id(event_range.out, events{lap}(:,2));
    %idx_lap      = [events{lap}(1,1), max(events{lap}(:,2))];
    idx_lap       = [min(time_in), max(time_out)];
    if sect_in == 13
    %for wheel section extract spikes when the wheel is moving    
        wheelNonZero    = find(D(lap).wh_speed~=0);
        if isempty(wheelNonZero)
            sprintf('Skipped lap %d without wheel run\n',lap)
            continue;
        end
        idx_lap         = [wheelNonZero(1), wheelNonZero(end)];
        D(lap).([name '_wheelNonZero']) = wheelNonZero;
    end
    
    X_lap        = X(idx_lap(1):idx_lap(2));
    Y_lap        = Y(idx_lap(1):idx_lap(2));
    acc_dst_lap  = [0; cumsum(sqrt((X_lap(2:end) - X_lap(1:end-1)).^2 + ...
                                   (Y_lap(2:end) - Y_lap(1:end-1)).^2))];
    speed_lap    = speed(idx_lap(1):idx_lap(2));
    wh_speed_lap = wh_speed(idx_lap(1):idx_lap(2));

    
    t_lap        = idx_lap(2) - idx_lap(1) + 1;
    cnt          = 0;
    %firing       = zeros(n_kept, t_lap); 
    spk_train    = zeros(n_kept, t_lap);
    for neu=1:n_cells
        if keep_neurons(neu)
            tmp              = zeros(1, t_lap); 
            cnt              = cnt + 1;
            
            %this re-filtering as introduced because get_spikes only uses
            %lap start times and idx_lap(:) might differ
            idx              = spk{neu}>=idx_lap(1) & spk{neu}<=idx_lap(2);
            %aligned to the start of the section            
            spikes_lap{cnt}      = spk{neu}(idx) - idx_lap(1) + 1;
            tmp(spikes_lap{cnt}) = 1; 
            %convolve the spike trains with a gauss filter 100 m
            %firing(cnt,:)    = Fs*conv(tmp,kernel, 'same');
            spk_train(cnt, :) = tmp; 
        end
    end        
   
    %Type of trial
    sec_lap                   = events{lap}-events{lap}(1,1)+1;
    sec_lap(sec_lap<0)        = 0;
    behav_lap                 = BehavType(lap);
    type_lap                  = TrialType(lap);
    side_lap                  = side_select(TrialType(lap));
    
    D(lap).trialId            = lap;
    D(lap).spikes             = spikes_lap;
    D(lap).clusterId          = find(keep_neurons);
    D(lap).X                  = X_lap;
    D(lap).Y                  = Y_lap;
    D(lap).speed              = speed_lap;
    D(lap).wh_speed           = wh_speed_lap;
    D(lap).sections           = sec_lap;
    D(lap).behav              = behav_lap;
    D(lap).type               = type_lap;
    D(lap).side               = side_lap;
    %D(lap).color              = color(TrialType(laps(lap)),:);
    D(lap).acc_dist           = acc_dst_lap;
    %D(lap).firing_rate        = firing;
    D(lap).spike_train        = spk_train;
    D(lap).duration           = idx_lap(2) - idx_lap(1);
    D(lap).start              = idx_lap(1);
    clear spikes *_lap tmp
end    

