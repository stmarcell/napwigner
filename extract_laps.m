function [D, varargout] = extract_laps(Fs, spk, keep_neurons, events, section_range, TrialType, BehavType, X, Y, speed, wh_speed, varargin)
% EXTRACT_LAPS Takes the HC-5 database and divides the vectors into laps.
% More formation about the database in
% crcns.org/files/data/hc-5/crcns-hc5-data-description.pdf
%
% INPUTS:
%   Fs: sampling frequency (Default: 1250)
%   spk: spike times per neuron
%   keep_neuron: whether spike train is to be generated (bool array)
%   X: position of the animal in the x-axis for the whole experiment
%   Y: position of the animal in the y-axis for the whole experiment
%   speed: speed of the animal in during the whole experiment
%   wh_speed: speed of the runing wheel
%   events: vector containing the time start/end for each lap and section
%           in the maze
%   event_range: struct with fields in and out containing section IDs that
%                describe entering and leaving the area of interest
%   TrialType: vector indicating the type of each lap, options are
%              1: right alternation, 2: left alt., 3: wrong right alt.
%              4: wrong left alt.
%   BehavType: vector indicating the behavior in each lap, options are
%              1: first, 2: regular, 3: uncertain
%   extra input variables of the same length as X are extracted per lap
%   into extra output variables
% 
%
% See also branch2, branch2_cleaned.m

% Revision:  Feb05: simplified things and variables
% Marcell Stippinger, 2016


numLaps         = length(events);
n_cells         = length(spk);
nXtra           = length(varargin);
if length(keep_neurons) == 1
    n_kept      = true(1,n_cells);
else
    n_kept      = sum(keep_neurons);
end
%kernel          = gausswin(0.1*Fs);
side_select     = [1 2 2 1];

% Extract spks when the rat is running in the sections [section_range]
D = repmat(struct([]), numLaps, 1);
for lap = 1:numLaps  
    
    % in some rare cases both left and right sections are visited
    % but non-visited section have entering and leaving time "0"
    [sect_in,  time_in]  = get_section_id(section_range.in, events{lap}(:,1));
    [sect_out, time_out] = get_section_id(section_range.out, events{lap}(:,2));
    %idx_lap      = [events{lap}(1,1), max(events{lap}(:,2))];
    idx_lap       = [min(time_in), max(time_out)];
    if sect_in == 13
    % for wheel section extract spikes when the wheel is moving
    % TODO: move test to a variable "condition" because it's more general
    %       apply it to all sections; verify whether we need a continuous
    %       interval or moments with moving wheel
        wheelNonZero    = find(wh_speed(idx_lap(1):idx_lap(2))~=0);
        if isempty(wheelNonZero)
            fprintf('Skipped lap %d without wheel run\n',lap);
            continue;
        end
        idx_diff_lap    = [wheelNonZero(1), wheelNonZero(end)];
        idx_lap(1)      = idx_lap(1) + idx_diff_lap(1) - 1;
        idx_lap(2)      = idx_lap(1) + idx_diff_lap(2) - idx_diff_lap(1);
        %D(lap).wheelNonZero = wheelNonZero;
    end
    
    X_lap        = X(idx_lap(1):idx_lap(2),1); %second axis needed for synth data
    Y_lap        = Y(idx_lap(1):idx_lap(2),1);
    acc_dst_lap  = [0; cumsum(sqrt((X_lap(2:end) - X_lap(1:end-1)).^2 + ...
                                   (Y_lap(2:end) - Y_lap(1:end-1)).^2))];
    speed_lap    = speed(idx_lap(1):idx_lap(2));
    wh_speed_lap = wh_speed(idx_lap(1):idx_lap(2));

    % Collect requested neurons
    t_lap        = idx_lap(2) - idx_lap(1) + 1;
    cnt          = 0;
    %firing       = zeros(n_kept, t_lap); 
    spk_train    = zeros(n_kept, t_lap);
    spikes_lap   = cell(sum(keep_neurons),1);
    for neu=1:n_cells
        if keep_neurons(neu)
            tmp              = zeros(1, t_lap); 
            cnt              = cnt + 1;
            
            %filtering was introduced because get_spikes only uses
            %lap start times and idx_lap(:) might differ
            idx              = spk{neu}>=idx_lap(1) & spk{neu}<=idx_lap(2);
            %align to the start of the section            
            spikes_lap{cnt}      = spk{neu}(idx) - idx_lap(1) + 1;
            tmp(spikes_lap{cnt}) = 1; 
            %convolve the spike trains with a gauss filter 100 ms
            %firing(cnt,:)    = Fs*conv(tmp, kernel, 'same');
            spk_train(cnt, :) = tmp; 
        end
    end
    
    for i=1:nXtra
        varargout{i}{lap} = varargin{i}(idx_lap(1):idx_lap(2)); %#ok<AGROW>
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
    D(lap).acc_dist           = acc_dst_lap;
    %D(lap).firing_rate        = firing;
    D(lap).spike_train        = spk_train; %repmat(spk_train,1,2);
    D(lap).duration           = idx_lap(2) - idx_lap(1);
    D(lap).start              = idx_lap(1);
    clear spikes *_lap tmp
end    

