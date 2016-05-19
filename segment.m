function [F,keep_neurons,varargout] = segment(D, name_var, bin_size, Fs, keep_neurons, firing_thr, maxTime, varargin)
%SEGMENT remove low firing rate neurons and segments in non-overlapping
%        windows
%        namevar: is the name of the field in the structure D which is to be
%        segmented.
%        Fs: sampling frequency
%        keep_neurons: a vector indicating neurons to include '1' and exclude '0'
%        firing_thr: minimum firing rate find which neurons should be kept
%        bin_size: size of the segmentation bin (s)
%        maxTime: maximum segmenting time (s)
%        extra input variables of the same length as D.X are segmented
%        into extra output variables
%
%Ruben Pinzon@2015

useSqrt             = 1; % square root tranform for pre-processing?    
bin_nsamples        = ceil(bin_size * Fs); % bin size (Seconds * Fs) = samples
data                = [D.(name_var)];
firing_rate         = mean(data,2) * Fs;
nLaps               = length(D);

if nargin < 5
    keep_neurons = 1;
end
if nargin < 6
    firing_thr = 0;
end
if nargin < 7
    maxTime = 0;
end
nXtra               = length(varargin);

if length(keep_neurons) == 1
    keep_neurons    = firing_rate >= firing_thr;
else
    if size(data,1) ~= length(keep_neurons);
    	error('Incompatible vector of neurons');
    end
    disp('Vector of neurons to remove provided');
    keep_neurons    = keep_neurons & (firing_rate >= firing_thr);
    %firing_thr      = NaN;
end

yDim                = sum(keep_neurons);

fprintf('%d neurons remained with firing rate above %2.2f Hz out of %d\n',...
                yDim, firing_thr, length(keep_neurons))



% Remove low firing rate neurons
Temp = struct('data',cell(1,nLaps)); %repmat(struct([]),length(D),1);
for itrial = 1 : nLaps
    Temp(itrial).data = D(itrial).(name_var)(keep_neurons,:);
end
for j = 1:nXtra
      varargout{j} = cell(1,nLaps); %#ok<AGROW>
end

F=struct('y',cell(1, length(Temp)));
%Extract bins for one trial, since NOT all the trials
%are of the same duration
for ilap = 1 : nLaps
    seq         = [];
    T           = floor(size(Temp(ilap).data, 2) / bin_nsamples);
    if maxTime ~= 0
       T_requested = floor(maxTime * Fs / bin_nsamples); 
       if T_requested > T
           error('Requested time larger than lenght of trial')
       end
       T = T_requested;
    end

    for j = 1:nXtra
          varargout{j}{ilap} = zeros(T,1);
    end
    
    seq.y           = nan(yDim, T);
    for t = 1:T
      iStart        = bin_nsamples * (t-1) + 1;
      iEnd          = bin_nsamples * t;
      %seq.y(:,t)    = sum(Temp(ilap).data(:, iStart:iEnd), 2)./bin_size;
      seq.y(:,t)    = mean(Temp(ilap).data(:, iStart:iEnd), 2).*Fs;
      
      for j = 1:nXtra
          varargout{j}{ilap}(t) = mean(varargin{j}{ilap}(iStart:iEnd));
      end
    end
    %normalization with square root transform
    if useSqrt
        seq.y       = sqrt(seq.y);
    end
    
    F(ilap).trialId     = D(ilap).trialId;
    F(ilap).clusterId   = D(ilap).clusterId(keep_neurons);
    if isfield(D, 'type')
        F(ilap).type    = D(ilap).type;
    else
        F(ilap).type    = D(ilap).trialType;
    end
    F(ilap).behav       = D(ilap).behav;
    F(ilap).side        = D(ilap).side;
    F(ilap).y           = seq.y;
    F(ilap).T           = T;
    
end
