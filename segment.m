function [F,keep_neurons] = segment(D, name_var, bin_size, Fs, keep_neurons, firing_thr, maxTime)
%SEGMENT remove low firing rate neurons and segments in non-overlapping
%        windows
%        namevar: is the name of the field in the structure D which is to be
%        segmented.
%        Fs: sampling frequency
%        keep_neurons: a vector indicating neurons to include '1' and exclude '0'
%        firing_thr: minimum firing rate find which neurons should be kept
%        bin_size: size of the segmentation bin
%        maxTime: maximum segmenting time s 
%
%Ruben Pinzon@2015

useSqrt             = 1; % square root tranform for pre-processing?    
bin_width           = ceil(bin_size * Fs); % bin size (Seconds * Fs) = samples
data                = [D.(name_var)];
firing_rate         = mean(data,2) * Fs;

if nargin < 5
    keep_neurons = 1
end
if nargin < 6
    firing_thr = 0
end
if nargin < 7
    maxTime = 0
end

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
Temp = repmat(struct([]),length(D),1);
for itrial = 1 : length(D)
    Temp(itrial).data = D(itrial).(name_var)(keep_neurons,:);
end

%Extract bins for one trial, since NOT all the trials
%are of the same duration
for ilap = 1 : length(Temp)
    seq         = [];
    T           = floor(size(Temp(ilap).data, 2) / bin_width);
    if maxTime ~= 0
       T_requested = floor(maxTime * Fs / bin_width); 
       if T_requested > T
           error('Requested time larger than lenght of trial')
       end
       T = T_requested;
    end
    
    seq.y           = nan(yDim, T);
    for t = 1:T
      iStart        = bin_width * (t-1) + 1;
      iEnd          = bin_width * t;
      seq.y(:,t)    = sum(Temp(ilap).data(:, iStart:iEnd), 2)./bin_size;
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
    F(ilap).y           = seq.y;
    F(ilap).T           = T;
    
end
