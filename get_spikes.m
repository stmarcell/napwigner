function [spk_per_neuron,varargout] = get_spikes(clusters, spikes, varargin)
% GET_SPIKES extract variables per neuron per lap
%
%   [spk, spk_lap] = get_spikes(clusters, spike)

nXtra = nargin - 2;

for j=1:max(clusters) %for each neuron  
    spk_per_neuron{j} = spikes(clusters==j);
    for v = 1 : nXtra
        varargout{v}{j} = varargin{v}(clusters==j);
    end
end
