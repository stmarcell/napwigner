function [spk_time_per_neuron,varargout] = get_spikes(clusters, spike_time, varargin)
% GET_SPIKES extract variables per neuron
%
%   [spk_time_per_neuron, ...] = get_spikes(clusters, spike_time, ...)
%
%   NOTE: indexing errors may occur if the neuron with the largest ID does
%   not fire

nXtra = nargin - 2;

spk_time_per_neuron = cell(1,max(clusters));
for j=1:max(clusters) %for each neuron  
    spk_time_per_neuron{j} = spike_time(clusters==j);
    for v = 1 : nXtra
        varargout{v}{j} = varargin{v}(clusters==j);
    end
end
