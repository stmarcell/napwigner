import numpy as np

def get_spikes(clusters, spikes, laps):
    # GET_SPIKES extract variables per neuron per lap
    #
    #   [spk, spk_lap] = get_spikes(clusters, Spike, laps)
    #   If more inputs ar provided, then
    #   [spk, spk_lap, X, X_lap] = get_spikes(clusters, Spike, laps, X)

    #exVars = 0;
    #if nargin>3; exVars = nargin - 3; end

    nclust = max(clusters)
    nlaps = len(laps)-1
    spk_per_neuron = np.empty((nclust,), dtype=object)
    spk_per_lap = np.empty((nlaps,nclust), dtype=object)
    for j in range(0,max(clusters)): #for each neuron
        spk_per_neuron[j] = np.copy(spikes[clusters==j]);
        #for v = 1 : exVars
        #   eval(['aux' num2str(v) '{j} = varargin{v+3}(clusters==j);']); 
        #end
        for k in range(0,len(laps)-1): # for each lap
            index = (spk_per_neuron[j]>=laps[k]) & (spk_per_neuron[j]<laps[k+1]);
            spk_per_lap[k,j] = spk_per_neuron[j][index];
            #for v = 1 : exVars
            #    eval(['aux_lap' num2str(v) '{k,j} = aux' num2str(v) '{j}(index);']); 
            #end
        #end
    #end

    #varargout{1} = spk_per_neuron;
    #varargout{2} = spk_per_lap;

    #for v = 1 : exVars
    #    varargout{2 + v} = eval(['aux' num2str(v)]);
    #    varargout{2 + exVars + v} = eval(['aux_lap' num2str(v)]);
    #
    #end 
    
    return spk_per_neuron, spk_per_lap