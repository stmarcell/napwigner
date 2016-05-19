import numpy as np

def add_spike_metadata(spikes, maze_events):
    # Find out the maze secion and the number of the lap based
    # on enter - leave events provided per lap:
    # maze_events[lap][section,event_type] where event_type is 
    # in [0, 1] for ente and leave respectively
    #
    #   [spk, spk_lap] = get_spikes(clusters, Spike, laps)
    #   If more inputs ar provided, then
    #   [spk, spk_lap, X, X_lap] = get_spikes(clusters, Spike, laps, X)

    #exVars = 0;
    #if nargin>3; exVars = nargin - 3; end

    spike_lap = np.zeros( len(spikes), dtype='uint8' ) # np.zeros( len(spikes), dtype=[('lap', 'uint8')])
    spike_sec = np.zeros( len(spikes), dtype='uint8' ) # np.zeros( len(spikes), dtype=[('section', 'uint8')])
    for lap in range(0,len(maze_events)):
        lap_events = maze_events[lap]
        id = lap_events[:,0]
        lap_enter = min(id[id>0])
        lap_leave = max(lap_events[:,1])
        print lap_enter, lap_leave
        spike_lap[(lap_enter<=spikes) & (spikes<=lap_leave)] = lap+1
        for section in range(0,len(lap_events)):
            sec_enter = lap_events[section,0]
            sec_leave = lap_events[section,1]
            spike_sec[(sec_enter<=spikes) & (spikes<=sec_leave)] = section+1
    
    return spike_lap, spike_sec