import numpy as np
import warnings

# The ordering of indices is he following, wherever possible
# Lap ID > Cluster ID > Time > Position > Others
# it determines the shape of the arrays.


# CONST based on http://code.activestate.com/recipes/65207-constants-in-python/?in=user-97991
# Put in const.py...:
class const(object):
    class ConstError(TypeError): pass
    def __setattr__(self,name,value):
        if name in self.__dict__:
            raise self.ConstError("Can't rebind const(%s)"%name)
        self.__dict__[name]=value
# you may also want to add the obvious __delattr__
#import sys
#sys.modules[__name__]=_const()


hc5 = const()
hc5.name            = "hc5"
hc5.TrialType_tx    = np.array([None, 'left', 'right', 'errorLeft', 'errorRight'])
hc5.SideType_tx     = np.array([None, 'left', 'right', 'right', 'left'])
hc5.BehavType_tx    = np.array(['first', 'regular', 'uncertain'])
hc5.Section_tx      = np.array([None, 'mid_arm', 'preturn', 'turn', 'turn', 'lat_arm', 'lat_arm',
                         'reward', 'reward', 'delay', 'delay', 'delay', 'delay', 'wheel'])
hc5.Side_tx         = np.array([None, 'center', 'center', 'right', 'left', 'right', 'left',
                         'right', 'left', 'center', 'center', 'center', 'center', 'center'])
hc5.Section_count   = len(hc5.Section_tx[1:])

ubb = const()
ubb.name            = "ubb"
ubb.TrialType_tx    = np.array([None, 'free_run'])
ubb.SideType_tx     = np.array([None, 'free_run'])
ubb.BehavType_tx    = np.array(['first', 'regular', 'uncertain'])
ubb.Section_tx      = np.array([None, 'field'])
ubb.Side_tx         = np.array([None, 'center'])
ubb.Section_count   = len(ubb.Section_tx[1:])

class PlaceField(object):
    pass
        
def cartesian(x,y):
    return np.column_stack((np.tile(x, len(y)), np.repeat(y, len(x))))
def poisson_continued(x, mu):
    from scipy.special import gamma
    return np.exp(-mu)*np.power(mu,x)/gamma(x+1)
def expand_dims(arr, axes):
    if type(axes) is int:
        return np.expand_dims(arr,axes)
    else:
        for axe in axes:
            arr = np.expand_dims(arr,axe)
        return arr
        
def _get_2d_field_mapper(par):
    numtics, posX, posY, spkX, spkY, spk_clust, clust_names, spatial_bin, freq = par
    return _get_2d_field(numtics, posX, posY, spkX, spkY, spk_clust, clust_names, spatial_bin, freq)

def _get_2d_field(numtics, posX, posY, spkX, spkY, spk_clust, clust_names, spatial_bin, freq):
    #Py_BEGIN_ALLOW_THREADS
    # From now on only read external variables
    pos_histo, sepx, sepy = np.histogram2d(posX, posY, spatial_bin)
    lenx, leny = pos_histo.shape
    nclust = len(clust_names)
    avgr = np.zeros((nclust))
    rate = np.empty((nclust,lenx,leny))
    info = np.empty((nclust,lenx,leny))
    for c in range(0,nclust):
        sel = spk_clust == clust_names[c]
        spk_histo, sepx, sepy = np.histogram2d(spkX[sel], spkY[sel], spatial_bin)
        avgr[c] = sum(sel) * freq / numtics
        rate[c,:,:] = np.divide(spk_histo,pos_histo) * freq
        info[c,:,:] = np.multiply(spk_histo,(np.log(rate[c,:,:])-
                                             np.log(sum(sel)/numtics)) ) * freq
    #Py_END_ALLOW_THREADS
    return avgr, rate, info, pos_histo
    
def _get_1d_field_mapper(par):
    numtics, posX, spkX, spk_clust, clust_names, spatial_bin, freq = par
    return _get_1d_field(numtics, posX, spkX, spk_clust, clust_names, spatial_bin, freq)

def _get_1d_field(numtics, posX, spkX, spk_clust, clust_names, spatial_bin, freq):
    #Py_BEGIN_ALLOW_THREADS
    # From now on only read external variables
    pos_histo, sepx = np.histogram(posX, spatial_bin)
    lenx = pos_histo.shape[0]
    nclust = len(clust_names)
    avgr = np.zeros((nclust))
    rate = np.empty((nclust,lenx))
    info = np.empty((nclust,lenx))
    for c in range(0,nclust):
        sel = spk_clust == clust_names[c]
        spk_histo, sepx = np.histogram(spkX[sel], spatial_bin)
        avgr[c] = sum(sel) * freq / numtics
        rate[c,:] = np.divide(spk_histo,pos_histo) * freq
        info[c,:] = np.multiply(spk_histo,(np.log(rate[c,:])-
                                             np.log(sum(sel)/numtics)) ) * freq
    #Py_END_ALLOW_THREADS
    return avgr, rate, info, pos_histo

def bayesian_decode(spk_count, place_rate, place_proba, clust_rate = None):
    '''Infer the position of the animal from spike counts, samples are time-binned spk count vectors.
       :param spk_count: spike counts, (clusters) x (samples)
       :param place_rate: spike rate at the given spot, (clusters) x [place bins]
       :param place_proba: probability of occuring at the spot, [place bins]
       :param clust_rate: average firing rate of the neuron, optional, (clusters)
       If clust_rate is omitted, it will be calculated by simply averaging place_rate.
       Rates must be in accordance with the time bins, poisson distribution assumed for spk_count.
       :return: inferred probabilities, (samples) x [place bins]
    '''
    from scipy.stats import poisson
    # Choose pdf calculating method
    pdf = poisson.pmf
    #pdf = poisson_continued
    # Shapes
    nclust, nsample = spk_count.shape
    place_dim = place_proba.ndim
    place_shape = place_proba.shape
    result_shape = np.concatenate(([nsample],place_shape))
    # Check broadcasting
    if nclust != len(place_rate):
        raise ValueError('the instantaneous rate and the mean_rate must have the same first dimension')
    if (clust_rate is not None) and nclust != len(clust_rate):
        raise ValueError('the instantaneous rate and the clust_rate must have the same first dimension')
    if place_rate.shape[1:] != place_shape:
        raise ValueError('the place probabilities must broadcast well to the place rates')
    # Likelihood function P(spike_count_vector|spatial_position)
    proba_cond_place = np.zeros(result_shape)
    for c in range(0,nclust):
        spk_aligned = expand_dims(spk_count[c],tuple(range(1,1+place_dim)))
        clust_proba = pdf(spk_aligned, place_rate[c])
        #clust_proba[clust_proba==0] = 1
        proba_cond_place += np.log(clust_proba)
    # Normalization factors for the sampling bins, (clusters) x (samples)
    if clust_rate is None:
        mean_param = np.nanmean(place_rate,axis=tuple(range(1,1+place_dim)))
    else:
        mean_param = clust_rate
    tot_proba = pdf(spk_count,mean_param[:,np.newaxis])
    #tot_proba[tot_proba==0] = 1
    # Normalization factors for the sampling bins, (samples,)
    tot_proba = np.sum(np.log(tot_proba), axis=0)
    # Bayesian decoding
    return proba_cond_place - tot_proba[:,np.newaxis] + np.log(place_proba)


class data_query(object):
    '''helper class to provide const values and access to fields created by sio.loadmat'''
    class ConstError(TypeError): pass
    def __setattr__(self,name,value):
        raise self.ConstError("Can modify element (%s)"%name)
    def __getattr__(self,name):
        return self.__dict__[name]
    def __init__(self, data):
        if type(data)==np.ndarray:
            for key in data.dtype.names:
                self.__dict__[key] = data[key].item()
        else:
            for key in data:
                self.__dict__[key] = data[key]
    def __getstate__(self):
        odict = self.__dict__.copy() # copy the dict since we change it
        #del odict['fh']              # remove filehandle entry
        return odict
    def __setstate__(self, dict):
        #fh = open(dict['file'])      # reopen file
        #count = dict['lineno']       # read from file...
        #while count:                 # until line count is restored
        #    fh.readline()
        #    count = count - 1
        self.__dict__.update(dict)   # update attributes
        #self.fh = fh                 # save the file object


class maze_query(const):
    '''calculaing and querying maze data'''
    def __init__(self, data, mode='hc5'):
        '''import data from given format'''
        if (mode=='hc5'):
            self.__init_hc5(data)
        elif (mode=='ubb'):
            self.__init_ubb(data)
        else:
            raise ValueError('Mode %s not implemented'%mode)
        self.__dict__['mode']=mode
            
        
    def __calc_boundaries(self,X,Y=None):
        incl = np.ones(len(X)).astype(bool)
        if Y is None:
            for point in self.exclude_points:
                incl = incl & ~(X == point[0])
            boundaries = [np.min(X[incl]), np.max(X[incl]), np.nan, np.nan]
        else:
            for point in self.exclude_points:
                incl = incl & ~((X == point[0]) & (Y == point[1]))
            boundaries = [np.min(X[incl]),
                np.max(X[incl]), np.min(Y[incl]), np.max(Y[incl])]
        return boundaries, incl
    
    
    def __init_hc5(self, datafile):
        '''import hc-5 data from dict in format provided by sio.loadmat'''
        import scipy.io as sio
        if type(datafile) == str:
            data = sio.loadmat(datafile, squeeze_me=True, struct_as_record=True)
        else:
            data = datafile
        for key in ['Par', 'Clu', 'Laps', 'Spike', 'Track']:
            self.__dict__[key] = data_query(data[key])
        # Data cleaning
        self.__dict__['exclude_points'] = [np.array([-1,-1]), np.array([0,0])]
        # Track
        b, v = self.__calc_boundaries(self.Track.X, self.Track.Y)
        self.Track.__dict__['Valid'] = v
        boundaries = { 'CalcMinX' : b[0], 'CalcMinY' : b[2],
                       'CalcMaxX' : b[1], 'CalcMaxY' : b[3]}
        self.Par.__dict__.update(boundaries)
        maze_unit_track = self.__guess_unit_hc5([self.Track.X, self.Track.Y])
        # Spike
        b, v = self.__calc_boundaries(self.Spike.X, self.Spike.Y)
        self.Spike.__dict__['Valid'] = v
        maze_unit_spike = self.__guess_unit_hc5([self.Spike.X, self.Spike.Y])
        # Speed
        maze_unit_speed = self.__guess_unit_hc5(self.Track.speed)
        maze_unit_wheel = self.__guess_unit_hc5([self.Laps.WhlSpeedCW,self.Laps.WhlSpeedCCW])
        # Meta
        if maze_unit_spike != maze_unit_track:
            raise ValueError("Mismatched length units")
        if maze_unit_speed != maze_unit_wheel:
            import warnings
            warnings.warn("Mismatched speed units")
        self.__dict__['unit'] = maze_unit_spike
        self.__dict__['maze'] = hc5

    def __init_ubb(self, datafile):
        '''import ubb data from dict in format'''
        import os.path
        dir, name = os.path.split(datafile)
        freq = 100
        nlaps = 12
        laplength = 600 / nlaps * freq
        nclust = 100
        enterleft = np.column_stack((np.arange(0,(nlaps-1)*laplength,laplength),
            np.arange(laplength-1,nlaps*laplength-1,laplength))).astype(int)
        self.__dict__['Par'] = data_query({ 'RatName': name, 'SamplingFrequency': freq,
            'nTimebins': nlaps * laplength, 'SyncOn': 0, 'SyncOff': nlaps * laplength,
            'NLap': nlaps, 'NMazeSect': 1, 'MazeSectEnterLeft': enterleft[:,np.newaxis,:],
            'totNch':nclust, 'NchanPerShank':1, 'BehavType':np.ones((nlaps,)).astype(int),
            'TrialType':np.zeros((nlaps,)).astype(int)})
        self.__dict__['Clu'] = data_query({ 'shank':np.array(range(1,1+nclust)),
            'localClu':np.ones((nclust,)).astype(int), 'totClu':np.array(range(1,1+nclust)),
            'isIntern':np.zeros((nclust,)).astype(int)})
        # Position
        trk = np.loadtxt(datafile.replace('spike','pos'),ndmin=2)
        dim = trk.shape[1]
        vol = np.ceil(np.power(np.max(trk),dim))
        self.__dict__['dim'] = dim
        self.__dict__['vol'] = vol
        x, y = trk[:,0], (trk[:,1] if dim>1 else np.zeros_like(trk[:,0]))
        s = np.linalg.norm(trk[1:,:]-trk[:-1,:],axis=1,ord=2)
        s = np.concatenate((s,[s[-1]])) * freq
        # eeg
        try:
            eeg = np.loadtxt(datafile.replace('spike','state'),ndmin=1)
        except:
            eeg = np.zeros_like(x)
        self.__dict__['Track'] = data_query({'X': x, 'Y': y, 'Z': trk,
            'linX': x, 'linY': y, 'linZ': trk, 'linDist': x,
            'speed': s, 'eeg': eeg})
        self.__dict__['Laps'] = data_query({ 'taskType': 'free_run',
            'WhlSpeedCW': np.zeros_like(x), 'WhlSpeedCCW': np.zeros_like(x)})
        # Spikes
        spk = np.loadtxt(datafile,skiprows=1,ndmin=2)
        tim = (spk[:,0] * freq).astype(int)-1
        self.__dict__['Spike'] = data_query({'res': tim, 'totclu': spk[:,1],
            'X': x[tim], 'Y': y[tim], 'Z': trk[tim],
            'linX': x[tim], 'linY': y[tim], 'linZ': trk[tim], 'linDist': x[tim],
            'speed': s[tim]})
        maze_unit_spike = 'm'
        maze_unit_track = 'm'
        if maze_unit_spike != maze_unit_track:
            raise ValueError("Mismatched units")
        self.__dict__['exclude_points'] = [np.array([np.nan, np.nan])]
        b, v = self.__calc_boundaries(self.Track.X, self.Track.Y)
        self.Track.__dict__['Valid'] = True
        boundaries = { 'CalcMinX' : b[0], 'CalcMinY' : b[2],
                       'CalcMaxX' : b[1], 'CalcMaxY' : b[3]}
        self.Par.__dict__.update(boundaries)
        self.Spike.__dict__['Valid'] = True
        self.__dict__['unit'] = maze_unit_spike
        self.__dict__['maze'] = ubb
        
    def __guess_unit_hc5(self, lengths):
        '''guess unit of length knowing that the maze is 1200mm long'''
        if np.max(lengths)>1000:
            return 'mm'
        else:
            return '0.1in'

    def __guess_speed_hc5(self, lengths):
        '''guess unit of length knowing that the maze is 1200mm long'''
        if np.max(lengths)>500:
            return 'mm'
        else:
            return '0.1in'
    
    def propose_bins(self, widthmm, dim=None):
        if dim is None:
            if self.mode == 'hc5':
                return np.arange(0,self.mm_to_len(1200+widthmm),self.mm_to_len(widthmm))
            elif self.mode == 'ubb':
                length = np.power(self.vol,1.0/self.dim)
                return np.arange(0,length+self.mm_to_len(widthmm),self.mm_to_len(widthmm))
            else:
                pass
        elif dim==1 :
            if self.mode == 'hc5':
                return np.arange(0,self.Par.CalcMaxprojDist+self.mm_to_len(widthmm),self.mm_to_len(widthmm))
            elif self.mode == 'ubb':
                length = np.power(self.vol,1.0/self.dim)
                return (np.arange(0,length+self.mm_to_len(widthmm),self.mm_to_len(widthmm)))
            else:
                pass
        else:
            raise ValueError('Custom dim not implemented yet.')
        
    def tic_to_time(self, arr):
        '''convert clock tics to sec'''
        freq = self.Par.SamplingFrequency
        return arr * (1.0 / freq)
        
    def time_to_tic(self, arr):
        '''convert sec to clock tics'''
        freq = self.Par.SamplingFrequency
        return np.floor(arr * freq)
        
    def len_to_mm(self, arr):
        '''convert inherent length to mm'''
        if self.unit=='mm':
            return arr
        elif self.unit=='0.1in':
            return arr*2.54
        elif self.unit=='m':
            return arr*1000
        else:
            raise ValueError("Unknown unit")
        
    def mm_to_len(self, arr):
        '''convert inherent length to mm'''
        if self.unit=='mm':
            return arr
        elif self.unit=='0.1in':
            return arr/2.54
        elif self.unit=='m':
            return arr/1000.0
        else:
            raise ValueError("Unknown unit")
    
    def get_section_ids(self, name):
        '''convert section name to section id in range(0,maze.Section_count)'''
        if (name == 'all'):
            sel = np.array(range(0,self.maze.Section_count))
        else:
            sel = self. maze.Section_tx[1:] == name
        if sel.dtype is np.dtype(bool) and ~sel.any():
            raise ValueError('Choose section id from %s'%self.maze.Section_tx)
        return sel
        
    def get_intervals(self, section_in, section_out):
        '''(nlaps) x |{start, end}| clock signals for given spatial course in each lap'''
        laps_events = self.Par.MazeSectEnterLeft
        section_id_in = self.get_section_ids(section_in) if (type(section_in)==str) else section_in
        section_id_out = self.get_section_ids(section_out) if (type(section_out)==str) else section_out
        ret = []
        for tics in laps_events:
            tic_in = tics[section_id_in]
            tic_in = min(tic_in[tic_in>0])
            tic_out = tics[section_id_out]
            tic_out = max(tic_out[tic_out>0])
            ret.append([tic_in,tic_out])
        return np.array(ret)

    def get_tics(self, section):
        "(nlaps) x (1) clock signals for given maze section in each lap"
        laps_events = self.Par.MazeSectEnterLeft
        section_id = self.get_section_ids(section) if (type(section)==str) else section
        ret = []
        for tics in laps_events:
            ret.append(tics[section_id])
        return np.array(ret)

    #def get_tics(self, section_in=None, section_out=None):
    #    if (section_out is None):
    #        return self.__get_tics1(section_in)
    #    else:
    #        return self.__get_tics2(section_in, section_out)

    def get_spikes(self, tic_start, tic_end):
        '''return spikes in format (clock signal, cell id) for given clock interval'''
        spk_tic = self.Spike.res
        spk_cluster = self.Spike.totclu
        sel = (spk_tic>=tic_start) & (spk_tic<tic_end)
        return spk_tic[sel], spk_cluster[sel]
    
    def get_spike_meta(self, tic_start, tic_end, field):
        '''return spike meta field for given clock interval'''
        spk_tic = self.Spike.res
        data = getattr(self.Spike,field)
        sel = (spk_tic>=tic_start) & (spk_tic<tic_end)
        return data[sel]
    
    def get_shifted_spikes(self, tic_start, tic_end):
        '''return spikes in format (clock signal-clock start, cell id) for given clock interval'''
        spk_tic, spk_cluster = self.get_spikes(tic_start,tic_end)
        return spk_tic-tic_start, spk_cluster
    
    def get_spike_trains(self, tic_start, tic_end, clusters = None):
        '''return fired (cell id) x (clock) as bool for given clock interval'''
        spk_tic, spk_cluster = self.get_shifted_spikes(tic_start,tic_end)
        if clusters is None:
            clusters = self.Clu.totClu
        spk_train = np.zeros((len(clusters),tic_end-tic_start))
        nclusters = len(clusters)
        for c in range(0,nclusters):
            sel = (spk_cluster == clusters[c])
            spk_train[c,spk_tic[sel]] = 1
        return spk_train
        
    def __bin_it(self, data, width = None, normalize = 'none'):
        '''calculate averages (field) x (time bin) for given bin width'''
        method = {'none': (lambda x: 1), 'tics': (lambda x: 1.0/x),
                  'time': (lambda x: 1.0/self.tic_to_time(x))}
        normalizer = method[normalize]
        nfields, ntics = data.shape
        if (width is None) or (ntics<=width):
            boundaries = np.array([0.0, ntics])
        else:
            boundaries = np.arange(0.0,ntics,1.0*width)
        nbins = len(boundaries)-1
        binned = np.zeros((nfields,nbins))
        for b in range(0,nbins):
            start, end = np.floor(boundaries[[b,b+1]])
            factor = normalizer(end-start)
            binned[:,b] = np.sum(data[:,start:end],axis=1) * factor
        #centers = (boundaries[1:]+boundaries[:-1])/2
        return boundaries, binned
        
    def get_spike_counts(self, tic_start, tic_end, width = None, smooth_tic = 0, clusters = None):
        '''count spikes in (cell id) x (time bin) as int for given clock interval,
           optionally the time series is smoothed with a sigma=smooth_tic gaussian filter
        '''
        spk_train = self.get_spike_trains(tic_start, tic_end, clusters)
        if (smooth_tic>0):
            from scipy.ndimage.filters import gaussian_filter1d
            spk_train = gaussian_filter1d(spk_train.astype(float),smooth_tic,axis=1)
        bin_tic, spk_count = self.__bin_it(spk_train, width, normalize='none')
        #bin_tim = self.tic_to_time(tic_start + bin_tic)
        return tic_start + bin_tic, spk_count

    def __get_field(self, tic_start, tic_end, field):
        '''return data field for given clock interval'''
        if field[0:3] == 'Whl':
            key = 'Laps'
        else:
            key = 'Track'
        data = getattr(getattr(self,key),field)
        return data[tic_start:tic_end]

    def get_fields(self, tic_start, tic_end, fields):
        '''return data field for given clock interval'''
        if type(fields)==str:
            fields = [fields]
        nfields = len(fields)
        ntics = tic_end-tic_start
        rawdata = np.zeros((nfields,ntics))
        for f in range(0,nfields):
            rawdata[f,:] = self.__get_field(tic_start,tic_end,fields[f])
        return rawdata
    
    def get_averages(self, tic_start, tic_end, width = None, fields = []):
        '''bin values for given clock interval'''
        rawdata = self.get_fields(tic_start,tic_end,fields)
        averages = self.__bin_it(rawdata, width, normalize='tics')
        return averages
 
 
    ###   P L A C E   F I E L D   M A N I P U L A T I O N   ###
 
    def __get_2d_place_fields(self, spatial_bin, scheme='%s'):
        '''bin values for given cluster'''
        # get lap info
        laps = self.get_intervals('all','all')
        nlaps = len(laps)
        clusters = self.Clu.totClu
        nclust = len(clusters)
        freq = 1.0 * self.Par.SamplingFrequency
        histo, sepx, sepy = np.histogram2d([], [], spatial_bin)
        lenx, leny = len(sepx)-1, len(sepy)-1
        # initialize
        avgr = np.zeros((nclust,nlaps))
        rate = np.empty((nlaps,nclust,lenx,leny))
        info = np.empty((nlaps,nclust,lenx,leny))
        p_x_ = np.empty((nlaps,lenx,leny))
        # process data per lap
        for l in range(0,nlaps):
            lap_start, lap_end, numtics = laps[l,0], laps[l,1], laps[l,1]-laps[l,0]
            posX = self.__get_field(lap_start, lap_end, scheme%'X')
            posY = self.__get_field(lap_start, lap_end, scheme%'Y')
            pos_histo, sepx, sepy = np.histogram2d(posX, posY, spatial_bin)
            spk_tic, spk_clust = self.get_spikes(lap_start, lap_end)
            spkX = self.get_spike_meta(lap_start, lap_end, scheme%'X')
            spkY = self.get_spike_meta(lap_start, lap_end, scheme%'Y')
            for c in range(0,nclust):
                sel = spk_clust == clusters[c]
                spk_histo, sepx, sepy = np.histogram2d(spkX[sel], spkY[sel], spatial_bin)
                avgr[c,l] = sum(sel) * freq / numtics
                rate[l,c,:,:] = np.divide(spk_histo,pos_histo) * freq
                info[l,c,:,:] = np.multiply(spk_histo,(np.log(rate[l,c])-
                                                   np.log(sum(sel)/numtics)) ) * freq
            p_x_[l,:,:] = pos_histo
        return avgr, rate, info, p_x_
        
    def __get_2d_place_fields_parallel(self, spatial_bin, scheme='%s'):
        '''bin values for given cluster, note: serializing data might impact performance'''
        
        #from joblib import Parallel, delayed
        import multiprocessing as mp
        import time
        import sys

        num_cores = mp.cpu_count()
        # get lap info
        laps = self.get_intervals('all','all')
        nlaps = len(laps)
        clusters = self.Clu.totClu
        nclust = len(clusters)
        freq = 1.0 * self.Par.SamplingFrequency
        histo, sepx, sepy = np.histogram2d([], [], spatial_bin)
        lenx, leny = len(sepx)-1, len(sepy)-1
        # initialize
        avgr = np.zeros((nclust,nlaps))
        rate = np.empty((nlaps,nclust,lenx,leny))
        info = np.empty((nlaps,nclust,lenx,leny))
        p_x_ = np.empty((nlaps,lenx,leny))
        # process data per lap
        jobs = []
        print ('start',time.clock())
        sys.stdout.flush()
        for l in range(0,nlaps):
            lap_start, lap_end, numtics = laps[l,0], laps[l,1], laps[l,1]-laps[l,0]
            posX = self.__get_field(lap_start, lap_end, scheme%'X')
            posY = self.__get_field(lap_start, lap_end, scheme%'Y')
            spk_tic, spk_clust = self.get_spikes(lap_start, lap_end)
            spkX = self.get_spike_meta(lap_start, lap_end, scheme%'X')
            spkY = self.get_spike_meta(lap_start, lap_end, scheme%'Y')
            #jobs.append(delayed(_get_2d_field)(
            #            numtics, posX, posY, spkX, spkY, spk_clust, nclust, spatial_bin, freq))
            jobs.append([numtics, posX, posY, spkX, spkY, spk_clust, clusters, spatial_bin, freq])
        print ('prep\'d',time.clock())
        sys.stdout.flush()
        #ret = Parallel(n_jobs=int(num_cores/2), verbose=100)(jobs)
        #OR ret = Parallel(n_jobs=-2, verbose=100)(jobs)
        p = mp.Pool(max(1,int(num_cores-1)));
        ret = p.map(_get_2d_field_mapper, jobs)
        print ('combine',time.clock())
        sys.stdout.flush()
        #ret = Parallel(n_jobs=num_cores)(jobs)
        for l in range(0,nlaps):
            a, r, i, h = ret[l]
            avgr[:,l] = a
            rate[l,:,:,:] = r
            info[l,:,:,:] = i
            p_x_[l,:,:] = h
        print ('fin',time.clock())
        return avgr, rate, info, p_x_
   
    def get_2d_place_field(self, spatial_bin, scheme='%s', parallel=True):
        '''do place field statitics for given cluster'''
        if parallel:
            avgr, rate, info, p_x_ = self.__get_2d_place_fields_parallel(spatial_bin, scheme)
        else:
            avgr, rate, info, p_x_ = self.__get_2d_place_fields(spatial_bin, scheme)
        # prepare output
        stat = PlaceField()
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=RuntimeWarning)
            stat.rate = avgr
            stat.mean = np.nanmean(rate, axis=0)
            stat.median = np.nanmedian(rate, axis=0)
            stat.std = np.nanstd(rate, axis=0)
            med = stat.median
            # one extra spike per lap does not necesserily indicate outlier (for zero median)
            stat.outlier = np.sum(np.less(rate, med - np.sqrt(med)) |
                                  np.greater(rate - 1./p_x_[:,np.newaxis,:,:], med + np.sqrt(med)),
                                  axis=0)
            stat.info_mean = np.nanmean(info, axis=0)
            stat.info_std = np.nanstd(info, axis=0)
            stat.p_x_ = np.nansum(p_x_, axis=0)
        return stat

    def __get_1d_place_fields(self, spatial_bin, scheme='%s'):
        '''bin values for given cluster'''
        # get lap info
        laps = self.get_intervals('all','all')
        nlaps = len(laps)
        clusters = self.Clu.totClu
        nclust = len(clusters)
        freq = 1.0 * self.Par.SamplingFrequency
        histo, sepx = np.histogram([], spatial_bin)
        lenx = len(sepx)-1
        # initialize
        avgr = np.zeros((nclust,nlaps))
        rate = np.empty((nlaps,nclust,lenx))
        info = np.empty((nlaps,nclust,lenx))
        p_x_ = np.empty((nlaps,lenx))
        # process data per lap
        for l in range(0,nlaps):
            lap_start, lap_end, numtics = laps[l,0], laps[l,1], laps[l,1]-laps[l,0]
            posX = self.__get_field(lap_start, lap_end, scheme%'Dist')
            pos_histo, sepx = np.histogram(posX, spatial_bin)
            spk_tic, spk_clust = self.get_spikes(lap_start, lap_end)
            spkX = self.get_spike_meta(lap_start, lap_end, scheme%'Dist')
            for c in range(0,nclust):
                sel = spk_clust == clusters[c]
                spk_histo, sepx = np.histogram2d(spkX[sel], spatial_bin)
                avgr[c,l] = sum(sel) * freq / numtics
                rate[l,c,:] = np.divide(spk_histo,pos_histo) * freq
                info[l,c,:] = np.multiply(spk_histo,(np.log(rate[l,c])-
                                                   np.log(sum(sel)/numtics)) ) * freq
            p_x_[l,:] = pos_histo
        return avgr, rate, info, p_x_
        
    def __get_1d_place_fields_parallel(self, spatial_bin, scheme='%s'):
        '''bin values for given cluster, note: serializing data might impact performance'''
        
        #from joblib import Parallel, delayed
        import multiprocessing as mp
        import time
        import sys

        num_cores = mp.cpu_count()
        # get lap info
        laps = self.get_intervals('all','all')
        nlaps = len(laps)
        clusters = self.Clu.totClu
        nclust = len(clusters)
        freq = 1.0 * self.Par.SamplingFrequency
        histo, sepx = np.histogram([], spatial_bin)
        lenx = len(sepx)-1
        # initialize
        avgr = np.zeros((nclust,nlaps))
        rate = np.empty((nlaps,nclust,lenx))
        info = np.empty((nlaps,nclust,lenx))
        p_x_ = np.empty((nlaps,lenx))
        # process data per lap
        jobs = []
        print ('start',time.clock())
        sys.stdout.flush()
        for l in range(0,nlaps):
            lap_start, lap_end, numtics = laps[l,0], laps[l,1], laps[l,1]-laps[l,0]
            posX = self.__get_field(lap_start, lap_end, scheme%'Dist')
            spk_tic, spk_clust = self.get_spikes(lap_start, lap_end)
            spkX = self.get_spike_meta(lap_start, lap_end, scheme%'Dist')
            #jobs.append(delayed(_get_1d_field)(
            #            numtics, posX, spkX, spk_clust, nclust, spatial_bin, freq))
            jobs.append([numtics, posX, spkX, spk_clust, clusters, spatial_bin, freq])
        print ('prep\'d',time.clock())
        sys.stdout.flush()
        #ret = Parallel(n_jobs=int(num_cores/2), verbose=100)(jobs)
        #OR ret = Parallel(n_jobs=-2, verbose=100)(jobs)
        p = mp.Pool(max(1,int(num_cores-1)));
        ret = p.map(_get_1d_field_mapper, jobs)
        print ('combine',time.clock())
        sys.stdout.flush()
        #ret = Parallel(n_jobs=num_cores)(jobs)
        for l in range(0,nlaps):
            a, r, i, h = ret[l]
            avgr[:,l] = a
            rate[l,:,:] = r
            info[l,:,:] = i
            p_x_[l,:] = h
        print ('fin',time.clock())
        return avgr, rate, info, p_x_
    
    def get_1d_place_field(self, spatial_bin, scheme='lin', parallel=True):
        '''do place field statitics for given cluster'''
        if parallel:
            avgr, rate, info, p_x_ = self.__get_1d_place_fields_parallel(spatial_bin, scheme)
        else:
            avgr, rate, info, p_x_ = self.__get_1d_place_fields(spatial_bin, scheme)
        # prepare output
        stat = PlaceField()
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=RuntimeWarning)
            stat.rate = avgr
            stat.mean = np.nanmean(rate, axis=0)
            stat.median = np.nanmedian(rate, axis=0)
            stat.std = np.nanstd(rate, axis=0)
            med = stat.median
            # one extra spike per lap does not necesserily indicate outlier (for zero median)
            stat.outlier = np.sum(np.less(rate, med - np.sqrt(med)) |
                                  np.greater(rate - 1./p_x_[:,np.newaxis,:], med + np.sqrt(med)),
                                  axis=0)
            stat.info_mean = np.nanmean(info, axis=0)
            stat.info_std = np.nanstd(info, axis=0)
            stat.p_x_ = np.nansum(p_x_, axis=0)
            print("check", p_x_.shape, stat.p_x_.shape)
        return stat


       
    ###  P L O T T I N G   ###  
    def plot_track(self,ax,start,end,smoothen=0):
        '''plot the animals position in the maze and its speed'''
        from scipy.ndimage.filters import gaussian_filter
        tim = self.tic_to_time(np.array(range(start,end)))
        pos = self.get_fields(start,end,['X','Y'])
        spe = self.get_fields(start,end,['speed'])
        if smoothen:
            spe = gaussian_filter(spe.astype(float),self.time_to_tic(smoothen))
        ax.set_xlabel("time (s)")
        ax.plot(tim,self.len_to_mm(pos).T,'k-',label='position')
        ax.plot(tim,self.len_to_mm(spe).T,'b-',label='speed',zorder=-1)
        ax.set_ylabel("position, speed (mm, mm/s)")

    def plot_activity(self,ax,start,end,smoothen=0):
        '''plot wheel speed and eeg'''
        from scipy.ndimage.filters import gaussian_filter
        tim = self.tic_to_time(np.array(range(start,end)))
        whl = np.abs(self.get_fields(start,end,['WhlSpeedCW']))-np.abs(self.get_fields(start,end,['WhlSpeedCCW']))
        if smoothen:
            whl = gaussian_filter(whl.astype(float),self.time_to_tic(smoothen))
        eeg = self.get_fields(start,end,['eeg'])
        eeg = (eeg-np.nanmin(eeg))/(np.nanmax(eeg)-np.nanmin(eeg))
        ax.set_xlabel("time (s)")
        ax.plot(tim,whl.T,'r-',label='wheel')
        ax.plot(tim,800+eeg.T*100,'g-',label='eeg',zorder=-1)
        
    def plot_spikes(self,ax,start,end):
        '''plot spikes with hlines'''
        from matplotlib.pyplot import cm
        totclu = self.Clu.totClu
        locclu = self.Clu.localClu
        ncol = max(max(locclu),10)
        color=cm.rainbow(np.linspace(0,1,ncol))
        spk_tic, spk_clu = self.get_spikes(start,end)
        spk_tim = self.tic_to_time(spk_tic)
        for clu in totclu:
            idx = clu-1
            sel = spk_clu == clu
            pin = locclu[idx]-1
            ax.plot(spk_tim[sel],spk_clu[sel],'|',markersize=4,color=color[pin])
        ax.set_xlabel("time (s)")
        ax.set_ylabel("cluster ID")
        
    def plot_counts(self,ax,start,end,width):
        '''show spike counts as a heatmap'''
        totclu = self.Clu.totClu
        spk_tic, spk_count = self.get_spike_counts(start,end,width)
        spk_tim = self.tic_to_time(spk_tic)
        ext = spk_tim[0], spk_tim[-1], min(totclu), max(totclu)
        ax.matshow(spk_count,extent=ext,origin='lower',aspect='auto')
        ax.xaxis.set_ticks_position('bottom')
        #ax.set_ylim([0,max(totclu)])
        ax.set_xlabel("time (s)")
        ax.set_ylabel("cluster ID")
        
    def visualize_spikes(self,start,end,figsize=(16,12)):
        '''plot above, spikes with hlines'''
        import matplotlib.pyplot as plt
        fig, ax1 = plt.subplots(figsize=figsize)
        ax2 = ax1.twinx()
        self.plot_spikes(ax1,start,end)
        self.plot_track(ax2,start,end,0.05)
        self.plot_activity(ax2,start,end)
        plt.legend(loc='upper right')
        fig.show()
        
    def visualize_counts(self,start,end,width,figsize=(16,12)):
        '''show above, spike counts as a heatmap'''
        import matplotlib.pyplot as plt
        fig, ax1 = plt.subplots(figsize=figsize)
        ax2 = ax1.twinx()
        self.plot_track(ax2,start,end,0.05)
        self.plot_activity(ax2,start,end)
        self.plot_counts(ax1,start,end,width)
        plt.legend(loc='upper right')
        fig.show()
        
    @staticmethod
    def seizmic_plots(self,X,Y,step=1,figsize=(16,12)):
        '''plot several channels with stepping, as common is seizmic diagrams'''
        import matplotlib.pyplot as plt
        fig, ax1 = plt.subplots(figsize=figsize)
        nDims = Y.shape[1]
        for i in range(0,nDims):
            ax1.plot(X,Y[:,i]+step*i)
        fig.show()
        
    
    ###  L I N E A R I Z I N G   T H E   T R A C K   ###
    @staticmethod
    #def __func_bifur(y, base, center, height, width):
    def __func_bifur(y, center, width, base, height):
        '''histogram of a fork'''
        return base + height*np.cos(2*np.pi*(y-center)/width)
    @staticmethod
    #def __func_plateau(x, base, x1, x2, height, width):
    def __func_plateau(x, x1, x2, width, base, height):
        '''histogram of small hill@x1 plateau small hill@(x2-2w) then big hill@x2'''
        return base + height*np.exp(-((x-x1)/width)**2) + (
            height*np.exp(-((x-x2+2*width)/width)**2) + 2*height*np.exp(-((x-x2)/width)**2) )
    
    @staticmethod
    def __fit_points(field, base=[None, None], debug=False):
        # projection of visited points
        histo_y=np.sum(field,axis=0)
        histo_x=np.sum(field,axis=1)
        # amplitude
        max_y, max_x = max(histo_y), max(histo_x)
        # grid axes
        base_y = np.array(range(0,len(histo_y))) if base[1] is None else base[1]
        base_x = np.array(range(0,len(histo_x))) if base[0] is None else base[0]
        # range definitions
        start_y, end_y = base_y[np.nonzero(histo_y)[0][[0,-1]]]
        start_x, end_x = base_x[np.nonzero(histo_x)[0][[0,-1]]]
        
        from scipy.optimize import curve_fit
        # triple equidistant hills in y
        pbifur, pcov = curve_fit(maze_query.__func_bifur, base_y, histo_y,
                                 bounds=([start_y, (end_y-start_y)/4, 0, max_y/4],
                                 [end_y, 3*(end_y-start_y)/4, max_y/2, 3*max_y/4]))
        # small hill, ___, small hill, big hill
        pplateau, pcov = curve_fit(maze_query.__func_plateau, base_x, histo_x,
                                 bounds=([start_x, (start_x+end_x)/2, start_x, 0, 0],
                                 [(start_x+end_x)/2, end_x, (3*start_x+end_x)/4, max_x/2, max_x/2]))
        
        if debug:
            import matplotlib.pyplot as plt
            f, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, sharex='col', sharey='row')
            ax1.matshow(field.T,extent=(base_x[0],base_x[-1],base_y[0],base_y[-1]),origin='lower')
            ax1.set_xlabel('X (default unit)')
            ax1.set_ylabel('Y (default unit)')
            ax1.xaxis.set_label_position('top')
            #ax1.xaxis.set_ticks_position('none') <-- not effective, imshow overrides to top
            print (pbifur)
            ax2.plot(histo_y,base_y)
            ax2.plot(maze_query.__func_bifur(base_y, *pbifur),base_y)
            ax2.plot(pbifur[2],pbifur[0],'rx')
            ax2.xaxis.set_ticks_position('top')
            ax2.set_xlabel('marginal histo(y)')
            ax2.xaxis.set_label_position('top')
            ax2.set_ylim(base_y[[0,-1]])
            print (pplateau)
            ax3.plot(base_x,histo_x)
            ax3.plot(base_x,maze_query.__func_plateau(base_x, *pplateau))
            #ax3.xaxis.set_ticks_position('none') <-- not effective, sharex overrides to bottom
            plt.setp(ax3.get_xticklabels(), visible=False)
            ax3.set_ylabel('marginal histo(x)')
            ax3.set_xlim(base_x[[0,-1]])
            ax4.axis('off')
            plt.show()

        return (pbifur[[0,1]], pplateau[[0,1]])

    
    def __fit_track_hc5(self, debug=False):
        div=self.propose_bins(25)
        posX = self.Track.X
        posY = self.Track.Y
        pos_histo, sepx, sepy = np.histogram2d(posX, posY, (div,div))
        pos_histo = np.array(pos_histo>0).astype(int)
        sepx = (sepx[1:]+sepx[:-1])/2
        sepy = (sepy[1:]+sepy[:-1])/2
        anchor_y, anchor_x = maze_query.__fit_points(pos_histo, [sepx, sepy], debug)
        self.__dict__['anchor_x'] = anchor_x
        self.__dict__['anchor_y'] = anchor_y
        
    def fit_track(self, debug=False):
        if (self.maze.name=='hc5'):
            self.__fit_track_hc5(debug)
        else:
            print ("No need to fit track")
            self.__dict__['anchor_x'] = np.array([0])
            self.__dict__['anchor_y'] = np.array([0])
        
    @staticmethod
    def __corner_proj(x, y, r=1, cx=0, cy=0, ord=None):
        c = np.column_stack((cx, cy))
        v = np.column_stack((x,y))-c
        l = np.linalg.norm(v,axis=1,ord=ord)
        p = c + r*np.divide(v,np.column_stack((l,l)))
        # TODO: transposing all vectors makes the formula for division simpler
        return p
    @staticmethod
    def __side_proj(x, r=1, cx=0):
        side = np.sign(x-cx)
        p = cx+side*r
        return p
    @staticmethod
    def __side_projx(x, y, r=1, cx=0, cy=None):
        p = np.column_stack((maze_query.__side_proj(x,r,cx), y))
        return p
    @staticmethod
    def __side_projy(x, y, r=1, cx=None, cy=0):
        p = np.column_stack((x, maze_query.__side_proj(y,r,cy)))
        return p
    @staticmethod
    def __nasty_proj(x, y, r=1, cx=0, cy=0):
        return maze_query.__corner_proj(x, y, r=r, cx=cx, cy=cy, ord=np.inf)
    @staticmethod
    def __rank(x, bins):
        n = len(bins)
        r = np.zeros(np.array(x).shape).astype(int)
        for i in range(0,n):
            r += np.array(x>bins[i]).astype(int)
        return r
    @staticmethod
    def __unique(seq):
        seen = set()
        seen_add = seen.add
        return [x for x in seq if not (x in seen or seen_add(x))]

    def __draw_track_hc5(self,res=1000):
        width = self.anchor_y[1]/2
        pivot_y = self.anchor_y[[0]] + np.array([2, 0, -2])*width
        pivot_x = self.anchor_x + np.array([0, 0])*width
        trk = []
        #  geometric convention (origin in the lower left corner
        #  /----U--->\      P----------P
        #  UL   1    UR     |          |
        #  +<---M----+      P----------P
        #  LL   2    LR     |          |
        #  \----L--->/      P----------P
        # mid
        trk.append((np.linspace(pivot_x[1],pivot_x[0],res),np.linspace(pivot_y[1],pivot_y[1],res)))
        # upper left
        trk.append((np.linspace(pivot_x[0],pivot_x[0],res),np.linspace(pivot_y[1],pivot_y[2],res)))
        # upper
        trk.append((np.linspace(pivot_x[0],pivot_x[1],res),np.linspace(pivot_y[2],pivot_y[2],res)))
        # upper right
        trk.append((np.linspace(pivot_x[1],pivot_x[1],res),np.linspace(pivot_y[1],pivot_y[2],res)))
        # lower left
        trk.append((np.linspace(pivot_x[0],pivot_x[0],res),np.linspace(pivot_y[1],pivot_y[0],res)))
        # lower
        trk.append((np.linspace(pivot_x[0],pivot_x[1],res),np.linspace(pivot_y[0],pivot_y[0],res)))
        # lower right
        trk.append((np.linspace(pivot_x[1],pivot_x[1],res),np.linspace(pivot_y[1],pivot_y[0],res)))
        return np.hstack(trk).T
        
    def __project_bins_hc5(self,spatial_bin,debug=False):
        trk = self.__draw_track()
        bin_x, bin_y = spatial_bin
        binned = zip(np.digitize(trk[:,0], bin_x), np.digitize(trk[:,1], bin_y))
        return np.array(maze_query.__unique(binned))-1
       

    def __project_track(self,x,y,valid,debug=False):
        width = self.anchor_y[1]/2
        pivot_y = self.anchor_y[[0]]
        cent_y = self.anchor_y[0] + np.array([-1, 1])*width
        pivot_x = self.anchor_x + np.array([1, -1])*width
        cent_x = self.anchor_x + np.array([1, -1])*width
        sec_x, sec_y = maze_query.__rank(x, pivot_x), maze_query.__rank(y, pivot_y)
        proj = np.column_stack((x,y))
        #  geometric convention (origin in the lower left corner
        #  /----------\      /----------\
        #  |UL  UM  UR|      | C      C |
        #  +----------+      +-P------P-+
        #  |LL  LM  LR|      | C      C |
        #  \----------/      \----------/
        # lower left
        sel = (sec_x==0) & (sec_y==0)
        proj[sel] = maze_query.__nasty_proj(x[sel],y[sel],width,cent_x[0],cent_y[0])
        # upper left
        sel = (sec_x==0) & (sec_y==1)
        proj[sel] = maze_query.__nasty_proj(x[sel],y[sel],width,cent_x[0],cent_y[1])
        # lower mid
        sel = (sec_x==1) & (sec_y==0)
        proj[sel] = maze_query.__side_projy(x[sel],y[sel],width,None,cent_y[0])
        # upper mid
        sel = (sec_x==1) & (sec_y==1)
        proj[sel] = maze_query.__side_projy(x[sel],y[sel],width,None,cent_y[1])
        # lower right
        sel = (sec_x==2) & (sec_y==0)
        proj[sel] = maze_query.__nasty_proj(x[sel],y[sel],width,cent_x[1],cent_y[0])
        # upper right
        sel = (sec_x==2) & (sec_y==1)
        proj[sel] = maze_query.__nasty_proj(x[sel],y[sel],width,cent_x[1],cent_y[1])
        proj[~valid] = self.exclude_points[0]

        if debug:
            import matplotlib.pyplot as plt
            plt.plot(x,y,',')
            plt.plot(proj[:,0],proj[:,1],'.')
            print(pivot_x,pivot_y)
            ps = cartesian(pivot_x,pivot_y)
            plt.plot(ps[:,0],ps[:,1],'x',label='pivot')
            ps = cartesian(cent_x,cent_y)
            plt.plot(ps[:,0],ps[:,1],'x',label='center')
            plt.legend(loc='upper right')
            plt.xlabel("X (default unit)")
            plt.ylabel("Y (default unit)")
        return proj
    

    def __linearize_track(self,x,y,valid,debug=False):
        width = self.anchor_y[1]/2
        pivot_y = self.anchor_y[[0]] + np.array([1, -1])*1.0
        pivot_x = self.anchor_x + np.array([1, -1])*1.0
        sec_x, sec_y = maze_query.__rank(x, pivot_x), maze_query.__rank(y, pivot_y)
        proj = np.ones_like(x)
        #  geometric convention (origin in the lower left corner
        # /----U-----\      /----------\
        # UL        UR      |          |
        # +----M-----+      +PP------PP+
        # LL        LR      |          |
        # \----L-----/      \----------/
        # mid
        dist = 0
        sel = (sec_x==1) & (sec_y==1)
        proj[sel] = dist + np.abs(x[sel]-pivot_x[1])
        # upper left
        dist += pivot_x[1]-pivot_x[0]
        sel = (sec_x==0) & (sec_y==2)
        proj[sel] = dist + np.abs(y[sel]-self.anchor_y[[0]])
        # upper
        dist += 2*width
        sel = (sec_x==1) & (sec_y==2)
        proj[sel] = dist + np.abs(x[sel]-pivot_x[0])
        # upper right
        dist += pivot_x[1]-pivot_x[0]
        sel = (sec_x==2) & (sec_y==2)
        proj[sel] = dist + np.abs(y[sel]-(self.anchor_y[[0]]+2*width))
        # lower left
        dist += 2*width
        sel = (sec_x==0) & (sec_y==0)
        proj[sel] = dist + np.abs(y[sel]-self.anchor_y[[0]])
        # lower
        dist += 2*width
        sel = (sec_x==1) & (sec_y==0)
        proj[sel] = dist + np.abs(x[sel]-pivot_x[0])
        # lower right
        dist += pivot_x[1]-pivot_x[0]
        sel = (sec_x==2) & (sec_y==0)
        proj[sel] = dist + np.abs(y[sel]-(self.anchor_y[[0]]-2*width))
        proj[~valid] = self.exclude_points[0][0]

        if debug:
            import matplotlib.pyplot as plt
            f, (ax1, ax2) = plt.subplots(1,2,sharey='row')
            ax1.plot(x,proj,'.')
            ax1.set_ylabel("Linearized position (default unit)")
            ax1.set_xlabel("X (default unit)")
            ax2.plot(y,proj,'.')
            ax2.set_xlabel("Y (default unit)")
            print(pivot_x,pivot_y)
        return proj
    
    def __linearize_tracks_hc5(self,debug=False):
        proj = self.__project_track(self.Track.X,self.Track.Y,self.Track.Valid,debug)
        self.Track.__dict__['projX'] = proj[:,0]
        self.Track.__dict__['projY'] = proj[:,1]
        proj = self.__linearize_track(proj[:,0],proj[:,1],self.Track.Valid,debug)
        self.Track.__dict__['projDist'] = proj[:]
        proj = self.__project_track(self.Spike.X,self.Spike.Y,self.Spike.Valid,False)
        self.Spike.__dict__['projX'] = proj[:,0]
        self.Spike.__dict__['projY'] = proj[:,1]
        proj = self.__linearize_track(proj[:,0],proj[:,1],self.Spike.Valid,False)
        self.Spike.__dict__['projDist'] = proj[:]

    def __linearize_tracks_ubb(self,debug=False):
        self.Track.__dict__['projX'] = self.Track.linX
        self.Track.__dict__['projY'] = self.Track.linY
        self.Track.__dict__['projZ'] = self.Track.linZ
        self.Track.__dict__['projDist'] = self.Track.linDist
        self.Spike.__dict__['projX'] = self.Spike.linX
        self.Spike.__dict__['projY'] = self.Spike.linY
        self.Spike.__dict__['projZ'] = self.Spike.linZ
        self.Spike.__dict__['projDist'] = self.Spike.linDist
        
    def linearize_tracks(self,debug=False):
        if (self.maze.name=='hc5'):
            self.__linearize_tracks_hc5(debug)
        else:
            self.__linearize_tracks_ubb(debug)
            print ("No need to linearize track")
        b, v = self.__calc_boundaries(self.Track.projDist)
        boundaries = { 'CalcMinprojDist' : b[0], 'CalcMaxprojDist' : b[1]}
        self.Par.__dict__.update(boundaries)
