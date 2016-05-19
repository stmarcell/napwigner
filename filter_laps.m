function keep = filter_laps(D,varargin)
%Filter the laps in thel struct D with firing rate count below 1 s.d of
%       the mean of all the laps.
%
%Ruben Pinzon @2015
type = 'spike_cnt';
if nargin > 1
    type = varargin{1};
end

if strcmp(type,'spike_cnt')
    n_trials    = length(D);
    cnt_total   = zeros(1, n_trials);
    for n = 1 : n_trials    
        cnt_total(n) = sum(sum(D(n).y))/D(n).T;        
    end

    %trials with low total firing count
    mu = mean(cnt_total);
    sd = std(cnt_total);

    keep = cnt_total >= (mu-sd); 
    
else 
    T = [D.T];
    keep = T >= (mean(T)-std(T));
end

fprintf('%d trials filter out: [%s]\n', sum(~keep), sprintf('%d ',[D(~keep).trialId]))
