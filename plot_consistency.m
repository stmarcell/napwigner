function result = plot_consistency(events, cvdata, spikes, Fs, settings, varargin)
%PLOTCONSISTENCY is an auxiliary function to plot the firing and log Likelihood
%           and display the real class and predicted classification of events.
%
%           INPUTS:
%           P           : DataHigh type struct of dimension (1 x num events)
%           models      : cell containing trained GPFA models with dimension (1 x num models)
%
%           OPTIONS:
%           scaleK      : scale GPFA kernel (scale the speed of internal dynamics)
%           scaleRate   : scale the firing rate
%           useAllTrials: evaluate both training and test trials
%
%           OUTPUT:
%           stats       : a struct with dimensions (1 x folds) including the fields
%                       conf_matrix, class_output, real_label, and posterior, which are
%                       the classification confusion matrix where positive samples correspond
%                       to right alternations whereas negative samples are left alternations;
%                       output of the classifier {1:right, 2:left}, real label, and the
%                       log posterior P(data|model).%
%see also branch2, branch2_cleaned.m
%Stippinger Marcell, 2016

nlaps = length(events);

% Construct blurring window.
gaussFilter = gausswin(3);
gaussFilter = gaussFilter / sum(gaussFilter); % Normalize.

quantiles = [0.025 0.25 0.50 0.75 0.975];

%tot_times = (1:length(varargin{1}))/Fs;
tot_likelihood = [];

fig = figure('position',[100,100,1536,512]);

ax1 = subplot(1,6,1:5);
for ilap = 1 : nlaps
    lap_likelihood = cvdata(ilap).likelihood;
    lap_tolerance = cvdata(ilap).tolerance;
    lap_rate = sum(settings.bin_size * spikes(ilap).y .^ 2,1);

    % Do the blur.
    %smoothedVector = conv(lap_rate, gaussFilter, 'same');
    smoothedVector = lap_rate;

    T = length(lap_likelihood);
    lap_time_grid = events{ilap}(1,1)/Fs + (1:T)*settings.bin_size;

    h1 = plot(lap_time_grid,lap_likelihood,'k-','DisplayName','logLike');
    hold on
    h2 = plot(lap_time_grid,smoothedVector,'g-','DisplayName','Spike Count');
    plot(lap_time_grid,lap_likelihood-0.5*lap_tolerance,'k-');
    plot(lap_time_grid,lap_likelihood+0.5*lap_tolerance,'k-');
    tot_likelihood = [tot_likelihood; lap_likelihood];
    if ilap == 1
        hvector = [h1, h2];
    end
end

xlabel('Time (s)');
ylabel('Log likelihood  &  Spike frequency');

tot_class = zeros(size(tot_likelihood));
result = struct('all_mean',nanmean(tot_likelihood),...
                'all_stdev',nanstd(tot_likelihood),...
                'all_quan',quantile(tot_likelihood,quantiles));

    %result.all_mean = nanmean(tot_likelihood);
    %result.all_stdev = nanstd(tot_likelihood);
    
if ~isempty(varargin)
    state = varargin{1};
    tot_times = (1:length(state)).*settings.bin_size;
    sel = find(state);
    h1 = plot(tot_times(sel),-5*ones(size(sel)),'c+','DisplayName','Sharp wave');
    tot_class(sel)=1;
    result.spw_quan = quantile(tot_likelihood(sel),quantiles);
    result.spw_mean = nanmean(tot_likelihood(sel));
    result.spw_stdev = nanstd(tot_likelihood(sel));
    
    if length(varargin) > 2
        win = varargin{3};
    else
        win = 1;
    end
    state = conv(state,ones(win,1),'same');
    sel = find(~state);
    result.spk_quan = quantile(tot_likelihood(sel),quantiles);
    result.spk_mean = nanmean(tot_likelihood(sel));
    result.spk_stdev = nanstd(tot_likelihood(sel));
end
if length(varargin) > 1
    state = varargin{1} & varargin{2};
    sel = find(state);
    h2 = plot(tot_times(sel),-10*ones(size(sel)),'bo','DisplayName','Replay');
    tot_class(sel)=2;
    result.replay_quan = quantile(tot_likelihood(sel),quantiles);
    result.replay_mean = nanmean(tot_likelihood(sel));
    result.replay_stdev = nanstd(tot_likelihood(sel));

    state = varargin{1} & ~varargin{2};
    sel = find(state);
    h3 = plot(tot_times(sel),-15*ones(size(sel)),'rx','DisplayName','Inconsistent');
    tot_class(sel)=3;
    result.incons_quan = quantile(tot_likelihood(sel),quantiles);
    result.incons_mean = nanmean(tot_likelihood(sel));
    result.incons_stdev = nanstd(tot_likelihood(sel));

end

%hAnnotation = get(object_handle,'Annotation');
%hLegendEntry = get(hAnnotation','LegendInformation');
%set(hLegendEntry,'IconDisplayStyle','off')

spwclass = fitcdiscr(tot_likelihood,tot_class);
pred_class = predict(spwclass,tot_likelihood);
result.cm = confusionmat(tot_class,pred_class);

sym = {'c+', 'bo', 'rx'};
for i = 1:3
    sel = pred_class==i;
    sel = find(sel);
    if ~isempty(sel)
        plot(tot_times(sel),-20-i*5*ones(size(sel)),sym{i});
    end
end

hold off
xlim(ax1,[0,240]);
ylim(ax1,[-300,50]);

ax2 = subplot(1,6,6);

errorbar(1,result.all_mean,result.all_stdev,'ks','DisplayName','All Spiking');
hold on
try
errorbar(2,result.spk_mean,result.spk_stdev,'k.','DisplayName','Std.Spiking');
errorbar(3,result.replay_mean,result.replay_stdev,'bo','DisplayName','Replay');
errorbar(4,result.incons_mean,result.incons_stdev,'rx','DisplayName','Inconsistent');
catch
errorbar(2,result.spw_mean,result.spw_stdev,'c+','DisplayName','Sharp waves');    
end
hold off

linkaxes([ax1,ax2],'y');
% WONTFIX: We need to place legend into subplot 2 because subplot 1 spans
% several fields and some Matlab callback function is broken, zoom throws.
hl = legend(ax2, [hvector, h1, h2, h3]);

