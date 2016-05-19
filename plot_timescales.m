function plot_timescales(model, color, name, varargin)
%PLOTLATENT is an auxiliary function to plot the timescales fitted
%           by the GPFA procedure
%
%
%Marcell Stippinger, 2016


bin_size_ms     = 25; % 20ms at 1250 Hz in samples, refer to gpfa_mod.m
extra_opts      = assignopts(who,varargin);

N = length(model);
for v = 1 : N;
   %subplot(1,N,v)
   for l = 1 : length(model{v}.params)
       gamma = model{v}.params{l}.gamma;
       tau = bin_size_ms./sqrt(gamma);
       %plot(tau, v, 'o', 'color', color(l,:), 'displayname', name{v}), hold on   
       plot(tau, ones(1,length(tau))*v, 'o', 'color', color(l,:)), hold on   
   end
end

hold off
ax = gca;
ax.YTick = 1:N;
ax.YTickLabel = name;
%xlim([0 size(Xorth,2)])
xlabel('tau (ms)')
ylim([0 N+1.5])
ylabel('model')
legend({'fold 1', 'fold 2', 'fold 3'})

