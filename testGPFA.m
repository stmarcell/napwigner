% create toy example
[grd, Y, X, S] = GPFA.toyExample();
[grd, X] = grd.normFactors(X);

% fit model
model = GPFA('Tolerance', 1e-6);
model = model.fit(Y, grd.p, 'hist');
% model = model.fit(Y, grd.p, S);
[model, Xest] = model.normFactors(Y);

figure()
[pc,score,latent,tsquare] = princomp(YY(:,:,1));
plot(pc(:,1),'x')
%% diagnostic plots
cc = model.C' * grd.C;
[~, ndx] = max(abs(cc));
N = size(Y, 3);
for i = 1 : model.p
    subplot(3, model.p, i)
    cla, hold on
    plot(X(i, :), 'k')
    plot(sign(cc(ndx(i), i)) * Xest(ndx(i), :), 'r')
    axis tight
    plot(repmat((1 : N - 1) * model.T, 2, 1), ylim, 'k')
    ylabel(sprintf('Factor %d', i))
    xlabel('Time')
    xlim([0 N * model.T])
end

for i = 1 : model.p
    subplot(3, model.p, model.p + i)
    cla, hold on
    plot(grd.C(:, i), 'k')
    plot(sign(cc(ndx(i), i)) * model.C(:, ndx(i)), 'r')
    ylabel(sprintf('Loading %d', i))
    xlabel('Neuron #')
    axis tight
end

legend({'Ground truth', 'Model fit'})

subplot(3, model.p, 2 * model.p + 1)
cla, hold on
plot(grd.D(:), model.D(:), '.k')
axis tight
xlabel('Stim term: ground truth')
ylabel('Stim term: model fit')