function gpfaCompare(fileToRead1, fileToRead2)
%gpfaCompare(fileToRead1, fileToRead2) compares the parameters of two
%stored models
%
%  Imports data from the specified files and plots histogram of the ratios
%  of d, C, R and tau. Note, that while neurons are increasingly orddered
%  and therefore identic sets are matched, latent dimensions in d and C
%  might be unmatched in the two models. The same, no effort is made to
%  match CV folds to contain the same laps. Irrespective of this these
%  plots can give a hint on the rescaling.
%
%  Usage example:
%    gpfaCompare('~/marcell/napwigner/work/spike_SPW_D2_L4/spike_SPW_D2_L4_trainedGPFA_05.mat',...
%                '~/marcell/napwigner/work/spike_RUN_D2_L4/spike_RUN_D2_L4_trainedGPFA_05.mat')
%
%  Marcell Stippinger, 27-Apr-2016

% Import the file
Model1 = load('-mat', fileToRead1, 'M');
Model2 = load('-mat', fileToRead2, 'M');
Params1 = Model1.M.all.params;
Params2 = Model2.M.all.params;

nfolds = size(Params1, 2);
if nfolds ~= size(Params2, 2)
    print 'Models have different number of folds'
end

% Create new variables in the base workspace from those fields.
ratio_d = cell(1,nfolds);
ratio_C = cell(1,nfolds);
ratio_R = cell(1,nfolds);
ratio_R_diag = cell(1,nfolds);
ratio_tau = cell(1,nfolds);
for i = 1:nfolds
    ratio_d{i} = Params1{i}.d ./ Params2{i}.d;
    ratio_C{i} = Params1{i}.C ./ Params2{i}.C;
    ratio_R{i} = Params1{i}.R ./ Params2{i}.R;
    ratio_R_diag{i} = diag(ratio_R{i});
    ratio_tau{i} = sqrt(Params1{i}.gamma ./ Params2{i}.gamma);
end

figure();
subplot(2,2,1); hold on;
histogram([ratio_d{:}],linspace(0,8,81));
plot(mean(cell2mat(ratio_d)),1,'rx',median(cell2mat(ratio_d)),2,'bo');
xlabel('ratio M_1/M_2'); ylabel('count'); title({'d: expected value','(per CV fold)'});
hold off;

subplot(2,2,2); hold on;
histogram([ratio_C{:}],linspace(0,8,81));
plot(mean(cell2mat(ratio_C)),1,'rx',median(cell2mat(ratio_C)),2,'bo');
xlabel('ratio M_1/M_2'); ylabel('count'); title({'C: transformation','(per CV fold per latent dim)'});
hold off;
% note the ordering of latent dimensions might differ

subplot(2,2,3); hold on;
histogram([ratio_R_diag{:}],linspace(0,8,81));
plot(mean(cell2mat(ratio_R_diag)),1,'rx',median(cell2mat(ratio_R_diag)),2,'bo');
xlabel('ratio M_1/M_2'); ylabel('count'); title({'R: variance','(per CV fold)'});
hold off;

subplot(2,2,4); hold on;
histogram([ratio_tau{:}],linspace(0,32,81));
plot(mean(cell2mat(ratio_tau)),1,'rx',median(cell2mat(ratio_tau)),2,'bo');
xlabel('ratio M_1/M_2'); ylabel('count'); title({'1/tau: inverse timescale','(per CV fold per latent dim)'});
legend('count','mean','median');
hold off;
