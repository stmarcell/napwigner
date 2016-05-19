%GPFA4SYNTHETICDIMENSIONALITY This script reads the fitted GPFA models and
%plots the log likelihood of the different latent dimensionality.
%
%        DESCRIPTION: This script can be used for dimesionality estimation.
%Version 1.0 Marcell Stippinger


clc, close all; %clear all;
run('gpfa4syntheticSectionSettings.m');

%========================Paramteres and variables==========================

basepath        = '~/marcell/_Data_ubi_hpc/';
workpath        = '~/marcell/napwigner/work/';
[files, roots, animals] = get_matFiles(basepath,'spike_.*\.dat');
name_save_file  = 'trainedGPFA';

maxdim = 30;

%%
% ========================================================================%
%============== (2)  Load / Save data per project  =======================%
%=========================================================================%

for ianimal = 1:length(animals)

    fprintf('\nSelecting %d: %s\n\n',ianimal,files{ianimal});

    project         = strrep(animals{ianimal},'.dat','');
    savepath        = [workpath project '/'];


    M = cell(maxdim,1);
    like_train = cell(maxdim,1);
    like_test = cell(maxdim,1);
    legend_ = cell(maxdim,1);

    for zDim = 1:maxdim
        fn          = [project '_' ...
                       name_save_file '_' sprintf('%02d',zDim) '.mat'];

        fprintf('Will load from %s\n', fn);
        try
            % Load saved model parameters and info
            tmp = load([savepath fn], 'M');
            M{zDim} = tmp.M.all;
            % Get the likelihood train as a matrix, with folds on one axis,
            % iteration steps on the other axis (transpose either cell contents
            % or cell indexing)
            % tr = cellfun(@transpose,M{zDim}.like_train,'un',0);
            like_train{zDim} = sum(cell2mat(M{zDim}.like_train'),1);
            like_test{zDim} = sum(M{zDim}.like_test,2);
            % alternative to sum: nanmean
        catch
            like_train{zDim} = nan(1,200); %sum(nan(max_iter,settings.n_folds),1);
            like_test{zDim} = nan;
            disp('file not found');
        end
        legend_{zDim} = sprintf('%02d', zDim);
    end

    like_train = cell2mat(like_train); %#ok<NASGU>
    like_test = cell2mat(like_test); %#ok<NASGU>
    fn              = [project '_' ...
                       name_save_file '_' 'dim' '.mat'];
    save([savepath fn], 'like_train', 'like_test', 'legend_');

end


%%
% ========================================================================%
%============== (3)     Load all data             ========================%
%=========================================================================%


name_save_file  = 'trainedGPFA';

like_train = {};
like_test = {};
legend_ = {};

for ianimal = 1:length(files)
    project         = strrep(animals{ianimal},'.dat','');
    savepath        = [workpath project '/'];

    fn              = [project '_' ...
                       name_save_file '_' 'dim' '.mat'];

    fprintf('Will load from %s\n', fn);
    try
        tmp = load([savepath fn]);
        like_train{end+1} = tmp.like_train(:,end); %#ok<SAGROW>
        like_test{end+1} = tmp.like_test(:,end); %#ok<SAGROW>
        
        legend_{end+1} = animals{ianimal}; %#ok<SAGROW>
    catch EM
        fprintf('file not found: %s', EM.message);
    end
end


%%
% ========================================================================%
%============== (4)   Plot dim estimates          ========================%
%=========================================================================%


fign = sprintf('ubi_hpc_all');
plot(cell2mat(like_test),'-');
legend(legend_{:},'Location','southeast');
xlabel('#latent dimensions for GPFA');
ylabel('cross-validated log likelihood');
savefig([fign '.fig']);
saveas(gcf,[fign '.png'],'png');
saveas(gcf,[fign '.pdf'],'pdf');

%======================== Plot per synthetic dim ==========================

files_per_dim = 5;
synthetic_dim = 5;
offset = 2;

for idim = 1:synthetic_dim
    fign = sprintf('ubi_hpc_%d',idim);
    start = files_per_dim*(idim-1)+offset;
    data = cell2mat(like_test(start:start+files_per_dim-1));
    idx = ~isnan(data(:,1));
    plot(find(idx),data(idx,:),'-');
    legend(legend_(start:start+files_per_dim-1),'Location','southeast');
    xlabel('#latent dimensions for GPFA');
    ylabel('cross-validated log likelihood');
    xlim([0,21]);
    savefig([fign '.fig']);
    saveas(gcf,[fign '.png'],'png');
    saveas(gcf,[fign '.pdf'],'pdf');
end

% dimensions in rows, files in columns
alike = cell2mat(like_test);

save([workpath 'dimensionality.dat'],'alike','-ascii');
