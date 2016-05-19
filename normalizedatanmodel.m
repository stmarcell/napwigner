function [trialsout, varargout] = normalizedatanmodel(trialsin, varargin)
%ZEROMEANMODEL auxiliary function to set the mean of the observed variables
%           to zero and their std.dev. to 1 in the input structure.
%
%  Additionally GPFA models can be provided, which will be adjusted to
%  prefer zero mean and unit variance inputs. The covariance sructure of
%  the input can be analyzed then.
%
% Marcell Stippinger 2016


% Data
ntrials = length(trialsin);
trialsout = trialsin;
for m = 1 : ntrials
    trialy = trialsin(m).y;
    nbins = size(trialy,2);
    %nclusters = size(trialy,1);
    meanfr = mean(trialy,2);
    sdevfr = sqrt(var(trialy,0,2));
    sdevfr(sdevfr==0) = 1;
    trialsout(m).y = (trialy - repmat(meanfr,1,nbins)) ./ repmat(sdevfr,1,nbins);
end

% Models (any number)
nmodels = length(varargin);
for m = 1 : nmodels
    model = varargin{m};
    nfolds = length(model.params);
    for f = 1: nfolds
        model.params{f}.d(:,:)=0;
        nobservables = length(model.params{f}.R);
        model.params{f}.R=ones(nobservables);
    end
    varargout{m} = model;
end

