function Xtats = gpfaConsistency(P, models, varargin)
%gpfaConsistency  Given a data struct P containing spike count vectors and
%GPFA models, this file computes the log likelihood using a sliding window
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
%                       log posterior P(data|model).
%
% See classGPFA too.
%Version 1.0 Stippinger Marcell, 2016
scaleK       = 1.0;
scaleRate    = 1.0;
useAllTrials = false;
width = 7;
assignopts(who,varargin);

folds        = length(models{1}.params);
scale        = (scaleK ~= 1.0) || (scaleRate ~= 1.0);
if scale
    fprintf('Scaling the GP Kernel with %2.2f, rates with %2.2f\n',scaleK, scaleRate);
end
if useAllTrials
   disp('Warning: Using all the trials for testing');
end

n_laps      = length(P);
v_laps      = [P.trialId];
%model_like  = cell(length(models), n_laps);
%model_tol   = cell(length(models), n_laps);
Xtats = struct([]);

pre = floor((width-1)/2.0);
post = ceil((width-1)/2.0);
disp([pre, post]);

for m = 1 : length(models)
    %likelikehood   = -Inf*ones(folds, n_laps);
    likelikehood   = cell(folds, n_laps);

    for ifold = 1 : folds
        
        if ~useAllTrials
            usedlaps    = models{m}.trainTrials{ifold};
            unseenP     = ones(1,n_laps);
            for u = 1 : length(usedlaps)
                %u_idx = find(v_laps == usedlaps(u));
                %unseenP(u_idx) = 0;
                unseenP(v_laps == usedlaps(u)) = 0;
            end
            unseenP = find(unseenP ==1);
        else
            unseenP = 1:n_laps;
        end
        
        for p = 1 : length(unseenP) 
        
            %select the model parameters from the fold#1 
            model = models{m}.params{ifold};
            lap   = unseenP(p);
            %rescale time scale of the GP if needed.
            if scale
               model.gamma = model.gamma .* (scaleK .^ 2);
               model.d = model.d .* sqrt(scaleRate);
               model.C = model.C .* sqrt(scaleRate);
            end
            
            lap_length = size(P(lap).y,2);
            
            Plap = P(lap);
            lap_likelihood = nan(lap_length,1);
            for c = 1 : lap_length
                sbin = max(c-pre,1);
                ebin = min(c+post,lap_length);
                Plap.T = ebin-sbin+1;
                Plap.y = P(lap).y(:,sbin:ebin);
            
                % sorry, exactInferenceWithLL does not return per lap
                % likelihoods, so we need to do it in a for loop
                [traj, ll] = exactInferenceWithLL(Plap, model,'getLL',1);
                lap_likelihood(c) = ll / Plap.T;
                
            end
            likelikehood{ifold,lap} = lap_likelihood;

        end
        %remove trials used during training
    end
    
    for ilap = 1: n_laps
        %model_like(m,ilap) = max(likelikehood);
        model_like = nanmean([likelikehood{:,ilap}],2);
        model_tol  = (nanmax([likelikehood{:,ilap}],[],2)-nanmin([likelikehood{:,ilap}],[],2))/(folds-1);
        
        %Xtats.real_label     = type;
        Xtats(m,ilap).likelihood     = model_like;
        Xtats(m,ilap).tolerance      = model_tol;
    end
    [Xtats(m,:).real_label]     = P.type;

end

