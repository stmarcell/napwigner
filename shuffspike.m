function P = shuffspike(P)
%SHUFFSPIKE auxiliary function to independently and randomly permute the
%           time samples/bins in the variable data of the input structure P.
%
% Marcell Stippinger 2016


for lap = 1 : length(P)
    n_bins    = P(lap).T;
    n_clu     = size(P(lap).y,1);
    
    for clu = 1 : n_clu
        idx       = randperm(n_bins); 

        P(lap).y(clu,:) = P(lap).y(clu,idx); 

        %P(lap).shuffled_bins = idx;
    end
end
