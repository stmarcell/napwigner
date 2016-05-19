function R = reshape_laps(D, max_length)
%RESHAPE_LAPS is a utility to split long laps into smaller units
%             improving the speed of GPFA
%
% Author:
% Marcell Stippinger, 2016

%TODO: implement overlapping splittings to simulate continuity over time

nlaps           = length(D);
separators      = cell(1,nlaps);
npieces         = 0;

for ilap = 1 : nlaps
    lap_length  = D(ilap).T;
    lap_pieces  = ceil(lap_length*1./max_length);
    npieces     = npieces + lap_pieces;

    separators{ilap} = round(linspace(1,lap_length+1,lap_pieces+1));
end

% copy only the field names from D
R               = struct(D(1:0));
npieces         = 0;

for ilap = 1 : nlaps
    lap_sepa    = separators{ilap};
    lap_pieces  = length(lap_sepa)-1;

    %R(npieces+1:npieces+lap_pieces) = repmat(D(ilap),lap_pieces);
    
    for j = 1 : lap_pieces
        R(npieces+j) = D(ilap);
        R(npieces+j).T = lap_sepa(j+1)-lap_sepa(j);
        R(npieces+j).y = D(ilap).y(:,lap_sepa(j):lap_sepa(j+1)-1);
    end
    npieces = npieces + lap_pieces;
end

