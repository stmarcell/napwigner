
function [leftT, rightT, failed_trial] = protoLap(XT, YT, length_run, trial, X_Run_Lap, Y_Run_Lap, int_at_maze, Fs, animal, isIntern)
% PROTOLAP Calcualtes the typical (mean) trial of the animal give all the
%          trials succeded.


numLaps = length(length_run);
N       = size(X_Run_Lap,2);
longest = max(length_run)*Fs;
xl      = zeros(1,longest) ;yl = zeros(1,longest); %vectors to calculate the mean trajectory
xr      = zeros(1,longest) ;yr = zeros(1,longest); 
r       = 0;l = 0;

failed_trial = [];
for ilap = 1 : numLaps
    duration    = linspace(0,length_run(ilap),longest); 
    %real animal position and interporaltion
    x       = XT(int_at_maze(ilap,1):int_at_maze(ilap,2));
    y       = YT(int_at_maze(ilap,1):int_at_maze(ilap,2));   
    xi      = spline(linspace(0,length_run(ilap),length(x)),x,duration);
    yi      = spline(linspace(0,length_run(ilap),length(y)),y,duration);        
    if strcmp(trial{ilap}, 'right') 
        xr = xr + xi; yr = yr + yi; r  = r + 1;
        %count spikes in the grids of the right arm        
    elseif strcmp(trial{ilap}, 'left') 
        xl  = xl + xi;yl = yl + yi; l  = l + 1;       
    else
        failed_trial =[failed_trial ilap]; %#ok<AGROW>
    end
    
    %plot(x, y,'color',[0.8, 0.8, 0.8]), hold on


end
win= 10;
leftT= [medfilt1(xl,win); medfilt1(yl,win)]'./l;
rightT= [medfilt1(xr,win); medfilt1(yr,win)]'./r;

% plot(leftT(:,1), leftT(:,2),'k','linewidth', 2), hold on
% plot(rightT(:,1), rightT(:,2),'k','linewidth', 2), 
% title(sprintf('Animal %s',animal))
% xlabel('x')

end


%MEDFILT1       One-dimensional median filter
%
%       y = MEDFILT(x)
%       y = MEDFILT(x, w)
%
%       median filter the signal with window of width W (default is 5).

% Copyright (C) 1995-2009, by Peter I. Corke
%
% This file is part of The Machine Vision Toolbox for Matlab (MVTB).
% 
% MVTB is free software: you can redistribute it and/or modify
% it under the terms of the GNU Lesser General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% MVTB is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU Lesser General Public License for more details.
% 
% You should have received a copy of the GNU Leser General Public License
% along with MVTB.  If not, see <http://www.gnu.org/licenses/>.
function m = medfilt1(s, w)
        if nargin == 1,
                w = 5;
        end
        
        s = s(:)';
        w2 = floor(w/2);
        w = 2*w2 + 1;

        n = length(s);
        m = zeros(w,n+w-1);
        s0 = s(1); sl = s(n);

        for i=0:(w-1),
                m(i+1,:) = [s0*ones(1,i) s sl*ones(1,w-i-1)];
        end
        m = median(m);
        m = m(w2+1:w2+n);

end