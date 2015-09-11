
function proto_out = protoLap(XT, YT, length_run, trial, int_at_maze, Fs, show, roiDims)
% PROTOLAP Calcualtes the typical (mean) trial of the animal give all the
%          trials succeded.


numLaps = length(length_run);
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
        color = [0.8 0.8 1];
        %count spikes in the grids of the right arm        
    elseif strcmp(trial{ilap}, 'left') 
        xl  = xl + xi;yl = yl + yi; l  = l + 1; 
        color = [1 0.8 0.8];
    else
        failed_trial =[failed_trial ilap]; %#ok<AGROW>
    end
    if show
        plot(x, y,'color',color), hold on
    end

end
win= 10;
leftT= [medfilt1(xl,win); medfilt1(yl,win)]'./l;
rightT= [medfilt1(xr,win); medfilt1(yr,win)]'./r;

[gridsL, centersL, bin_mmL]      = get_grids(leftT,  1, show, roiDims, [1 0.6 0.6]);
[gridsR, centersR, bin_mmR]      = get_grids(rightT, 1, show, roiDims, [0.6 0.6 1]);

if show
    plot(leftT(:,1), leftT(:,2),'r','linewidth', 2), hold on
    plot(rightT(:,1), rightT(:,2),'b','linewidth', 2), 
    xlabel('x')
end
%[leftT, rightT, failed_trial, gridsL, gridsR, bin_mm]
proto_out.left_traj = leftT;
proto_out.right_traj = rightT;
proto_out.failures = failed_trial;
proto_out.rois_le = gridsL;
proto_out.rois_ri = gridsR;
proto_out.centers_le =  centersL;
proto_out.centers_ri =  centersR;
proto_out.bin_size_le =  bin_mmL;
proto_out.bin_size_ri =  bin_mmR;

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

function [grids, centers, bin_mm] = get_grids(X, connect, show, roiDims, color)
%   GETGRIDS creates rectagular segments along a instructive path
%
%       GRIDS = get_grids(X, segments, connect, show, roiDims)
%               X is the coordinates of the instructive path along which
%               the rectagles will be palced
%               segments in the number of rectangles to fit in th elenght
%               of the accumulated distance spanned by X
%               connect (1 or 0) is an option to connect adyacent
%               rectangles or leave them disconnected.
%               show (1 or 0) to show the rectagles and the instructive
%               path. roiDims = [width lenght] of rectangles
%
%Ruben Pinzon
%version 1.0@2015
    segments = roiDims(1);
    deltaX   = diff(X(:,1));
    deltaY   = diff(X(:,2));
    accdist  = cumsum(sqrt(deltaX.^2 + deltaY.^2));
    accdist  = accdist - accdist(1);
    bin_mm   = accdist(end)/(segments);
    fprintf('Distance traveled by the animal %3.3f mm\n', accdist(end))   
    roiDims(1) = bin_mm;

    grids    = zeros(2, 5, segments);
    centers  = zeros(2, segments);
    border_old = 1;
    for ibin = 1 : segments
        
        border_new = find(diff(accdist <= ibin*bin_mm));
        if isempty(border_new)
            border_new = length(X);
        end
%         plot(X([border_old, border_new],1),X([border_old, border_new],2),'k'), hold on
        
        deltaXX = X(border_new,1) - X(border_old,1);
        deltaYY = X(border_new,2) - X(border_old,2);
        offsetX  = sum(X([border_old border_new],1))/2;
        offsetY  = sum(X([border_old border_new],2))/2;
        angle    = atan2(deltaYY,deltaXX);
        
        rotation = [cos(angle) -sin(angle); 
                  sin(angle)  cos(angle)];
        translation = repmat([offsetX; offsetY],1,5); 
        
        ROI    = rotation*getROI(0, 0, roiDims)' + translation;
        
%         ROI    = rotation*ROI' + translation;

        if ibin > 1 && connect
            ROI(:,[4, 3]) = oldROI;
        end
        grids(:, :, ibin) = ROI;
        
        centers(:,ibin) = polygonCentroid(ROI(1,:),ROI(2, :));
        oldROI          = ROI(:,1:2);
        border_old      = border_new;
        if show 
            plot(ROI(1,:),ROI(2,:), 'color', color)
            text(centers(1,ibin), centers(2,ibin), num2str(ibin),'color','k','fontsize',9)
        end
        %count spikes in the grids of the central arm
    end
end

function ROI = getROI(xo, yo, roiDims)
%GETROI creates a rectagle centered at xo yo and roiDims = [width lenght]
ROI = [     
       xo + roiDims(1)/2,  yo + roiDims(2)/2;
       xo + roiDims(1)/2,  yo - roiDims(2)/2;
       xo - roiDims(1)/2,  yo - roiDims(2)/2;
       xo - roiDims(1)/2,  yo + roiDims(2)/2;
       xo + roiDims(1)/2,  yo + roiDims(2)/2
       ];
end

function [centroid, area] = polygonCentroid(varargin)
    %POLYGONCENTROID Compute the centroid (center of mass) of a polygon
    %
    %   CENTROID = polygonCentroid(POLY)
    %   CENTROID = polygonCentroid(PTX, PTY)
    %   Computes center of mass of a polygon defined by POLY. POLY is a N-by-2
    %   array of double containing coordinates of vertices.
    %
    %   [CENTROID AREA] = polygonCentroid(POLY)
    %   Also returns the (signed) area of the polygon. 
    %
    %   Example
    %     % Draws the centroid of a paper hen
    %     x = [0 10 20  0 -10 -20 -10 -10  0];
    %     y = [0  0 10 10  20  10  10  0 -10];
    %     poly = [x' y'];
    %     centro = polygonCentroid(poly);
    %     drawPolygon(poly);
    %     hold on; axis equal;
    %     drawPoint(centro, 'bo');
    % 
    %   References
    %   algo adapted from P. Bourke web page
    %
    %   See also:
    %   polygons2d, polygonArea, drawPolygon
    %

    %   ---------
    %   author : David Legland 
    %   INRA - TPV URPOI - BIA IMASTE
    %   created the 05/05/2004.
    %

    % Algorithme P. Bourke, vectorized version

    % HISTORY
    % 2012.02.24 vectorize code


    % parse input arguments
    if nargin == 1
        var = varargin{1};
        px = var(:,1);
        py = var(:,2);
    elseif nargin == 2
        px = varargin{1};
        py = varargin{2};
    end

    % vertex indices
    N = length(px);
    iNext = [2:N 1];

    % compute cross products
    common = px .* py(iNext) - px(iNext) .* py;
    sx = sum((px + px(iNext)) .* common);
    sy = sum((py + py(iNext)) .* common);

    % area and centroid
    area = sum(common) / 2;
    centroid = [sx sy] / 6 / area;
end