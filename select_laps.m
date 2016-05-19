function laps = select_laps(BehavType, TrialType, goal, side, side_alt)
%SELECTLAPS select laps IDs for training based on criteria
%
%        BehavType: the behavior of the animal, see Typebehav_tx
%        TrialType: the behavior of the animal, see Typetrial_tx
%        SideType: the behavior of the animal, see Typeside_tx
%        goal: what is to be leaarned, options are: run, wheel, sharp.
%        side: 0 for both, 1 for left, 2 for right
%        side_alt: side of next lap
%
% Marcell Stippinger, 2016

side_select     = [1 2 2 1];
SideType        = side_select(TrialType);
if nargin < 4
    side        = 0;
end
if nargin < 5
    side_alt    = 0;
end

goal_mask = false(1, length(SideType));
side_mask = false(1, length(SideType));
if strcmp('run', goal)
    goal_mask = (BehavType == 2);
    if side ~= 0
        side_mask = (SideType == side);
    else
        side_mask = SideType > 0;
    end
elseif strcmp('wheel', goal)
    BehavFuture = [BehavType(2:end) 0];
    goal_mask = (BehavType == 2) & (BehavFuture == 2);
    SideFuture = [SideType(2:end) 0];
    if side ~= 0
        if side_alt == 0
        	  side_alt = 3 - side;
        end
        side_mask = (SideType == side) & (SideFuture == side_alt);
    else
        side_mask = (SideType > 0) & (SideFuture > 0);
    end
else
    error('Unknown goal');
end

laps = find(goal_mask & side_mask);
