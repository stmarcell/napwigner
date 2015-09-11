function [rate, duration] = normFiringRate(XT, YT, X_Run_Lap, Y_Run_Lap, int_at_maze,...
                               gridsR, gridsL,  ctrNeuron, trial, verbose)
%NORMFIRINGRATE Produced a normalized unsmoothed firing rate by couting the
%               number of spikes on N segments (spatial bins) and
%               normalized by the time the animal spent in that bin.
segments            = size(gridsR,3);
[numLaps, N]        = size(X_Run_Lap);
rate                = zeros(N, segments, numLaps); 
duration            = zeros(numLaps, segments);

for ilap = 1 : numLaps
   %postion of the animal per lap 
   xt          = XT(int_at_maze(ilap,1):int_at_maze(ilap,2));
   yt          = YT(int_at_maze(ilap,1):int_at_maze(ilap,2));
   
   if strcmp(trial{ilap}, 'right') || strcmp(trial{ilap}, 'errorLeft')
       grid = gridsR;
   else
       grid = gridsL;
   end
   for iroi = 1 : segments
       ROI = grid(:,:,iroi);
       insideTrack     = inpolygon(xt,yt,ROI(1,:),ROI(2, :));        
       duration(ilap,iroi)  = sum(insideTrack)+1;
   end
   for icell = 1 : N     
       
       show = 0;
       if icell == ctrNeuron(1) && ilap==ctrNeuron(2);
           figure(ilap), hold on;
           show = 1;
       end
       
       % position of animal at each spike
       x            = X_Run_Lap{ilap,icell};
       y            = Y_Run_Lap{ilap,icell};
       
       %depeding on the arm of the maze, the respective bins are use to
       %count the spikes (gridsR, gridsL), centers and duration (c, d) are
       %only computed for one cell since all are the same during one lap.
       if strcmp(trial{ilap}, 'right') || strcmp(trial{ilap}, 'errorLeft')
           cnt= countROIspks(x, y, xt, yt, gridsR, show);
       else
           cnt= countROIspks(x, y, xt, yt, gridsL, show); 
       end
      %time normalization 
       if length(x) - sum(cnt) > 10
          fprintf('WARNING: Counting spikes cell %d, lap %d, true %d, counted %d\n', icell, ilap,length(x), sum(cnt)) 
       end
          
       rate(icell, :, ilap) = cnt./duration(ilap, :);
       if icell == ctrNeuron(1) && ilap==ctrNeuron(2);
            fprintf('Counting => Lap %d, cell %d, norm. f. rate %3.3f \n'...
                    ,ilap, icell, sum(rate(icell, :, ilap)))
       end
       if verbose; 
           fprintf('Lap %d, cell %d, spikes %d, counted %d\n', ilap, icell,length(x), sum(cnt));
       end

   end
   %save the mean position of the animal inside the bin    
end

end

function count= countROIspks(x, y, xt, yt, grid, show)
    %COUNTROISPKS Counts the spikes inside a grid partition
    %   
    %   COUNT = countROIspks(x, y, grid, show)
    %           x and y are the coordinates of the spike event
    %           grid with dimension 2 x 5 x S, containes the
    %           partitions of the area in S segments. Show enables the
    %           ploting of the countign and ROI regions along with the spks
    %           xt, yt are the spatial coordinates of the animal
    %           to calculate the time spent on each
    %           bin.
    %
    %Rube Pinzon
    %version 1.0 2015
    
    xcopy   = x;
    ycopy   = y;
    numROIs = size(grid,3);
    count   = zeros(1, numROIs);
    for iroi = 1 : numROIs
        ROI = grid(:,:,iroi);
        insideIndex = inpolygon(x,y,ROI(1,:),ROI(2, :));
        %to avoid the indefinite value, add +1 to all the bins        
        
        %remove counted spikes
        x(insideIndex) = [];
        y(insideIndex) = [];
        count(iroi)    = sum(insideIndex); 
                            
        if show
            text(mean(ROI(1,:)),mean(ROI(2,:)),num2str(count(iroi)), 'color', 'r')
            plot(ROI(1,:),ROI(2,:), 'r'), hold on
            plot(xt,yt)
            plot(x,y,'xb')
        end   
    end

end




