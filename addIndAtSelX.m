% addIndAtSelX.m
%
% Function that takes in specific x-coordinate value and adds
%  the corresponding index from the given indices vector
% Called by interactGetLegReversals(), meant to add selected point
%  from max/min indices list
% Converts x-coordinate value to index by using input x-coords to index
%  mapping vector, taking window around that index and finding max/min
%  within that window
% Adds max/min to indices list
%
% INPUTS:
%   inds - vector of indices to add/remove from
%   xCoords - vector of x coordinates, allows mapping between indices
%       (indices of this vector match indices of inds)
%   legPos - vector of leg positions
%   thisX - x coordinate of point to add/remove
%   maxMinStr - string 'max' or 'min' for whether point is max or min
%
% OUTPUTS:
%   thisInd - the index corresponding to thisX that was added
%   inds - updated inds vector, with thisInd added
%
% CREATED: 9/30/21 - HHY
%
% UPDATED:
%   9/30/21 - HHY
%
function [thisInd, inds] = addIndAtSelX(inds, xCoords, legPos, ...
    thisX, maxMinStr)

    % some parameters
    % size of window before and after identified index to look for max/min
    winLen = 5; 

    % get index closest to thisX
    [~, closeInd] = min(abs(xCoords - thisX));
    
    % get window to look for max/min
    winStartInd = closeInd - winLen;
    % set to start of trial if start of window ends up before start of
    %  trial
    if (winStartInd < 1)
        winStartInd = 1;
    end
    winEndInd = closeInd + winLen;
    % set to end of trial if end of window ends up after end of trial
    if (winEndInd > length(xCoords))
        winEndInd = length(xCoords);
    end
    
    % leg positions within window
    locLegPos = legPos(winStartInd:winEndInd);
    
    % find index of max/min within window
    switch maxMinStr
        case 'max'
            [~,locInd] = max(locLegPos);
        case 'min'
            [~,locInd] = min(locLegPos);
    end
    
    % convert local index to index as defined by xCoords and legPos
    thisInd = winStartInd + locInd - 1;
    
    % add this index into inds
    inds = [inds; thisInd];
    inds = sort(inds); % re-sort into ascending order
end