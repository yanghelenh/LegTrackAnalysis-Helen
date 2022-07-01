% convertNotMoveLogToBouts.m
% 
% Function that takes in logical indicatign when fly is not moving and
%  converts it a list of not moving bout start indices, end indices, and
%  all not moving indices
%
% Helper function for: findLegFtCmbNotMove() and findFlyNotMovingFt()
%
% INPUTS:
%   logVec - logical vector, true when fly not moving
%
% OUTPUTS:
%   startInd - start indices of not moving bouts
%   endInd - end indices of not moving bouts
%   allInd - all indices of not moving bouts
%
% CREATED: 7/1/22 - HHY
%
% UPDATED:
%   7/1/22 - HHY
%
function [startInd, endInd, allInd] = convertNotMoveLogToBouts(logVec)
    % get start and end indices of not move bouts
    startInd = find(diff(logVec) > 0.9) + 1;
    endInd = find(diff(logVec) < -0.9);

    % add start of trial if fly is not moving at start
    if (logVec(1))
        startInd = [1 startInd];
    end

    % add end of trial if fly is not moving at end
    if (logVec(end))
        endInd = [endInd length(legNotMoveLog)];
    end

    % convert notMoveLogical to indices
    allInd = find(logVec);
end