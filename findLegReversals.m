% findLegReversals.m
%
% Function to find when leg position in x reaches max and min, i.e.
%  reverses direction
% Algorithm: smooth leg position aggressively using sliding window; using
%  overlapping windows, find index of max, min; if index of max/min isn't 
%  the first or last index of the window, flag the window; find the max/min
%  within flagged windows in raw leg position data
% Then, refine calls using 1st derivative of leg positions: must exceed
%  specific thresholds of movement around max/min
%
% Called by interactGetLegReversals()
%
% INPUTS:
%   legXPos - column vector of single leg's X position
%   legXVel - column vector of same leg's X velocity
%   legRevParams - struct of parameters for getting these leg reversals
%       movAvgWinLen - length of window, in frames, for moving average
%       maxminWinLen - length of window, in frames, for finding max/min
%       adjThresh - threshold for adjacent indicies
%       maxPosVelThresh - threshold of what 1st derivative (vel) should 
%           exceed in the positive direction for a position maximum
%       maxNegVelThresh - in negative direction
%       minPosVelThresh - in pos dir, for leg pos min
%       minNegVelThresh - in neg dir, for leg pos min
%       numNegVelFrames - num of frames before and after max/min to check 
%           for negative velocity thresh
%       numPosVelFrames - num frames to check for positive velocity thresh
%
% OUTPUTS:
%   maxInds - column vector of indicies corresponding to X position maxima
%   minInds - column vector of indicies corresponding to X position minima
%
% CREATED: 9/29/21 - HHY
%
% UPDATED:
%   9/29/21 - HHY
%
function [maxInds, minInds] = findLegReversals(legXPos, legXVel, ...
    legRevParams)

end
