% interactGetLegReversals.m
%
% Function for getting max and min in leg positions, when leg reverses
%  direction
% Brings up GUI, 1 leg at a time, plotting leg X position, x's for computed
%  min and max. Starts with sliders to adjust parameters. When user clicks
%  button, changes to mode to select points to delete or add.
%
% INPUTS:
%   legTrack - struct of leg tracking data, output of preprocessLegTrack
%   moveNotMove - struct of moving/not-moving bout indicies, starts/ends
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
%   legIDs - struct of parameters, IDs legs
%       ind - indicies of legs, into raw position matricies, carried
%           throughout
%       names - names of each leg, matching with ind
%
% OUTPUTS:
%   maxInds - indicies (frames) of leg position maxima, all legs in 1
%       vector
%   minInds - indicies of leg position minima
%   maxWhichLeg - number indicating which leg max indicies belong to, same
%       size & matched with maxInds
%   minWhichLeg - which legs for min indicies
%
% CREATED: 9/29/21 - HHY
%
% UPDATED:
%   9/29/21 - HHY
%
function [maxInds, minInds, maxWhichLeg, minWhichLeg] = ...
    interactGetLegReversals(legTrack, moveNotMove, legRevParams, legIDs)

end
