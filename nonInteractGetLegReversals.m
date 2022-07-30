% nonInteractGetLegReversals.m
%
% Function for getting max and min in leg positions, when leg reverses
%  direction
% Unlike interactGetLegReversals(), is not interactive. Returns max and min
%  given input parameters
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
%   maxIndsAll - indicies (frames) of leg position maxima, all legs in 1
%       vector
%   minIndsAll - indicies of leg position minima
%   maxWhichLeg - number indicating which leg max indicies belong to, same
%       size & matched with maxIndsAll
%   minWhichLeg - which legs for min indicies
%
% CREATED: 7/27/22 - HHY
%
% UPDATED:
%   7/27/22 - HHY
% 
function [maxIndsAll, minIndsAll, maxWhichLeg, minWhichLeg] = ...
    nonInteractGetLegReversals(legTrack, moveNotMove, legRevParams, legIDs)

    % initialize maxIndsAll, minIndsAll, maxWhichLeg, minWhichLeg
    maxIndsAll = [];
    minIndsAll = [];
    maxWhichLeg = [];
    minWhichLeg = [];

    % loop through all legs
    for i = 1:length(legIDs.ind)
        thisLegInd = legIDs.ind(i);

        % get max and min for this leg
        [maxInds, minInds] = findLegReversals(...
            legTrack.srnLegX(:,thisLegInd), ...
            legTrack.legXVel(:,thisLegInd), moveNotMove.legNotMoveInd, ...
            legRevParams);

        % return values for this, concatenate into vector for all legs
        maxIndsAll = [maxIndsAll; maxInds];
        minIndsAll = [minIndsAll; minInds];
        
        % which leg marker, generate vectors
        thisMaxWhichLeg = ones(size(maxInds)) * legIDs.ind(i);
        thisMinWhichLeg = ones(size(minInds)) * legIDs.ind(i);
        % concatenate into vector for all legs
        maxWhichLeg = [maxWhichLeg; thisMaxWhichLeg];
        minWhichLeg = [minWhichLeg; thisMinWhichLeg];
    end
end