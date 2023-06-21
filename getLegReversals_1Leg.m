% getLegReversals_1Leg.m
% 
% Helper function for getLegReversals() and interactGetLegReversals(), to
% extract leg reversals for one leg. Replacement for findLegReversals().
% Uses MATLAB's findpeaks() instead of custom-written leg reversal finding
%
% INPUTS:
%   legXPos - column vector of single leg's X position
%   notMoveInd - indicies for when fly isn't moving, column vector
%   legRevParams - struct of parameters for getting these leg reversals
%       minProm - parameter MinPeakProminence of findpeaks
%       minDist - parameter MinPeakDistance of findpeaks
%
% OUTPUTS:
%   maxIndices - column vector of indicies of to X position maxima
%   minIndices - column vector of indicies cof to X position minima
%
% CREATED: 6/21/23 - HHY
%
% UPDATED:
%   6/21/23 - HHY
%
function [maxIndices, minIndices] = getLegReversals_1Leg(legXPos, ...
    notMoveInd, legRevParams)

    % get max ind for this leg
    [~, maxPeakInd] = findpeaks(legXPos, ...
        'MinPeakProminence', legRevParams.minProm, ...
        'MinPeakDistance', legRevParams.minDist);
    % remove all ind from when fly not moving
    % logical, true when maxInd is also not moving
    maxNotMovLog = ismember(maxPeakInd, notMoveInd);
    % remove all not moving maxInd
    maxIndices = maxPeakInd(~maxNotMovLog);

    % get min ind for this leg (find peaks in negative X pos)
    [~, minPeakInd] = findpeaks(-1*legXPos, ...
        'MinPeakProminence', legRevParams.minProm, ...
        'MinPeakDistance', legRevParams.minDist);
    % remove all ind from when fly not moving
    % logical, true when maxInd is also not moving
    minNotMovLog = ismember(minPeakInd, notMoveInd);
    % remove all not moving maxInd
    minIndices = minPeakInd(~minNotMovLog);
end