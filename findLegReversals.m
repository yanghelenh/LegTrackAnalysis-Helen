% findLegReversals.m
%
% Function to find when leg position in x reaches max and min, i.e.
%  reverses direction
% Algorithm: smooth leg position aggressively using sliding window; using
%  overlapping windows, find index of max, min; if index of max/min isn't 
%  the first or last index of the window, flag the window; find the max/min
%  within flagged windows in raw leg position data. Remove max/min in
%  not-moving bouts
% Then, refine calls using 1st derivative of leg positions: must exceed
%  specific thresholds of movement around max/min
%
% Called by interactGetLegReversals()
%
% INPUTS:
%   legXPos - column vector of single leg's X position
%   legXVel - column vector of same leg's X velocity
%   notMoveInd - indicies for when fly isn't moving, column vector
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
%   maxIndices - column vector of indicies of to X position maxima
%   minIndices - column vector of indicies cof to X position minima
%
% CREATED: 9/29/21 - HHY
%
% UPDATED:
%   9/30/21 - HHY
%
function [maxIndices, minIndices] = findLegReversals(legXPos, legXVel, ...
    notMoveInd, legRevParams)

    % length of trial, in frames
    numFrames = length(legXPos);

    % get moving average for leg X position
    legXMovAvg = movmean(legXPos, legRevParams.movAvgWinLen);
    
    
    % generate matrix from smoothed leg position: turn each position into
    %  window, will be maxminWinLen shorter than original
    
    % preallocate
    legXWins = zeros(length(legXMovAvg)-legRevParams.maxminWinLen, ...
        legRevParams.maxminWinLen);
    
    % populate legXWins
    for i = 1:size(legXWins, 1)
        legXWins(i,:) = legXMovAvg(i:i+legRevParams.maxminWinLen-1);
    end
    
    % for each window, find index of max and min
    [~, legXWinsMaxWinInds] = max(legXWins, [], 2);
    [~, legXWinsMinWinInds] = min(legXWins, [], 2);
    
    % find indicies where max/min is not 1 or window length (i.e. has local
    %  max/min within it
    legXWinsMaxInds = find(legXWinsMaxWinInds ~= 1 & ...
        legXWinsMaxWinInds ~= legRevParams.maxminWinLen);
    legXWinsMinInds = find(legXWinsMinWinInds ~= 1 & ...
        legXWinsMinWinInds ~= legRevParams.maxminWinLen);
    
    % merge adjacent windows (within legRevParams.adjThresh of each other)
    %  return start and end indicies of windows
    [maxWinsStartInd, ~, maxWinsEndInd, ~] = mergeAdjVecVals(...
        legXWinsMaxInds, legRevParams.adjThresh);
    [minWinsStartInd, ~, minWinsEndInd, ~] = mergeAdjVecVals(...
        legXWinsMinInds, legRevParams.adjThresh);
    
    
    % for each window, find the index in the position vector where the
    %  max/min actually occured
    
    % preallocate
    maxInds = zeros(length(maxWinsEndInd), 1);
    minInds = zeros(length(minWinsEndInd), 1);
    
    % find max within each window
    for i = 1:length(maxWinsEndInd)
        % get start index of window
        startInd = maxWinsStartInd(i);
        
        % get end index of window: windows defined by start index, account
        % for window size
        endInd = maxWinsEndInd(i) + legRevParams.maxminWinLen - 1;
        
        % get position values in this window
        winXPos = legXPos(startInd:endInd);
        
        % get max within window
        [~, thisMaxInd] = max(winXPos);
        
        % save max, convert to index within leg position vector
        maxInds(i) = thisMaxInd + startInd - 1;
    end
    
    % find min within each window
    for i = 1:length(minWinsEndInd)
        % get start index of window
        startInd = minWinsStartInd(i);
        
        % get end index of window: windows defined by start index, account
        % for window size
        endInd = minWinsEndInd(i) + legRevParams.maxminWinLen - 1;
        
        % get position values in this window
        winXPos = legXPos(startInd:endInd);
        
        % get min within window
        [~, thisMinInd] = min(winXPos);
        
        % save min, convert to index within leg position vector
        minInds(i) = thisMinInd + startInd - 1;
    end
    
    % remove indicies within not-moving bouts
    % logical, true when maxInd is also not moving
    maxNotMovLog = ismember(maxInds, notMoveInd);
    % remove all not moving maxInd
    maxIndsMov = maxInds(~maxNotMovLog);
    
    minNotMovLog = ismember(minInds, notMoveInd);
    minIndsMov = minInds(~minNotMovLog);
    

    
    % refine calls with 1st derivative thresholds
    
    % intialize: indicies into maxIndsMov of maxima to remove
    maxRmvInd = [];
    
    % loop through all maxima
    for i = 1:length(maxIndsMov)
        % window before maximum
        bfWinStart = maxIndsMov(i) - legRevParams.numPosVelFrames;
        bfWinEnd = maxIndsMov(i) - 1; % index right before
        % window after maximum
        afWinStart = maxIndsMov(i) + 1; % index right before
        afWinEnd = maxIndsMov(i) + legRevParams.numNegVelFrames;
        
        % if maximum occurs at first or last frame of the trial, only use 1
        %  window to check if maximum is valid
        if (maxIndsMov(i) == 1) % max is first frame of trial
            % only check that after window is valid
            afAvgVel = mean(legXVel(afWinStart:afWinEnd));
            % if invalid, save index
            if (afAvgVel > legRevParams.maxNegVelThresh)
                maxRmvInd = [maxRmvInd; i];
            end
            % continue to next iteration of loop
            continue;
        elseif (maxIndsMov(i) == numFrames) % max is last frame of trial
            % only check that before window is valid
            bfAvgVel = mean(legXVel(bfWinStart:bfWinEnd));
            % if invalid, save index
            if (bfAvgVel < legRevParams.maxPosVelThresh)
                maxRmvInd = [maxRmvInd; i];
            end
            % continue to next iteration of loop
            continue;
        end
    
        % check that windows are valid (contained b/w 1 and number of frames)
        %  if max is w/in numVelFrames of start or end of trial, shorten window
        %  accordingly (but can't be length zero, b/c taken care of with check
        %  on whether max is at start or end of trial

        % set before window to start with trial start
        if (bfWinStart < 1) 
            bfWinStart = 1;
        end
        % set after window to end with trial end
        if (afWinEnd > numFrames)
            afWinEnd = numFrames;
        end

        % get average velocities within window
        bfAvgVel = mean(legXVel(bfWinStart:bfWinEnd));
        afAvgVel = mean(legXVel(afWinStart:afWinEnd));

        % check that 1st derivative meets threshold for both windows, if not,
        %  add index to those to remove
        if (afAvgVel > legRevParams.maxNegVelThresh) || ...
                (bfAvgVel < legRevParams.maxPosVelThresh)
            maxRmvInd = [maxRmvInd; i];
        end
    end

    % remove all the indicies flagged for maxima
    maxIndices = maxIndsMov;
    maxIndices(maxRmvInd) = [];



    % for minima

    % initialize array of indicies to remove
    minRmvInd = [];

    % loop through all maxima
    for i = 1:length(minIndsMov)
        % window before minimum
        bfWinStart = minIndsMov(i) - legRevParams.numNegVelFrames;
        bfWinEnd = minIndsMov(i) - 1; % index right before 
        % window after minimum
        afWinStart = minIndsMov(i) + 1; % index right after
        afWinEnd = minIndsMov(i) + legRevParams.numPosVelFrames;

        % if minimum occurs at first or last frame of the trial, only use 1
        %  window to check if minimum is valid
        if (minIndsMov(i) == 1) % min is first frame of trial
            % only check that after window is valid
            afAvgVel = mean(legXVel(afWinStart:afWinEnd));
            % if invalid, save index
            if (afAvgVel < legRevParams.minPosVelThresh)
                minRmvInd = [minRmvInd; i];
            end
            % continue to next iteration of loop
            continue;
        elseif (minIndsMov(i) == numFrames) % min is last frame of trial
            % only check that before window is valid
            bfAvgVel = mean(legXVel(bfWinStart:bfWinEnd));
            % if invalid, save index
            if (bfAvgVel > legRevParams.minNegVelThresh)
                minRmvInd = [minRmvInd; i];
            end
            % continue to next iteration of loop
            continue;
        end

        % check that windows are valid (contained b/w 1 and number of frames)
        %  if max is w/in numVelFrames of start or end of trial, shorten window
        %  accordingly (but can't be length zero, b/c taken care of with check
        %  on whether max is at start or end of trial

        % set before window to start with trial start
        if (bfWinStart < 1) 
            bfWinStart = 1;
        end
        % set after window to end with trial end
        if (afWinEnd > numFrames)
            afWinEnd = numFrames;
        end

        % get average velocities within window
        bfAvgVel = mean(legXVel(bfWinStart:bfWinEnd));
        afAvgVel = mean(legXVel(afWinStart:afWinEnd));

        % check that 1st derivative meets threshold for both windows, if not,
        %  add index to those to remove
        if (afAvgVel < legRevParams.minPosVelThresh) || ...
                (bfAvgVel > legRevParams.minNegVelThresh)
            minRmvInd = [minRmvInd; i];
        end
    end

    % remove all the indicies flagged for minima
    minIndices = minIndsMov;
    minIndices(minRmvInd) = [];

end
