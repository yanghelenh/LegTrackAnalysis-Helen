% minMaxPosToSteps.m
%
% Function that takes indicies of min and max leg positions and converts
%  them into steps.
% Define step start as leg in back-most position (maxInds), as from test
%  data, that is where leg movements start when fly starts moving
% Filtered to ensure steps are w/in same moving bout
%
% INPUTS:
%   maxInds - indices for when leg positions at maxima, n x 1 vector
%       including all legs
%   minInds - indices for when leg positions at minima, n x 1 vector
%       including all legs
%   maxWhichLeg - n x 1 vector indicating which leg each max in maxInds
%       belongs to
%   minWhichLeg - n x 1 vector indicating which leg each min in minInds
%       belongs to
%   notMoveBout - m x 2 matrix of not moving bout starts (col 1) and ends
%       (col 2)
%   moveBout - q x 2 maxtrix of moving bout starts (col 1) and ends (col 2)
%
% OUTPUTS:
%   stepInds - p x 3 matrix of steps, with each step as row [start mid end]
%       defined by indices
%   stepWhichLeg - p x 1 vector indicating which leg each step belongs to
%
% CREATED: 10/1/21 - HHY
%
% UPDATED:
%   10/1/21 - HHY
%
function [stepInds, stepWhichLeg] = minMaxPosToSteps(maxInds, minInds, ...
    maxWhichLeg, minWhichLeg, notMoveBout, moveBout)

    % maximum possible number of steps, number of max or min
    if (length(maxInds) >= length(minInds))
        maxStepNum = length(maxInds);
    else
        maxStepNum = length(minInds);
    end
    
    % preallocate step arrays
    stepInds = zeros(maxStepNum, 3);
    stepWhichLeg = zeros(maxStepNum, 1);
    
    % initialize step counter
    counter = 1;
    
    % loop through all maxInds elements
    for i = 1:length(maxInds)
        
        % frame index of this potential start point (maximum)
        thisPotStart = maxInds(i);
        
        % get leg for this potential start point
        thisLeg = maxWhichLeg(i);
        
        % consider max and min inds only for this leg
        thisLegMaxInds = maxInds(maxWhichLeg == thisLeg);
        thisLegMinInds = minInds(minWhichLeg == thisLeg);
        
        % index of this potential start point into thisLegMaxInds
        thisPotStartInd = find(thisLegMaxInds == thisPotStart,1,'first');
        
        
        % check if this maxInd is within a not moving bout (shouldn't be,
        %  but could have been added manually; check)
        
        % index into notMoveBout for start point immediately preceding this
        %  maxInd
        thisNotMoveStartInd = find(thisPotStart > notMoveBout(:,1), 1, ...
            'last');
        % if there could be a not moving bout that this belongs to (empty
        %  if this maxInd falls earlier than start of first not moving
        %  bout; i.e. if fly is walking at start of trial)
        if ~(isempty(thisNotMoveStartInd))
            % get frame index value for start of not moving bout
            thisNotMoveStart = notMoveBout(thisNotMoveStartInd, 1);
            % get frame index value for end of not moving bout (same index,
            %  starts and ends paired)
            thisNotMoveEnd = notMoveBout(thisNotMoveStartInd, 2);
            
            % check if this potential start point falls within this not
            %  moving bout (is value less than thisNotMoveEnd?)
            % if yes, continue loop without executing rest
            if (thisPotStart < thisNotMoveEnd)
                continue;
            end
        end
        
        
        % get start and end indices of moving bout that this start point
        %  belongs to
        
        % index into moveBout for start point immediately preceding this
        %  maxInd
        thisMoveStartInd = find(thisPotStart >= moveBout(:,1), 1, 'last');
        % frame index of move bout start
        thisMoveStart = moveBout(thisMoveStartInd, 1);
        % frame index of move bout end; paired
        thisMoveEnd = moveBout(thisMoveStartInd, 2);
        
        
        % using this potential start point, find indices that define step
        %  mid point (from minInds) and end point (from maxInds)
        stepStartPt = thisPotStart;
        
        % check whether this step start point (as a max), is followed by
        %  another max or by a min
        % should be min, but if not, don't use this start point
        
        % check if there is next max point
        % if this start index is the last max point for this leg
        if (thisPotStartInd == length(thisLegMaxInds))
            % there must not be end point for this step, so don't count it
            continue;
        else % there are additional points
            nextMaxPt = thisLegMaxInds(thisPotStartInd + 1);
            nextMinPt = thisLegMinInds(find(thisLegMinInds > stepStartPt,...
                1, 'first'));
            % if max follows the start point instead of min, don't use this
            %  start point
            if (nextMaxPt < nextMinPt)
                continue;
            end
        end
        
        
        % find mid point - minInds value that immediate follows step start
        %  point
        % index into thisLegMinInds of this point
        midPtInd = find(thisLegMinInds > stepStartPt, 1, 'first');
        % if there is no such point, this step is cut off by the end of the
        %  trial
        if isempty(midPtInd)
            continue;
        end
        
        % convert this index into minInds into frame index of step mid
        %  point
        stepMidPt = thisLegMinInds(midPtInd);
        
        % check that this midpoint is within the same moving bout as the
        %  step start point
        if (stepMidPt < thisMoveStart) || (stepMidPt > thisMoveEnd)
            continue;
        end
        
        % check whether this step mid point (as min) is followed by another
        %  min or by a max
        % should be a max, but if not, don't use this start point
        % check whether there is a subsequent min point (if not, don't use)
        if (midPtInd ~= length(thisLegMinInds))
            nextMinPt = thisLegMinInds(midPtInd + 1);
            nextMaxPt = thisLegMaxInds(thisPotStartInd + 1);
            % check if next point is max
            if (nextMaxPt > nextMinPt)
                continue;
            end
        end
        
        
        % find end point - maxInds value that immediately follows step mid
        %  point index into thisLegMaxInds of this point
        endPtInd = find(thisLegMaxInds > stepMidPt, 1, 'first');
        % if there is no such point, this full step is not in trial,
        %  continue loop without executing rest
        if isempty(endPtInd)
            continue;
        end
        
        % convert this index into maxInds into frame index of step end
        %  point
        stepEndPt = thisLegMaxInds(endPtInd);
        
        % check that this endpoint is within the same moving bout as the
        %  step start point
        if (stepEndPt < thisMoveStart) || (stepEndPt > thisMoveEnd)
            continue;
        end
        
        
        % if this part of the code is reached, the step start, mid, and end
        %  points are all valid; add step to step arrays
        % step frame indices
        stepInds(counter,:) = [stepStartPt, stepMidPt, stepEndPt];
        % which leg
        stepWhichLeg(counter) = thisLeg;
        
        % update counter
        counter = counter + 1;        
    end
    
    % stepInds and stepWhichLeg shorter than when initialized b/c not all
    %  max/min inds become steps; clear excess zeros
    stepInds = stepInds(1:(counter - 1), :);
    stepWhichLeg = stepWhichLeg(1:(counter - 1));
end
