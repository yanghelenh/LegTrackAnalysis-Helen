% getLegPhaseFromSteps.m
%
% Function that takes in leg steps and converts them into leg phase.
%
% For each step (max-min-max), first max is 0 deg, min is 180 deg, second
%  max is 360 deg. X positions in the first half step (between the first
%  max and the min) map linearly to 0-180 deg. X positions in the second
%  half step (between the min and the second max) map linearly to 180-360
%  deg.
% For not moving, NaN
% Where there are gaps between steps (i.e. second max of step is not first
%  max of next step, for those frames, start point is second max of
%  previous step and end point is first max of next step. If start > end,
%  map linearly 0 to 180 deg. If end > start, map linearly 180 to 360 deg.
%  Same applies when there's a gap bewteen the end of a not moving bout and
%  the start of a step
%
% INPUTS:
%   legStepInds - n x 3 matrix of steps for all legs, each row is
%     max-min-max indices for one step
%   stepWhichLeg - n x 1 vector for which leg each step belongs to
%   legXPos - m x 6 matrix of leg X positions, where m = number of frames
%       in trial
%   notMoveBout - o x 2 matrix of start (col 1) and end (col 2) indices of
%       not moving bouts
%
% OUTPUTS:
%   legPhase - m x 6 matrix of leg phases, 1 column for each leg
%
% CREATED: 10/12/21 - HHY
%
% UPDATED:
%   10/22/21 - HHY
%   10/26/21 - HHY - add correction for when calculated phase <0 or >360
%   6/26/23 - HHY - bug fix for correction of phase <0 or >360
%
function legPhase = getLegPhaseFromSteps(legStepInds, stepWhichLeg, ...
    legXPos, notMoveBout)

    % some constants
    NOT_MOVE_VAL = NaN;
    START_VAL = 0; % step start phase, in degrees
    MID_VAL = 180; % step mid phase, in degrees
    END_VAL = 360; % step end phase, in degrees
    
    % preallocate legPhase output matrix
    legPhase = zeros(size(legXPos));
    
     
    numLegs = size(legXPos, 2); % number of legs
    numFrames = size(legXPos, 1); % number of frames in trial
    
    % loop through all legs
    for i = 1:numLegs
        
        % get steps corresponding to this leg only
        thisLegStepInds = legStepInds(stepWhichLeg == i, :);
        
        % number of steps for this leg
        numSteps = size(thisLegStepInds, 1);
        
        % track last index assigned a phase, to appropriately handle gaps
        lastInd = 0; % initialize at 0, start index
        
        % loop through all steps
        for j = 1:numSteps
            stepStartInd = thisLegStepInds(j,1);
            
            % check if last index is same as this step's start index (true
            %  if this step was followed by a previous step)
            % if there is a gap
            if (lastInd ~= stepStartInd)
                % get start and end of gap
                gapStart = lastInd + 1;
                gapEnd = stepStartInd - 1;
                
                % get the indices of the gap
                gapInds = gapStart:gapEnd;
                
                % check if this gap has a not moving portion within it
                nmbChkInd = find(notMoveBout(:,1) >= gapStart, 1, 'first');

                % end index of possible not move bout within gap
                nmbEndInd = notMoveBout(nmbChkInd, 2);
                
                % if end index of not move bout is within gap
                if (nmbEndInd <= gapEnd)
                    % some portion of gap is accounted for by not moving
                    %  bout, consider leg not moving during whole duration
                    %  of gap (since legs start moving at different times
                    %  on move start and not move bout call is approximate)
                    legPhase(gapInds, i) = NOT_MOVE_VAL;
                
                % if the gap is at the start of a trial but doesn't contain
                %  a not moving bout, it must be before the leg starts
                %  moving for that moving bout - consider leg as not moving
                elseif (lastInd == 0)
                    legPhase(gapInds, i) = NOT_MOVE_VAL;
                    
                % if not, gap is within a moving bout (missed step call)
                % treat like single missed step
                else
                    % find minimum position within gap
                    [~, gapMinIdx] = min(legXPos(gapInds,i));
                    % convert to index into trial
                    thisMinInd = gapMinIdx + gapStart - 1;
                    
                    % get phases during gap
                    
                    % for gap start ind to minInd
                    startPos = legXPos(gapStart, i);
                    midPos = legXPos(thisMinInd, i);
                    % vector of leg positions during this part of gap
                    thisGapLegPos = legXPos(gapStart:thisMinInd, i);
                    
                    % leg phases for first part of gap
                    phaseStart2Mid = get1stHalfStepPhase(startPos, ...
                        midPos, thisGapLegPos);
                    % save these into leg phase vector
                    legPhase(gapStart:thisMinInd, i) = phaseStart2Mid;
                    
                    % for minInd to gap end ind
                    endPos = legXPos(gapEnd, i);
                    % vector of leg positions during this part of gap
                    thisGapLegPos = legXPos(thisMinInd:gapEnd, i);
                    
                    % leg phases for second part of gap
                    phaseMid2End = get2ndHalfStepPhase(midPos, endPos, ...
                        thisGapLegPos);
                    % save these into leg phase vector
                    legPhase(thisMinInd:gapEnd, i) = phaseMid2End;   
                end
            end
            
            % get phases for this step
            
            % step indices
            stepMidInd = thisLegStepInds(j,2); 
            stepEndInd = thisLegStepInds(j,3);
            
            % leg position at these indices
            stepStartPos = legXPos(stepStartInd, i);
            stepMidPos = legXPos(stepMidInd, i);
            stepEndPos = legXPos(stepEndInd, i);
            
            % this step's leg positions
            stepPosStart2Mid = legXPos(stepStartInd:stepMidInd, i);
            stepPosMid2End = legXPos(stepMidInd:stepEndInd, i);
            
            % get phases for this step
            phaseStart2Mid = get1stHalfStepPhase(stepStartPos, ...
                stepMidPos, stepPosStart2Mid);
            phaseMid2End = get2ndHalfStepPhase(stepMidPos, stepEndPos, ...
                stepPosMid2End);
            
            % save these phases into leg phase vector
            legPhase(stepStartInd:stepMidInd, i) = phaseStart2Mid;
            legPhase(stepMidInd:stepEndInd, i) = phaseMid2End;
            
            
            % update lastInd counter
            lastInd = stepEndInd;
        end
        
        % fill in any gaps between end of last step and end of trial
        % treat as not moving
        if (lastInd ~= numFrames)
            % start and end indices of end gap
            endGapStartInd = lastInd + 1;
            endGapEndInd = numFrames;
            
            legPhase(endGapStartInd:endGapEndInd, i) = NOT_MOVE_VAL;
        end
        
        
        % correct for any phase values < 0 or > 360
        
        % find indices where phase invalid
        invalInd = find((legPhase(:,i) < START_VAL) | ...
            (legPhase(:,i) > END_VAL));

        if ~isempty(invalInd)
            
            % convert indices into bouts (find times when single invalid pt vs
            %  many in a row)
            % returned indices are into legPhase
            [invalStartInd, invalEndInd, invalDur] = findBouts(invalInd);
            
            % for all times when it's a single invalid point, correct by
            %  converting to mean of +/- numMeanPts points
            % for times when it's multiple invalid points in a row, correct by
            %  converting to START_VAL or END_VAL
            % loop through all invalid bouts
            numMeanPts = 2; % use this number of points before and after
            for k = 1:length(invalDur)
                % single invalid point
                if (invalDur(k) == 0)
                    % get points before and after this point
                    bfInds = ...
                        (invalStartInd(k)-numMeanPts):(invalStartInd(k)-1);
                    aftInds = ...
                        (invalStartInd(k)+1):(invalStartInd(k)+numMeanPts);
                    
                    % check that before and after inds are valid
                    if (any(bfInds < 1))
                        bfInds = [];
                    end
                    
                    if(any(aftInds > numFrames))
                        aftInds = [];
                    end
                    
                    % get values of points before and after
                    bfVals = legPhase(bfInds,i);
                    aftVals = legPhase(aftInds,i);
                    
                    % check that phases of before and after points are valid
                    %  (not NaN)
                    if (any(isnan(bfVals)))
                        bfVals = [];
                    end
                    
                    if(any(isnan(aftVals)))
                        aftVals = [];
                    end
                    
                    % get mean phase of points
                    %  if empty vector, will be NaN (keep that)
                    meanPtVal = mean([bfVals; aftVals]);
                    
                    % replace invalid point with mean val
                    legPhase(invalStartInd(k),i) = meanPtVal;
                
                % multiple invalid points in a row
                else
                    % get mean value of invalid points in bout
                    meanInvalVal = mean(...
                        legPhase(invalStartInd(k):invalEndInd(k),i));
                    
                    % replace these invalid points with START_VAL or END_VAL,
                    %  depending on which direction they're invalid
                    if (meanInvalVal > END_VAL)
                        legPhase(invalStartInd(k):invalEndInd(k),i) = ...
                            repmat(END_VAL, invalDur(k) + 1, 1);
                    elseif (meanInvalVal < START_VAL)
                        legPhase(invalStartInd(k):invalEndInd(k),i) = ...
                            repmat(START_VAL, invalDur(k) + 1, 1);
                    else % should never get here, but NaN if weird
                        legPhase(invalStartInd(k):invalEndInd(k),i) = ...
                            NaN(invalDur(k) + 1, 1);
                    end
                end    
            end
        end
        
        
    end  
end