% findBoutStepIndsBall.m
%
% Helper function for saveBallLegStepParamCond_bouts() that takes in bout 
%  peak, start, and end times (output of findCondYawVelPeaksFT() helper 
%  function) as well as step index info, and returns aligned steps, as 
%  indices into legStep matrices.
% Also, returns whether bout peak is during stance or swing for each leg
%
% Adaptation of findBoutStepInds(), which runs on free-walking data. Input
%  is times and not indices because leg and FicTrac are not on the same
%  time scale, unlike in free-walking data
%
% INPUTS:
%   legSteps - struct of legStep info
%   peakTimes - times of bout peaks
%   boutStartTimes - times of bout starts
%   boutEndTimes - indices of bout ends
%   legT - timing vector for leg frames, for this whole trial
%   maxNumSteps - number of steps to each side of peak to consider as part
%       of bout (max bout length is this x2 + 1)
%   legIDs - struct defining indices corresponding to each leg
%
% OUTPUTS:
%   boutStepInd - maxNumSteps*2+1 x numLegs x numBouts size matrix defining
%       aligned bouts as indices into legSteps matrices
%   boutPkSwingStance - numLegs x numBouts size matrix indicating whether
%       the peak of each bout is during swing or stance for each leg (-1
%       for swing, +1 for stance)
%   rmInd - indices of bouts removed
%
% CREATED: 6/22/23 - HHY
%
% UPDATED:
%   6/22/23 - HHY
%   8/7/23 - HHY - output removed bouts ind
%
function [boutStepInd, boutPkSwingStance, rmInd] = findBoutStepIndsBall(legSteps, ...
    peakTimes, boutStartTimes, boutEndTimes, legT, maxNumSteps, legIDs)

    % convert bout peak, start, and end times into indices on leg timescale
    % preallocate
    peakInd = zeros(size(peakTimes));
    boutStartInd = zeros(size(boutStartTimes));
    boutEndInd = zeros(size(boutEndTimes));

    % to keep track of indices to remove: we can't find a corresponding
    %  legT index for any of the 3 parameters associated with each bout
    rmInd = [];

    % get closest leg index for peaks
    for i = 1:length(peakTimes)
        thisPeakInd = find(legT <= peakTimes(i), 1, 'last');
        if (~isempty(thisPeakInd))
            peakInd(i) = thisPeakInd;
        else
            rmInd = [rmInd; i];
        end
    end

    % get closest leg index for starts
    for i = 1:length(boutStartTimes)
        thisStartInd = find(legT <= boutStartTimes(i), 1, 'last');
        if (~isempty(thisStartInd))
            boutStartInd(i) = thisStartInd;
        else
            rmInd = [rmInd; i];
        end
    end

    % get closest leg index for starts
    for i = 1:length(boutEndTimes)
        thisEndInd = find(legT <= boutEndTimes(i), 1, 'last');
        if (~isempty(thisEndInd))
            boutEndInd(i) = thisEndInd;
        else
            rmInd = [rmInd; i];
        end
    end

    % remove any bouts where we're missing leg time indices
    rmInd = unique(rmInd); % no repeats
    % remove these bouts
    peakInd(rmInd) = [];
    boutStartInd(rmInd) = [];
    boutEndInd(rmInd) = [];

    % some params
    numLegs = length(legIDs.ind);
    numBouts = length(peakInd);

    % preallocate
    boutStepInd = nan(maxNumSteps * 2 + 1, numLegs, numBouts);
    boutPkSwingStance = nan(numLegs, numBouts);

    % loop through all peaks
    for i = 1:length(peakInd)
        thisPeakInd = peakInd(i); 

        % loop through all legs
        for j = 1:numLegs
            % logical for steps for this leg only
            thisLegLog = legSteps.stepWhichLeg == legIDs.ind(j);

            % get index of step during which peak occurs
            % if peak is on step boundary, takes later step (based on which
            %  one is equal)
            % also ensure it's for this leg only
            stepAtPeak = find((legSteps.stepInds(:,1)<=thisPeakInd) & ...
                (legSteps.stepInds(:,3)>thisPeakInd) & thisLegLog);

            % check that a step is present
            % if a step is not present, skip to next leg
            if (isempty(stepAtPeak))
                continue;
            end

            % determine whether peak falls during step's swing or stance
            % if the peak is during the first half
            if (thisPeakInd < legSteps.stepInds(stepAtPeak, 2))
                thisSwSt = legSteps.stepSwingStance(stepAtPeak, 1);
            % otherwise, in second half    
            else
                thisSwSt = legSteps.stepSwingStance(stepAtPeak, 2);
            end
            boutPkSwingStance(j,i) = thisSwSt;

            % place step at peak into output matrix
            boutStepInd(maxNumSteps + 1, j, i) = stepAtPeak;

            % get steps prior to peak step, place in output matrix
            for k = 1:maxNumSteps
                % adjacent steps should have adjacent indices in stepInds
                %  and end of one step should equal start of another in
                %  frame indices
                adjStepInd = stepAtPeak - k;

                % if the step index is before trial start, end update
                if (adjStepInd < 1)
                    break;
                % check that adjacent step still belongs to same leg    
                elseif (legSteps.stepWhichLeg(adjStepInd)~=legIDs.ind(j))
                    break;
                % check that end frame of this step matches start of next one
                % if not, end update
                elseif (legSteps.stepInds(adjStepInd,3)~=...
                        legSteps.stepInds(adjStepInd+1,1))
                    break;
                % check that the start of this step happens after the start 
                %  the bout
                elseif (legSteps.stepInds(adjStepInd,1) > ...
                        boutStartInd(i))
                    break;
                % if valid, add to matrix    
                else
                    boutStepInd(maxNumSteps + 1 - k,j,i) = adjStepInd;
                end
            end

            % get steps after peak step, place in output matrix
            for k = 1:maxNumSteps
                % adjacent steps should have adjacent indices in stepInds
                %  and end of one step should equal start of another in
                %  frame indices
                adjStepInd = stepAtPeak + k;

                % if the step index is after trial end, end update
                if (adjStepInd > length(legSteps.stepWhichLeg))
                    break;
                % check that adjacent step still belongs to same leg 
                elseif (legSteps.stepWhichLeg(adjStepInd)~=legIDs.ind(j))
                    break;
                % check that end of last frame matches start of this one
                % otherwise, end update
                elseif (legSteps.stepInds(adjStepInd,1) ~= ...
                        legSteps.stepInds(adjStepInd - 1,3))
                    break;
                % check that the end of this step happens before the end of
                %  the bout
                elseif (legSteps.stepInds(adjStepInd,3) < ...
                        boutEndInd(i))
                    break;
                % if valid, add to matrix
                else
                    boutStepInd(maxNumSteps + 1 + k, j, i) = adjStepInd;
                end
            end

        end
    end
end