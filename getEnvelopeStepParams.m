% getEnvelopeStepParams.m
%
% Helper function that takes previously extracted AEP and PEP values for
%  each step and uses spline interpolation between them to extract the
%  envelope of the signal, which is an estimate of the continous AEP/PEP
%  value. 
% Then, uses these values to extract continous step length measure
% 
% INPUTS:
%   legTrack - struct of leg tracking data
%   legSteps - struct of computed step parameters
%   moveNotMove - struct of moving/not moving times and indices
%
% OUTPUTS:
%   legStepsCont - struct of output data, with fields:
%     AEPX - numTimePts x numLegs matrix for interpolated AEP X values
%     PEPX - numTimePts x numLegs matrix for interpolated PEP X values
%     AEPY - numTimePts x numLegs matrix for interpolated AEP Y values
%     PEPY - numTimePts x numLegs matrix for interpolated PEP Y values
%     stepLengthX - numTimePts x numLegs matrix for interpolated step
%       length in Y values, as PEPX - AEPX
%     stepLengthY - numTimePts x numLegs matrix for interpolated step
%       length in X values, as PEPY - AEPY
%     stepLength - numTimePts x numLegs matrix for interpolated step
%       length in XY plane, as distance between (AEPX, AEPY) and (PEPX,
%       PEPY) for each time point
%     stepDirection - numTimePts x numLegs matrix for interpolated step
%       direction in XY plane, as angle between (AEPX, AEPY) and (PEPX,
%       PEPY) for each time point
%     t - numTimePts vector for time points of these interpolated values
%       (same as legTrack.t)
%     legIDs - struct of leg info
%       names - names of each leg corresponding to indices (R1-3, L1-3)
%       ind - indices of each leg
%
% CREATED: 8/3/23 - HHY
%
% UPDATED:
%   8/3/23 - HHY
%
function legStepsCont = getEnvelopeStepParams(legTrack, legSteps, ...
    moveNotMove)
    % get leg info (deal with slightly different names for free walk data)
    legIDs = legSteps.legIDs;
    if (isfield(legIDs, 'name'))
        legIDs.names = legIDs.name;
        legIDs = rmfield(legIDs,'name');
    end

    numTimePts = length(legTrack.t);
    numLegs = length(legIDs.ind);

    % initialize output
    legStepsCont.AEPX = zeros(numTimePts, numLegs);
    legStepsCont.PEPX = zeros(numTimePts, numLegs);
    legStepsCont.AEPY = zeros(numTimePts, numLegs);
    legStepsCont.PEPY = zeros(numTimePts, numLegs);
    legStepsCont.stepLengthX = zeros(numTimePts, numLegs);
    legStepsCont.stepLengthY = zeros(numTimePts, numLegs);
    legStepsCont.stepLength = zeros(numTimePts, numLegs);
    legStepsCont.stepDirection = zeros(numTimePts, numLegs);
    legStepsCont.legIDs = legIDs;

    if (isrow(legTrack.t))
        legTrack.t = legTrack.t';
    end

    legStepsCont.t = legTrack.t;

    endTrialFlag = false;

    % loop through all legs
    for i = 1:numLegs
        % logical for steps belonging to this leg
        thisLegLog = legSteps.stepWhichLeg == legIDs.ind(i);
        % this leg step times
        thisLegTimes = legSteps.stepT(thisLegLog,:);

        if (isempty(thisLegTimes))
            endTrialFlag = true;
            break;
        end

        % this leg's AEP step time
        AEPtimes = thisLegTimes(:,2);
        % this leg's PEP step time (account for step being min-max-min)
        PEPtimes = [thisLegTimes(:,1); thisLegTimes(end,3)];

        % this leg's AEPs
        thisAEPX = legSteps.stepAEPX(thisLegLog,1);
        thisAEPY = legSteps.stepAEPY(thisLegLog,1);
        % this leg's PEPs (account for adding last min)
        thisPEPX = [legSteps.stepPEPX(thisLegLog,1); ...
            legSteps.stepPEPX(end,2)];
        thisPEPY = [legSteps.stepPEPY(thisLegLog,1); ...
            legSteps.stepPEPY(end,2)];

        % median for AEPs and PEPs
        medAEPX = median(thisAEPX);
        medAEPY = median(thisAEPY);
        medPEPX = median(thisPEPX);
        medPEPY = median(thisPEPY);

        % replace not moving times with median values - this prevents
        %  spline interpolation from giving wildly 
        % handle not moving times appropriately for free and ball walking
        %  data
        if (isfield(moveNotMove,'notMoveInd'))
            notMoveTimes = legTrack.t(moveNotMove.notMoveInd);
        elseif (isfield(moveNotMove, 'legNotMoveInd'))
            notMoveTimes = legTrack.t(moveNotMove.legNotMoveInd);
        end

        thisAEPX = [thisAEPX; repmat(medAEPX,length(notMoveTimes),1)];
        thisAEPY = [thisAEPY; repmat(medAEPY,length(notMoveTimes),1)];
        thisPEPX = [thisPEPX; repmat(medPEPX,length(notMoveTimes),1)];
        thisPEPY = [thisPEPY; repmat(medPEPY,length(notMoveTimes),1)];

        % times for interpolation
        AEPtimes = [AEPtimes; notMoveTimes];
        PEPtimes = [PEPtimes; notMoveTimes];

        % get repeats in times, remove
        [AEPtimes, AEPrepInd, ~] = unique(AEPtimes,'stable');
        [PEPtimes, PEPrepInd, ~] = unique(PEPtimes, 'stable');

        % remove repeats from AEP and PEP
        thisAEPX(~AEPrepInd) = [];
        thisAEPY(~AEPrepInd) = [];
        thisPEPX(~PEPrepInd) = [];
        thisPEPY(~PEPrepInd) = [];
        
        % get interpolated AEP and PEP for this leg, save to output
        legStepsCont.AEPX(:,legIDs.ind(i)) = ...
            interp1(AEPtimes, thisAEPX, legTrack.t, 'spline', nan);
        legStepsCont.AEPY(:,legIDs.ind(i)) = ...
            interp1(AEPtimes, thisAEPY, legTrack.t, 'spline', nan);
        legStepsCont.PEPX(:,legIDs.ind(i)) = ...
            interp1(PEPtimes, thisPEPX, legTrack.t, 'spline', nan);
        legStepsCont.PEPY(:,legIDs.ind(i)) = ...
            interp1(PEPtimes, thisPEPY, legTrack.t, 'spline', nan);

        % get step length X and Y for this leg
        legStepsCont.stepLengthX(:,legIDs.ind(i)) = ...
            legStepsCont.PEPX(:,legIDs.ind(i)) - ...
            legStepsCont.AEPX(:,legIDs.ind(i));
        legStepsCont.stepLengthY(:,legIDs.ind(i)) = ...
            legStepsCont.PEPY(:,legIDs.ind(i)) - ...
            legStepsCont.AEPY(:,legIDs.ind(i));

        % get step length in XY plane
        legStepsCont.stepLength(:,legIDs.ind(i)) = ...
            distBtw2Pts(legStepsCont.AEPX(:,legIDs.ind(i)),...
            legStepsCont.AEPY(:,legIDs.ind(i)), ...
            legStepsCont.PEPX(:,legIDs.ind(i)),...
            legStepsCont.PEPY(:,legIDs.ind(i)));

        legStepsCont.stepDirection(:,legIDs.ind(i)) = ...
            findAngle2Pts(legStepsCont.AEPX(:,legIDs.ind(i)),...
            legStepsCont.AEPY(:,legIDs.ind(i)), ...
            legStepsCont.PEPX(:,legIDs.ind(i)),...
            legStepsCont.PEPY(:,legIDs.ind(i)));
    end   

    if (endTrialFlag)
        legStepsCont = [];
    end
end