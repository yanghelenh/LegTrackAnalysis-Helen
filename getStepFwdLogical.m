% getStepFwdLogical.m
%
% Function that takes in legSteps struct, a specific leg, and sample times
%  and returns a logical for whether the leg was moving as in forward
%  walking (swing brings the the leg more anterior) or backwards walking
%  (swing brings the leg more posterior). True for forward.
%
% INPUTS:
%   legSteps - struct of leg parameter values, output of
%       extractLegStepsFromPData()
%   sampTime - times at which to return values
%   whichLeg - 'R1', 'R2', 'R3', 'L1', 'L2', 'L3' string specifying which
%       leg
%
% OUTPUTS:
%   outLog - column vector, logical of whether the fly is walking forward.
%     False if fly isn't walking
%
% CREATED: 8/10/23 - HHY
%
% UPDATED:
%   8/10/23 - HHY
%
function outLog = getStepFwdLogical(legSteps, sampTime, whichLeg)
    % get index for this leg
    thisLegInd = legSteps.legIDs.ind(strcmpi(legSteps.legIDs.names, ...
        whichLeg));

    % get swing/stance calls belonging to this leg
    thisSwingStance = legSteps.stepSwingStance(...
        legSteps.stepWhichLeg == thisLegInd, :);
    % get times belonging to this leg
    thisStepT = legSteps.stepT(legSteps.stepWhichLeg == thisLegInd, :);

    % initialize output, false for everything
    outLog = false(size(sampTime));
    % loop through all steps, set appropriate times of outLog to true if
    %  fwd step
    for i = 1:size(thisSwingStance,1)
        % [-1 1] is forward step, [1 -1] is backwards step
        % check based on first column value
        if (thisSwingStance(i,1) == -1) % if fwd
            thisStartTime = thisStepT(i,1);
            thisEndTime = thisStepT(i,3);

            % times for this step, convert to logical into sampTime
            thisFwdLog = (sampTime >= thisStartTime) & ...
                (sampTime <= thisEndTime);

            % update output logical
            outLog(thisFwdLog) = true;
        end
    end
end