% getFrameSwingStance.m
%
% Function to convert step-wise swing stance calls to frame-wise swing
%  stance calls.
% For each frame, for each leg, swing = -1; stance = 1; not moving = 0;
%  no step call = -2.
% For time between end of not moving bout and start of a step, call as
%  stance.
%
% INPUTS:
%   legSteps - struct of leg step data, output from processLegTrack(),
%     fields needed for this function:
%       stepInds
%       stepWhichLeg
%       stepSwingStance - n x 2 matrix for swing/stance calls for each half
%           step; 1 for stance, -1 for swing
%   frameTimes - time stamps for each frame
%   notMoveInd - indices for each not moving frame
%   legIDs - struct of parameters, IDs legs
%       ind - indicies of legs, into raw position matricies, carried
%           throughout
%       names - names of each leg, matching with ind
%
% OUTPUTS:
%   framesSwingStance - frames x #legs matrix of calls
%
% CREATED: 10/3/21 - HHY
%
% UPDATED:
%   10/3/21 - HHY
%
function framesSwingStance = getFrameSwingStance(legSteps, frameTimes, ...
    notMoveInd, legIDs)

    NOT_MOVE_VAL = 0; % value assigned when fly is not moving during frame
    NO_STEP_CALL = -2; % value when there is no step call

    % preallocate (frames x #legs); with -2 for no step call, as others
    %  better defined and this captures remainder
    framesSwingStance = ones(length(frameTimes), length(legIDs.ind))...
        * NO_STEP_CALL;

    % loop through legs
    for i = 1:length(legIDs.ind)
        thisLeg = legIDs.ind(i);

        % this leg's steps
        thisLegStepInds = legSteps.stepInds(...
            legSteps.stepWhichLeg == thisLeg, :);
        % this leg's swing/stance calls
        thisLegSwingStance = legSteps.stepSwingStance(...
            legSteps.stepWhichLeg == thisLeg, :);

        % loop through steps, assign swing/stance values for those indicies
        for j = 1:size(thisLegStepInds,1)        
            % first half step (start to mid)
            framesSwingStance(thisLegStepInds(j,1):thisLegStepInds(j,2), ...
                thisLeg) = thisLegSwingStance(j,1);
            % second half step (mid to end)
            framesSwingStance(thisLegStepInds(j,2):thisLegStepInds(j,3), ...
                thisLeg) = thisLegSwingStance(j,2);
        end

        % assign not moving value
        framesSwingStance(notMoveInd, thisLeg) = NOT_MOVE_VAL;
    end
end