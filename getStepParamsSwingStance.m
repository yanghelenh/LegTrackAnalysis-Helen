% getStepParamsSwingStance.m
%
% Function to get step parameters during swing and stance phases, instead
%  of during 1st and 2nd half steps
%
% INPUTS:
%   legSteps - struct of leg step data, output from processLegTrack(), with
%     fields:
%       maxIndsAll
%       minIndsAll
%       maxWhichLeg
%       minWhichLeg
%       userSelVal
%       stepInds
%       stepWhichLeg
%       stepLengths
%       stepDirections
%       stepDurations
%       stepSpeeds
%       stepVelX
%       stepVelY
%       stepSwingStance - n x 2 matrix for swing/stance calls for each half
%           step; 1 for stance, -1 for swing
%       swingStanceMethod
%
% OUTPUTS:
%   stanceStepParams - struct of step parameters during stance, with
%     fields:
%       stepLengths
%       stepDirections
%       stepDurations
%       stepSpeeds
%       stepVelX
%       stepVelY
%   swingStepParams - struct of step parameters during swing, with fields:
%       stepLengths
%       stepDirections
%       stepDurations
%       stepSpeeds
%       stepVelX
%       stepVelY
%
% CREATED: 10/2/21 - HHY
%
% UPDATED:
%   10/2/21 - HHY
%
function [stanceStepParams, swingStepParams] = getStepParamsSwingStance(...
    legSteps)

    % stance 
    stanceStepParams.stepLength = groupStepParamsBySwingStance(...
        legSteps.stepLengths, legSteps.stepSwingStance, 1);
    stanceStepParams.stepDirections = groupStepParamsBySwingStance(...
        legSteps.stepDirections, legSteps.stepSwingStance, 1);
    stanceStepParams.stepDurations = groupStepParamsBySwingStance(...
        legSteps.stepDurations, legSteps.stepSwingStance, 1);
    stanceStepParams.stepSpeeds = groupStepParamsBySwingStance(...
        legSteps.stepSpeeds, legSteps.stepSwingStance, 1);
    stanceStepParams.stepVelX = groupStepParamsBySwingStance(...
        legSteps.stepVelX, legSteps.stepSwingStance, 1);
    stanceStepParams.stepVelY = groupStepParamsBySwingStance(...
        legSteps.stepVelY, legSteps.stepSwingStance, 1);
    stanceStepParams.stepFtFwd = groupStepParamsBySwingStance(...
        legSteps.stepFtFwd, legSteps.stepSwingStance, 1);
    stanceStepParams.stepFtLat = groupStepParamsBySwingStance(...
        legSteps.stepFtLat, legSteps.stepSwingStance, 1);
    stanceStepParams.stepFtYaw = groupStepParamsBySwingStance(...
        legSteps.stepFtYaw, legSteps.stepSwingStance, 1);

    % get stepInds and leg assignments (different array size than other params)
    stanceStepParams.stepInds = ...
        legSteps.stepInds(legSteps.stepSwingStance(:,1)==1,1:2);
    stanceStepParams.stepInds = [stanceStepParams.stepInds; ...
        legSteps.stepInds(legSteps.stepSwingStance(:,2)==1,2:3)];

    stanceStepParams.whichLeg = ...
        legSteps.stepWhichLeg(legSteps.stepSwingStance(:,1)==1);
    stanceStepParams.whichLeg = [stanceStepParams.whichLeg; ...
        legSteps.stepWhichLeg(legSteps.stepSwingStance(:,2)==1)];


    % swing
    swingStepParams.stepLength = groupStepParamsBySwingStance(...
        legSteps.stepLengths, legSteps.stepSwingStance, -1);
    swingStepParams.stepDirections = groupStepParamsBySwingStance(...
        legSteps.stepDirections, legSteps.stepSwingStance, -1);
    swingStepParams.stepDurations = groupStepParamsBySwingStance(...
        legSteps.stepDurations, legSteps.stepSwingStance, -1);
    swingStepParams.stepSpeeds = groupStepParamsBySwingStance(...
        legSteps.stepSpeeds, legSteps.stepSwingStance, -1);
    swingStepParams.stepVelX = groupStepParamsBySwingStance(...
        legSteps.stepVelX, legSteps.stepSwingStance, -1);
    swingStepParams.stepVelY = groupStepParamsBySwingStance(...
        legSteps.stepVelY, legSteps.stepSwingStance, -1);
    swingStepParams.stepFtFwd = groupStepParamsBySwingStance(...
        legSteps.stepFtFwd, legSteps.stepSwingStance, -1);
    swingStepParams.stepFtLat = groupStepParamsBySwingStance(...
        legSteps.stepFtLat, legSteps.stepSwingStance, -1);
    swingStepParams.stepFtYaw = groupStepParamsBySwingStance(...
        legSteps.stepFtYaw, legSteps.stepSwingStance, -1);

    
    swingStepParams.stepInds = ...
        legSteps.stepInds(legSteps.stepSwingStance(:,1)==-1,1:2);
    swingStepParams.stepInds = [swingStepParams.stepInds; ...
        legSteps.stepInds(legSteps.stepSwingStance(:,2)==-1,2:3)];

    swingStepParams.whichLeg = ...
        legSteps.stepWhichLeg(legSteps.stepSwingStance(:,1)==-1);
    swingStepParams.whichLeg = [swingStepParams.whichLeg; ...
        legSteps.stepWhichLeg(legSteps.stepSwingStance(:,2)==-1)];
end