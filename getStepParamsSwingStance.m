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
%       stepXLengths
%       stepYLengths
%       stepDirections
%       stepDurations
%       stepSpeeds
%       stepVelX
%       stepVelY
%       stepAEPX
%       stepAEPY
%       stepPEPX
%       stepPEPY
%       stepT
%       stepFtFwd
%       stepFtLat
%       stepFtYaw
%       stepSwingStance - n x 2 matrix for swing/stance calls for each half
%           step; 1 for stance, -1 for swing
%       swingStanceMethod
%
% OUTPUTS:
%   stanceStepParams - struct of step parameters during stance, with
%     fields:
%       stepInds
%       stepWhichLeg
%       stepLengths
%       stepXLengths
%       stepYLengths
%       stepDirections
%       stepDurations
%       stepSpeeds
%       stepVelX
%       stepVelY
%       stepAEPX
%       stepAEPY
%       stepPEPX
%       stepPEPY
%       stepT
%       stepFtFwd
%       stepFtLat
%       stepFtYaw
%   swingStepParams - struct of step parameters during swing, with fields:
%       stepInds
%       stepWhichLeg
%       stepLengths
%       stepXLengths
%       stepYLengths
%       stepDirections
%       stepDurations
%       stepSpeeds
%       stepVelX
%       stepVelY
%       stepAEPX - calculated, but AEP and PEP applies really stance only
%       stepAEPY
%       stepPEPX
%       stepPEPY
%       stepT
%       stepFtFwd
%       stepFtLat
%       stepFtYaw
%
% CREATED: 10/2/21 - HHY
%
% UPDATED:
%   10/2/21 - HHY
%   7/6/22 - HHY - instead of manually specifying all legStep parameters,
%       loop through all fields, exempting named ones
%
function [stanceStepParams, swingStepParams] = getStepParamsSwingStance(...
    legSteps)

    % legSteps fields that we're not feeding into
    %  groupStepParamsBySwingStance
    irregularParamNames = {'maxIndsAll', 'minIndsAll', 'maxWhichLeg',...
        'minWhichLeg','userSelVal', 'stepInds', 'stepWhichLeg', ...
        'stepSwingStance', 'swingStanceMethod', 'stepT'};

    % get all field names for legSteps
    legStepsFieldNames = fieldnames(legSteps);

    % loop through all field names, separate parameters by swing and stance
    %  for all applicable parameters
    for i = 1:length(legStepsFieldNames)
        thisName = legStepsFieldNames{i};
        % if this is not in the list of parameters we're not feeding into
        %  groupStepParamsBySwingStance()
        if ~any(strcmp(thisName,irregularParamNames))
            % stance
            stanceStepParams.(thisName) = groupStepParamsBySwingStance(...
                legSteps.(thisName), legSteps.stepSwingStance,1);
            % swing
            swingStepParams.(thisName) = groupStepParamsBySwingStance(...
                legSteps.(thisName), legSteps.stepSwingStance, -1);
        end
    end

    % get stepInds and leg assignments (different array size than other params)
    % stance
    stanceStepParams.stepInds = ...
        legSteps.stepInds(legSteps.stepSwingStance(:,1)==1,1:2);
    stanceStepParams.stepInds = [stanceStepParams.stepInds; ...
        legSteps.stepInds(legSteps.stepSwingStance(:,2)==1,2:3)];

    stanceStepParams.whichLeg = ...
        legSteps.stepWhichLeg(legSteps.stepSwingStance(:,1)==1);
    stanceStepParams.whichLeg = [stanceStepParams.whichLeg; ...
        legSteps.stepWhichLeg(legSteps.stepSwingStance(:,2)==1)];

    stanceStepParams.stepT = ...
        legSteps.stepT(legSteps.stepSwingStance(:,1)==1,1:2);
    stanceStepParams.stepT = [stanceStepParams.stepT; ...
        legSteps.stepT(legSteps.stepSwingStance(:,2)==1,2:3)];
    
    % swing
    swingStepParams.stepInds = ...
        legSteps.stepInds(legSteps.stepSwingStance(:,1)==-1,1:2);
    swingStepParams.stepInds = [swingStepParams.stepInds; ...
        legSteps.stepInds(legSteps.stepSwingStance(:,2)==-1,2:3)];

    swingStepParams.whichLeg = ...
        legSteps.stepWhichLeg(legSteps.stepSwingStance(:,1)==-1);
    swingStepParams.whichLeg = [swingStepParams.whichLeg; ...
        legSteps.stepWhichLeg(legSteps.stepSwingStance(:,2)==-1)];

    swingStepParams.stepT = ...
        legSteps.stepT(legSteps.stepSwingStance(:,1)==-1,1:2);
    swingStepParams.stepT = [swingStepParams.stepT; ...
        legSteps.stepT(legSteps.stepSwingStance(:,2)==-1,2:3)];
end