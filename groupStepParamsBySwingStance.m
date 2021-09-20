% groupStepParamsBySwingStance.m
%
% Function that takes in step parameter, swing/stance calls, whether to
%  return swing/stance and returns parameter only for that type of leg
%  movement
%
% INPUTS:
%   stepParam - step parameter values for each half step, as n x 2 matrix
%   stepSwingStance - swing/stance calls for each half step, as n x 2
%       matrix matched to stepParam
%   whichType - scalar value indicating which type to keep (standard is 1
%       for stance, -1 for swing)
%
% OUTPUTS:
%   swingStanceParam - step parameter values for half step, as m x 1 vector
%       where m is the number of half steps that match the type criteria
%
% CREATED: 9/9/21 - HHY
%
% UPDATED:
%   9/9/21 - HHY
%
function swingStanceParam = groupStepParamsBySwingStance(stepParam, ...
    stepSwingStance, whichType)

    % matches in first half step
    swingStanceParam = stepParam(stepSwingStance(:,1)==whichType, 1);
    % matches in second half step
    swingStanceParam = [swingStanceParam; ...
        stepParam(stepSwingStance(:,2)==whichType, 2)];
end