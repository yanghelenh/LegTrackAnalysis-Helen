% groupStepParamsBySwingStance.m
%
% Function that takes in step parameter, swing/stance calls, whether to
%  return swing/stance and returns parameter only for that type of leg
%  movement
%
% INPUTS:
%   stepParam - step parameter values for each half step, as n x 1, 2, or 3
%       matrix
%   stepSwingStance - swing/stance calls for each half step, as n x 1, 2,
%       or 3 matrix matched to stepParam
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
%   10/4/22 - HHY - was not working correctly. re-wrote
%   10/11/22 - HHY - modify to work on all step params, regardless of
%       number of elements (1, 2, or 3)
%
function swingStanceParam = groupStepParamsBySwingStance(stepParam, ...
    stepSwingStance, whichType)

%     % matches in first half step
%     swingStanceParam = stepParam(stepSwingStance(:,1)==whichType, 1);
%     % matches in second half step
%     swingStanceParam = [swingStanceParam; ...
%         stepParam(stepSwingStance(:,2)==whichType, 2)];

    % different options depending on whether stepParam has 1, 2 or 3
    %  elements
    switch (size(stepParam,2))
        % if only 1 option, same regardless of swing or stance
        case 1 
            swingStanceParam = stepParam;
        % 2 elements, one for each half step, chose appropriate one
        case 2
            swingStanceParam = zeros(size(stepParam,1),1);
            for i =1:size(stepParam,1)
                swingStanceParam(i) = ...
                    stepParam(i,stepSwingStance(i,:) == whichType);
            end
        % 3 elements mark step start, mid, end; grab 2 adjacent points for
        %  start and end for appropriate half step
        case 3
            swingStanceParam = zeros(size(stepParam,1),2);
            for i =1:size(stepParam,1)
                % whether swing/stance match for 1st half step or 2nd
                switch (find(stepSwingStance(i,:) == whichType))
                    % 1st half step, take 1st 2 points
                    case 1
                        swingStanceParam(i,:) = stepParam(i,1:2);
                    % 2nd half step, take 2nd 2 points
                    case 2
                        swingStanceParam(i,:) = stepParam(i,2:3);
                end
            end 
    end
end