% determineStepDirection.m
%
% Function that determines the angle moved in the xy plane for each leg
%  during each step.
% 0 deg points straight up toward head; 90 deg points right; -90 points
%  left; 180/-180 points down abdomen
%
% INPUTS:
%   legXPos - matrix of leg x positions over time; num frames x num tracked
%       pts
%   legYPos- matrix of leg y positions over time; num frames x num tracked
%       pts
%   stepInd - struct output of determineStepStartEndInd, specifies start
%       and end indicies of each half step
%   legInd - indicies of leg tracked pts
%
% OUTPUTS:
%   stepDir - struct of step directions, in degrees; each of below fields
%     has field for each leg
%       back2Front - for steps where leg moves front to back, as specified
%           by stepInd
%       front2Back - for steps where leg moves back to front, as specified
%           by stepInd
%
% CREATED: 11/19/20 - HHY
%
% UPDATED: 11/19/20 - HHY
%
function stepDir = determineStepDirection(legXPos, legYPos, stepInd, ...
    legInd)

    % get field names
    stepDirNames = fieldnames(stepInd);
    
    % get subfield names
    legNames = fieldnames(stepInd.(stepDirNames{1}));

    % direction for each half step, as defined by start and end pts
    % preallocate
    for i = 1:length(legNames)
        for j = 1:length(stepDirNames)
            stepDir.(stepDirNames{j}).(legNames{i}) = ...
                zeros(size(stepInd.(stepDirNames{j}).(legNames{i}),1),1);
        end
    end

    % compute direction moved for each leg
    for i = 1:length(legNames)
        for j = 1:length(stepDirNames)
        stepDir.(stepDirNames{j}).(legNames{i}) = findAngle2Pts(...
            legXPos(stepInd.(stepDirNames{j}).(legNames{i})(:,1),...
            legInd(i)), ...
            legYPos(stepInd.(stepDirNames{j}).(legNames{i})(:,1),...
            legInd(i)), ...
            legXPos(stepInd.(stepDirNames{j}).(legNames{i})(:,2),...
            legInd(i)), ...
            legYPos(stepInd.(stepDirNames{j}).(legNames{i})(:,2),...
            legInd(i)));
        end
    end

end