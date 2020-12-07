% determineStepLength.m
%
% Function that determines length of each step (distance moved in xy plane)
% In units of body lengths
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
%   stepLen - struct of step directions, in degrees; each of below fields
%     has field for each leg
%       back2Front - for steps where leg moves front to back, as specified
%           by stepInd
%       front2Back - for steps where leg moves back to front, as specified
%           by stepInd
%
% CREATED: 11/23/20 - HHY
%
% UPDATED: 11/23/20 - HHY
%
function stepLen = determineStepLength(legXPos, legYPos, stepInd, legInd)

    % get field names
    stepLenNames = fieldnames(stepInd);
    
    % get subfield names
    legNames = fieldnames(stepInd.(stepLenNames{1}));

    % direction for each half step, as defined by start and end pts
    % preallocate
    for i = 1:length(legNames)
        for j = 1:length(stepLenNames)
            stepLen.(stepLenNames{j}).(legNames{i}) = ...
                zeros(size(stepInd.(stepLenNames{j}).(legNames{i}),1),1);
        end
    end

    % compute direction moved for each leg
    for i = 1:length(legNames)
        for j = 1:length(stepLenNames)
        stepLen.(stepLenNames{j}).(legNames{i}) = distBtw2Pts(...
            legXPos(stepInd.(stepLenNames{j}).(legNames{i})(:,1),...
            legInd(i)), ...
            legYPos(stepInd.(stepLenNames{j}).(legNames{i})(:,1),...
            legInd(i)), ...
            legXPos(stepInd.(stepLenNames{j}).(legNames{i})(:,2),...
            legInd(i)), ...
            legYPos(stepInd.(stepLenNames{j}).(legNames{i})(:,2),...
            legInd(i)));
        end
    end

end