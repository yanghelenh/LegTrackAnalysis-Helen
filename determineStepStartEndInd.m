% determineStepStartEndInd.m
%
% Function that takes in zero crossing indicies (from findVelZeroXing())
%  and returns the indicies of the start and end of each corresponding half
%  step
%
% INPUTS:
%   zeroXing - struct from findVelZeroXing()
%   legInd - indicies of leg tracked pts
%
% OUTPUTS:
%   stepInd - struct of step indicies, with fields:
%       back2Front - has field for each leg, for half steps where leg moves
%           back to front
%       front2Back - for half steps where leg moves front to back
%
% CREATED: 11/18/20 - HHY
%
% UPDATED:
%   11/18/20 - HHY
%   11/23/20 - HHY
%

function stepInd = determineStepStartEndInd(zeroXing, legInd)

    % get field names for zeroXing struct
    legNames = fieldnames(zeroXing.pos2Neg);

    % initialize with empty vector, necessary for appending step indicies
    for i = 1:length(legNames)
        stepInd.front2Back.(legNames{i}) = [];
        stepInd.back2Front.(legNames{i}) = [];
    end

    % neg2pos is front-most position; pos2neg is back-most position
    for i = 1:length(legInd)
        for j = 1:length(zeroXing.pos2Neg.(legNames{i}))
            curInd = zeroXing.pos2Neg.(legNames{i})(j);
            midInd = zeroXing.neg2Pos.(legNames{i})(find(...
                zeroXing.neg2Pos.(legNames{i}) > curInd, 1, 'first'));
            endInd = zeroXing.pos2Neg.(legNames{i})(find(...
                zeroXing.pos2Neg.(legNames{i}) > curInd, 1, 'first'))...
                - 1;

            % save into struct
            if (~isempty(midInd) && ~isempty(endInd))
                stepInd.back2Front.(legNames{i}) = ...
                    [stepInd.back2Front.(legNames{i}); ...
                    curInd (midInd-1)];
                stepInd.front2Back.(legNames{i}) = ...
                    [stepInd.front2Back.(legNames{i}); ...
                    midInd endInd];
            end
        end
    end

end