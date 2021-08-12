% findBouts.m
%
% Function that returns start and end of bouts given indicies of all
%  elements that belong to the bouts. Also returns duration in samples of
%  all bouts
%
% INPUTS:
%   inBoutInd - vector listing all indicies that are part of the bout
%
% OUTPUTS:
%   boutStartInd - vector of start indicies of all the bouts
%   boutEndInd - vector of end indicies of all the bouts
%   boutDur - vector of duration of each bout (in samples)
%
% CREATED: 10/26/20 - HHY
%
% UPDATED:
%   10/26/20 - HHY
%
function [boutStartInd, boutEndInd, boutDur] = findBouts(inBoutInd)

    % make sure inBoutInd is a column vector
    if ~iscolumn(inBoutInd)
        inBoutInd = inBoutInd';
    end

    boutIndChange = diff(inBoutInd); % difference between indicies

    % when difference greater than 1, start of new bout
%     boutStartInd = inBoutInd(find(boutIndChange > 1) + 1);
    
    boutStartInd = inBoutInd(find(boutIndChange > 1) + 1);
    % append first index for start of first bout
    boutStartInd = [inBoutInd(1); boutStartInd]; 
    
    % end of bout; offset by 1 unit from starts (no +1)
    boutEndInd = inBoutInd(find(boutIndChange > 1)); 
    % append end value for end of last bout
    boutEndInd = [boutEndInd; inBoutInd(end)]; 
    
    % compute bout duration in samples by subtracting starts and ends
    boutDur = boutEndInd - boutStartInd;
end