% findFlyNotMovingFt.m
%
% Function that takes in fly's smoothed, normalized total speed, a
%  threshold on that speed, and a minimum bout length to extract when the
%  fly is moving or not moving
%
% Helper function for interactGetNotMovingIndWFt()
%
% INPUTS:
%   totSpdNormSmo - smoothed, normalized total speed
%   thresh - threshold on total speed, below which not moving
%   minBoutLenSamp - minimum bout length of moving or not moving bout, in
%       samples
%
% OUTPUTS:
%   notMoveLogical - logical for each index of totSpdNormSmo, 1 when fly
%       not moving, 0 when fly moving
%   notMoveStartInd - start indices for when fly not moving
%   notMoveEndInd - end indices for when fly not moving
%
% CREATED: 6/30/22 - HHY
%
% UPDATED:
%   6/30/22 - HHY
%

function [notMoveLogical, notMoveStartInd, notMoveEndInd] = ...
    findFlyNotMovingFt(totSpdNormSmo, thresh, minBoutLenSamp)


    % use threshold to get moving/not moving bouts
    moveLogical = totSpdNormSmo > thresh;
    
    % indicies where fly transitions between moving and not moving
    transInd = find(diff(moveLogical)) + 1;
    
    % add edges to transitions, so first and last bouts are included
    if transInd(end) == length(moveLogical)
        transInd = [1 transInd];
    else
        transInd = [1 transInd (length(moveLogical)+1)];
    end
    
    % duration of each bout of movement and not movement
    boutDur = diff(transInd);
    
    % cutoff for min bout length, in samples
    minBoutLenSamp = round(minBoutLen * sampRate);
    
    % deal with bouts that are too short (less than minBoutLen) - merge
    % with other, adjacent short bouts or merge into longer sequence
    whichBout = 1;
    while whichBout<=length(boutDur)
        % bout is too short
        if (boutDur(whichBout) < minBoutLenSamp)
            % index into transInd for start of bout
            boutStartInd = whichBout;
            
            % continue until bouts stop being too short
            for k = (whichBout+1):length(boutDur)
                if (boutDur(k) < minBoutLenSamp)
                    whichBout = k;
                else
                    break;
                end
            end
            
            % index into moveLogical for bout transitions
            boutStartMLInd = transInd(boutStartInd);
            % index into tranInd for end of bout
            boutEndInd = whichBout +1;
            % equivalent for index into moveLogical
            boutEndMLInd = transInd(boutEndInd) - 1;
            
            % is this a moving or not moving bout
            boutMoveLogical = moveLogical(boutStartMLInd);
            
            % multiple short bouts, to be merged
            if (whichBout ~= boutStartInd)
                % assign all of these short bouts to 1 longer bout, of the
                % same type as the first bout
                moveLogical(boutStartMLInd:boutEndMLInd) = boutMoveLogical;
            % one short bout, type change   
            else
                moveLogical(boutStartMLInd:boutEndMLInd) = ~boutMoveLogical;
            end
        end
        whichBout = whichBout + 1;
    end
    
    % plot total speed with moving and not moving portions highlighted,
    % with scrolling
    
    % generate patches
    moveStartInd = find(diff(moveLogical) > 0.9) + 1;
    moveStartTimes = t(moveStartInd);
%     moveStartTimes = t((diff(moveLogical) > 0.9) + 1);
    moveEndInd = find(diff(moveLogical) < -0.9);
    moveEndTimes = t(moveEndInd);

    % add start of trial if fly is moving at start
    if (moveLogical(1))
        moveStartTimes = [t(1) moveStartTimes];
        moveStartInd = [1 moveStartInd];
    end
    % add end of trial if fly is moving at end
    if (moveLogical(end))
        moveEndTimes = [moveEndTimes t(end)];
        moveEndInd = [moveEndInd length(moveEndInd)+1];
    end

end