% findFlyNotMovingFt.m
%
% Function that takes in fly's smoothed, normalized total speed, a
%  threshold on that speed, and a minimum bout length to extract when the
%  fly is moving or not moving.
% Merges bouts that are too short and assigns 
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
%   notMoveInd - logical for each index of totSpdNormSmo, 1 when fly
%       not moving, 0 when fly moving
%   notMoveStartInd - start indices for when fly not moving
%   notMoveEndInd - end indices for when fly not moving
%
% CREATED: 6/30/22 - HHY
%
% UPDATED:
%   7/1/22 - HHY
%

function [notMoveInd, notMoveStartInd, notMoveEndInd] = ...
    findFlyNotMovingFt(totSpdNormSmo, thresh, minBoutLenSamp)

    % use threshold to get not moving indices
    notMoveLogical = totSpdNormSmo <= thresh;
    
    % indicies where fly transitions between moving and not moving
    transInd = find(diff(notMoveLogical)) + 1;
    
    % add edges to transitions, so first and last bouts are included
    if transInd(end) == length(notMveLogical)
        transInd = [1 transInd];
    else
        transInd = [1 transInd (length(notMoveLogical)+1)];
    end
    
    % duration of each bout of movement and not movement
    boutDur = diff(transInd);
    
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
            
            % index into notMoveLogical for bout transitions
            boutStartMLInd = transInd(boutStartInd);
            % index into transInd for end of bout
            boutEndInd = whichBout +1;
            % equivalent for index into notMoveLogical
            boutEndMLInd = transInd(boutEndInd) - 1;
            
            % is this a moving or not moving bout
            boutMoveLogical = notMoveLogical(boutStartMLInd);
            
            % multiple short bouts, to be merged
            if (whichBout ~= boutStartInd)
                % assign all of these short bouts to 1 longer bout, of the
                %  predominant type (more samples)
                % get number of moving and not moving samples during this
                %  window
                numNotMove = sum(notMoveLogical(boutStartMLInd:boutEndMLInd));
                numMove = sum(~notMoveLogical(boutStartMLInd:boutEndMLInd));

                % if greater than or equal number of not moving samples,
                %  assign all bouts to not moving
                if (numNotMove >= numMove)
                    notMoveLogical(boutStartMLInd:boutEndMLInd) = true;
                else % otherwise, moving
                    notMoveLogical(boutStartMLInd:boutEndMLInd) = false;
                end
            % one short bout, type change   
            else
                notMoveLogical(boutStartMLInd:boutEndMLInd) = ~boutMoveLogical;
            end
        end
        whichBout = whichBout + 1;
    end
    
    % get start and end indices of not move bouts
    notMoveStartInd = find(diff(notMoveLogical) > 0.9) + 1;
    notMoveEndInd = find(diff(notMoveLogical) < -0.9);

    % add start of trial if fly is not moving at start
    if (notMoveLogical(1))
        notMoveStartInd = [1 notMoveStartInd];
    end
    % add end of trial if fly is not moving at end
    if (notMoveLogical(end))
        notMoveEndInd = [notMoveEndInd length(notMoveLogical)];
    end

    % convert notMoveLogical to indices
    notMoveInd = find(notMoveLogical);
end