% mergeAdjVecVals.m
%
% Function that takes in vector with a sequence of numbers that are
%  sometimes consecutive and sometimes not and returns the values in the
%  original vector at the start point and end point of those sequences and
%  the corresponding indicies into the original vector.
%
% INPUTS:
%   seqVec - vector that contains numerical sequence
%   maxDiff - maximum difference between adjacent indicies
%
% OUTPUTS:
%   startptVal - values (from seqVec) of start points of sequences
%   startptInd - indicies into seqVec of start points of sequences
%   endptVal - values (from seqVec) of end points of sequences
%   endptInd - indicies into seqVec of end points of sequences
%
% CREATED: 8/9/21 - HHY
%
% UPDATED:
%   8/9/21 - HHY
%   8/10/21 - HHY - fix bug in how start points of sequences being defined
%       (needs to be relative to original vector)
%   8/11/21 - HHY - fix bug in how start points of sequences being
%       calculated when maxDiff > 1 means some sequences are merged
%   3/4/23 - HHY - really strange bug where empty values passed to this
%       helper function; deal with empty seqVec
%
function [startptVal, startptInd, endptVal, endptInd] = mergeAdjVecVals(...
    seqVec, maxDiff)

    if ~isempty(seqVec)
        % check if seqVec is row vector; if yes, make into column and remember
        if (isrow(seqVec))
            seqVec = seqVec';
            undoTrans = true;
        else
            undoTrans = false;
        end
        
        % find endpts of sequences
        diffVec=diff(seqVec); % difference between elements
        seqInd=find([diffVec; inf]>maxDiff); % ind within sequence
        seqLen=diff([0; seqInd]); % length of the sequences
        endptInd=cumsum(seqLen); % endpoints of the sequences
        
        % start point values (seqLen is defined in indicies of seqVec)
        %  fixed 8/11/21
        startptVal = seqVec(endptInd - seqLen + 1);
        
        % end point values
        endptVal = seqVec(endptInd);
        
        % get start point indicies from start point values
        % preallocate
        startptInd = zeros(size(startptVal));
        % loop through all startptVals, find corresponding index into secVec
        for i = 1:length(startptVal)
            % possible indicies, when equal
            possInd = find(seqVec == startptVal(i));
            if (i > 1)
                % indicies in order, so actual one is next possible one
                startptInd(i) = possInd(...
                    find(possInd > startptInd(i-1), 1, 'first'));
            else
                % has to be first if first index of startptVal
                startptInd(i) = possInd(1);
            end
        end
        
        % if it was row vector, convert back into row vector
        if undoTrans
            endptInd = endptInd';
            endptVal = endptVal';
            startptInd = startptInd';
            startptVal = startptVal';
        end
    else
        startptVal = [];
        startptInd = [];
        endptVal = [];
        endptInd = [];
    end
end