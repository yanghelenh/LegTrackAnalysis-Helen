% rmvIndAtSelX.m
%
% Function that takes in specific x-coordinate value and removes
%  the corresponding index from the given indices vector
% Called by interactGetLegReversals(), meant to remove selected point
%  from max/min indices list
% Converts x-coordinate value to index by using input x-coords to index
%  mapping vector, then finds nearest index in indices vector, and
%  removes it
%
% INPUTS:
%   inds - vector of indices to add/remove from
%   xCoords - vector of x coordinates, allows mapping between indices
%       (indices of this vector match indices of inds)
%   thisX - x coordinate of point to add/remove
%
% OUTPUTS:
%   thisInd - the index corresponding to thisX that was added/removed
%   inds - updated inds vector, with thisInd added/removed
%
% CREATED: 9/30/21 - HHY
%
% UPDATED:
%   9/30/21 - HHY
%
function [thisInd, inds] = rmvIndAtSelX(inds, xCoords, thisX) 

    % get index closest to thisX, within full set of indices of trial
    [~, closeInd] = min(abs(xCoords - thisX));
    
    % get index closest to closeInd, within list of indices in inds
    [~, ptToRmv] = min(abs(inds - closeInd));
    
    thisInd = inds(ptToRmv);
    
    % remove thisInd from inds
    inds(ptToRmv) = [];
end