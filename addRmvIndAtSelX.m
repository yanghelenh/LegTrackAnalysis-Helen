% addRmvIndAtSelX.m
%
% Function that takes in specific x-coordinate value and adds or removes
%  the corresponding index from the given indices vector
% Called by interactGetLegReversals(), meant to add/remove selected point
%  from max/min indices list
% Converts x-coordinate value to index by using input x-coords to index
%  mapping vector, taking window around that index and finding max/min
%  within that window
% Adds/removes nearest value from indices list
%
% INPUTS:
%   inds - vector of indices to add/remove from
%   xCoords - vector of x coordinates, allows mapping between indices
%       (indices of this vector match indices of inds)
%   thisX - x coordinate of point to add/remove
%   addRmvStr - string 'add' or 'remove' for whether to add or remove pt
%   maxMinStr - string 'max' or 'min' for whether point is max or min
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
function [thisInd, inds] = addRmvIndAtSelX(inds, xCoords, thisX, ...
    addRmvStr, maxMinStr)

end