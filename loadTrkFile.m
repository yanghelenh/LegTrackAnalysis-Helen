% loadTrkFile.m
%
% Function to load Animal Part Tracker (APT) .trk output file. Returns leg
%  X and Y coordinates, in pixels.
% For current orientation of fly, X is the axis parallel to the long axis
%  of the fly. Y is axis perpendicular to long axis of fly.
%
% INPUT:
%   trkFilePath - full path to .trk file
%
% OUTPUT:
%   legX - a number of frames x number of points matrix of X-coordinates
%   legY - a number of frames x number of points matrix of Y-coordinates
%
% CREATED: 11/16/20 - HHY
%
% UPDATED: 
%   11/16/20 - HHY
%
function [legX, legY] = loadTrkFile(trkFilePath)

    % need to load .trk file as .mat file; pTrk is the data we want
    load(trkFilePath, '-mat', 'pTrk');
    
    % get x and y coordinates separately
    legX = squeeze(pTrk(:,1,:))';
    legY = squeeze(pTrk(:,2,:))';
end