% preprocessLegTrack.m
%
% Function that takes in full path to .trk file and returns preprocessed
%  leg tracking data (leg position and velocity) in struct.
% .trk file is output from APT
%
% NOTE: if I decide to add filtering on leg position here, ala Anipose,
%  this is where it should be added
%
% INPUTS:
%   trkFilepath - full path to .trk file 
%   frameTimes - time at which each leg video frame occurred
%   refPts - struct of parameters defining reference tracked points
%   smoParams - struct of parameters for smoothing leg velocity
%
% OUTPUTS:
%   legTrack - struct with processed data, fields:
%       trkFilepath - full path to .trk file
%       trkFilename - name of .trk file
%       legX - raw leg positions in X (front-back axis), straight from .trk
%       legY - raw leg positions in Y (front-back axis), straight from .trk
%       numFrames - number of leg tracking frames
%       numPoints - number of points being tracked
%       t - time vector for leg position/velocity data, time of each leg
%           video frame
%       srnLegX - leg positions in X after alignment and normalization
%       srnLegY - leg positions in Y after alignment and normalization
%       legXVel - instantaneous leg velocity in X, with Gaussian process
%           smoothing
%       legYVel - instantaneous leg velocity in Y, with Gaussian process
%           smoothing
%       refPts - struct of parameters defining reference tracked points
%       smoParams - struct of parameters for smoothing leg velocity
%
% CREATED: 9/24/21 - HHY
%
% UPDATED:
%   9/24/21 - HHY
%
function legTrack = preprocessLegTrack(trkFilepath, frameTimes, ...
    refPts, smoParams)

    % load .trk file; get leg positions
    [legTrack.legX, legTrack.legY] = loadTrkFile(trkFilepath);

    % some parameters about .trk file
    % number of frames in this trial
    legTrack.numFrames = size(legTrack.legX, 1); 
    % number of tracked points
    legTrack.numPoints = size(legTrack.legX, 2);
    
    % frame times
    legTrack.t = frameTimes;

    % get leg positions after alignment to fly midpoint and normalization 
    %  to body length units
    [legTrack.srnLegX, legTrack.srnLegY] = shiftRotateNormalizeLegPos(...
        legX, legY, refPts);

    % get leg velocities, with Gaussian process smoothing on position
    legTrack.legXVel = findLegVel(legTrack.srnLegX, smoParams);
    legTrack.legYVel = findLegVel(legTrack.srnLegY, smoParams);
    
    % copy over parameter structs
    legTrack.refPts = refPts;
    legTrack.smoParams = smoParams;

end