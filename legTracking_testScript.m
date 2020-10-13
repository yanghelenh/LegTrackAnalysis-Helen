% legTracking_testScript.m
%
% Script for testing leg tracking processing
% TEMPORARY
%
% CREATED: 10/12/20 - HHY
%
% UPDATED: 
%   10/12/20 - HHY
%

% file paths
trialPath = ...
    '/Volumes/Samsung_2tb/Wilson Lab/RAW DATA/190507/fly01/fov02/trial01';

legTrackPath = ...
    '/Users/hyang/Dropbox (HMS)/APTlegTracking';
legTrackFile = '190507_fly01_fov02_trial01_legVid_caLegTrack_v1_mdn.trk';


curDir = pwd; % remember current directory

% load in leg tracking data from .trk file; need to load as .mat
cd(legTrackPath);
load(legTrackFile, '-mat', 'pTrk');

% load in leg video data, FicTrac data
cd(trialPath);
load('legVidDat.mat', 'legVidFrameTimes');
load('fictracDat.mat');
ficTracTimes = t;
