% legTrackEphysAnalysisScript.m
%
% Quick and dirty script for starting to look at leg tracking and ephys
%  responses together

%% analysis constants
% indicies for specific body parts
% fit line to head and 3 pts defining midpt b/w legs
refPts.lineFitInd = 7:10; 
refPts.midPtInd = 9;
refPts.headPtInd = 7;
refPts.abdPtInd = 11;

% parameters for determining whether fly is moving
notMoveParams.padLen = int32(50); % pad length in samples
notMoveParams.sigma = int32(10); % in samples
notMoveParams.medFiltNumSamps = 10;
notMoveParams.zeroVelThreshMed = 0.004;
notMoveParams.zeroVelThresh = 0.01;
notMoveParams.movePosVelThresh = 0.01;
notMoveParams.moveNegVelThresh = -0.02;
notMoveParams.minBoutLen = 5; % in samples
notMoveParams.r2LegInd = 2;
notMoveParams.l2LegInd = 5;

%% 200821_fly01_cell01_trial01
% full path to .trk leg tracking data
trkFilePath = ['/Users/hyang/Dropbox (HMS)/APTlegTracking/'...
    'Tracking_v1_mdn_labeled2798_ephys/' ...
    '200821_fly01_cell01_trial01_legVid_crop_caLegTrack_v1_mdn.trk'];
pDataFilePath = ['/Users/hyang/Dropbox (HMS)/EphysAnalysis-Helen/pData/'...
    '200821_fly01_cell01_trial01.mat'];

%% 200826_fly01_cell01_trial01
% full path to .trk leg tracking data
trkFilePath = ['/Users/hyang/Dropbox (HMS)/APTlegTracking/'...
    'Tracking_v1_mdn_labeled2798_ephys/' ...
    '200826_fly01_cell01_trial01_legVid_crop_caLegTrack_v1_mdn.trk'];
pDataFilePath = ['/Users/hyang/Dropbox (HMS)/EphysAnalysis-Helen/pData/'...
    '200826_fly01_cell01_trial01.mat'];

%% pre-process leg tracking data

% load .trk file; get leg positions
[legX, legY] = loadTrkFile(trkFilePath);

% get leg positions with after alignment to fly midpoint and normalization 
%  to body length units
[srnLegX, srnLegY] = shiftRotateNormalizeLegPos(legX, legY, refPts);

% find indicies of bouts when the fly isn't moving, based on 2 midlegs not
%  moving
zeroVelInd = findFlyNotMovingMidlegs(srnLegX, srnLegY, notMoveParams);


%% load ephys data