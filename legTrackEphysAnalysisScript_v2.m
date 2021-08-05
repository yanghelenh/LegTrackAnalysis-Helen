% legTrackEphysAnalysisScript_v2.m
%
% Quick and dirty script for starting to look at leg tracking and ephys
%  responses together
%
% 8/5/21 - HHY

%% analysis constants
% indicies for specific body parts
% fit line to head and 3 pts defining midpt b/w legs
refPts.lineFitInd = 7:10; 
refPts.midPtInd = 9;
refPts.headPtInd = 7;
refPts.abdPtInd = 11;

% parameters for smoothing for determining leg velocities
smoParams.padLen = int32(50); % pad length in samples
smoParams.sigma = int32(10); % in samples

% parameters for determining whether fly is moving
notMoveParams.medFiltNumSamps = 10;
notMoveParams.zeroXVelThreshMed = 0.004;
notMoveParams.zeroXVelThresh = 0.01;
notMoveParams.movePosXVelThresh = 0.01;
notMoveParams.moveNegXVelThresh = -0.02;
notMoveParams.zeroYVelThreshMed = 0.008;
notMoveParams.zeroYVelThresh = 0.02;
notMoveParams.movePosYVelThresh = 0.01;
notMoveParams.moveNegYVelThresh = -0.02;
notMoveParams.minBoutLen = 10; % in samples
notMoveParams.r2LegInd = 2;
notMoveParams.l2LegInd = 5;

% parameters for finding velocity zero-crossings
zeroXingParams.legInd = 1:6;
zeroXingParams.legNames = {'R1', 'R2', 'R3', 'L1', 'L2', 'L3'}; 
zeroXingParams.minStepDur = 2;

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

%% load in pData
load(pDataFilePath, 'ephysData', 'ephysSpikes', 'fictrac', ...
    'fictracProc', 'fictracParams', 'leg');

%% pre-process leg tracking data

% load .trk file; get leg positions
[legX, legY] = loadTrkFile(trkFilePath);

% get leg positions after alignment to fly midpoint and normalization 
%  to body length units
[srnLegX, srnLegY] = shiftRotateNormalizeLegPos(legX, legY, refPts);

% get leg velocities, with Gaussian process smoothing on position
legXVel = findLegVel(srnLegX, smoParams);
legYVel = findLegVel(srnLegY, smoParams);

% find indicies of bouts when the fly isn't moving, based on 2 midlegs not
%  moving in direction parallel to long axis of fly
% zeroVelInd = findFlyNotMovingMidlegs(legXVel, notMoveParams);

% find indicies of bouts when fly isn't moving, based on 2 midlegs not
%  moving in x and y direction (updated from just x direction motion)
zeroVelInd = findFlyNotMovingMidlegsXY(legXVel, legYVel, notMoveParams);

% find zero crossings in velocity, to find min and max leg positions, for
%  determining steps
zeroXing = findVelZeroXing(legXVel, zeroVelInd, zeroXingParams);

% get start and end indicies of each half step
stepInd = determineStepStartEndInd(zeroXing, zeroXingParams.legInd);

% call swing and stance for all legs, using movement direction in direction
%  parallel to long axis of fly (x)
% [legStance, legSwing, legSwingStanceNotMove] = ...
%     callSwingStanceFrames(fictrac, leg.frameTimes, legXVel, ...
%     zeroVelInd, zeroXingParams.legInd);

% call swing stance for all legs, using movement in both x and y directions
[legStance, legSwing, legSwingStanceNotMove] = ...
    callSwingStanceFrames_v2(fictrac, leg.frameTimes, legXVel, ...
    legYVel, zeroVelInd, zeroXingParams.legInd);

% phase of leg position, in x direction
legXPhase = determineLegPhase(srnLegX, zeroXingParams.legInd);

% phase of leg position, in y direction
legYPhase = determineLegPhase(srnLegY, zeroXingParams.legInd);

% leg direction during step
stepDirs = determineStepDirection(srnLegX, srnLegY, stepInd, ...
    zeroXingParams.legInd);

% length of each step
stepLengths = determineStepLength(srnLegX, srnLegY, stepInd, ...
    zeroXingParams.legInd);