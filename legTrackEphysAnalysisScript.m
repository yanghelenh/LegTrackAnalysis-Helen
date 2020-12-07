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

% parameters for smoothing for determining leg velocities
smoParams.padLen = int32(50); % pad length in samples
smoParams.sigma = int32(10); % in samples

% parameters for determining whether fly is moving
notMoveParams.medFiltNumSamps = 10;
notMoveParams.zeroVelThreshMed = 0.004;
notMoveParams.zeroVelThresh = 0.01;
notMoveParams.movePosVelThresh = 0.01;
notMoveParams.moveNegVelThresh = -0.02;
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
zeroVelInd = findFlyNotMovingMidlegs(legXVel, notMoveParams);

% find zero crossings in velocity, to find min and max leg positions, for
%  determining steps
zeroXing = findVelZeroXing(legXVel, zeroVelInd, zeroXingParams);

% get start and end indicies of each half step
stepInd = determineStepStartEndInd(zeroXing, zeroXingParams.legInd);

% call swing and stance for all legs
[legStance, legSwing, legSwingStanceNotMove] = ...
    callSwingStanceFrames(fictrac, leg.frameTimes, legXVel, ...
    zeroVelInd, zeroXingParams.legInd);

% phase of leg position, in x direction
legXPhase = determineLegPhase(srnLegX, zeroXingParams.legInd);

% leg direction during step
stepDirs = determineStepDirection(srnLegX, srnLegY, stepInd, ...
    zeroXingParams.legInd);

% length of each step
stepLengths = determineStepLength(srnLegX, srnLegY, stepInd, ...
    zeroXingParams.legInd);

%% 
winLen = 0.25; % duration of window, in sec
winShift = 0.025; % offset between windows


winEndTime = leg.frameTimes(end) - winLen;

fracStanceAll = [];
% for each leg, for each time window, find amount of time spent in stance
for i = 1:length(zeroXingParams.legInd)
    winStartTime = leg.frameTimes(1); % start time of window
    fracStanceLeg = [];
    winStartTimesAll = [];
    % walk through whole trial
    while (winStartTime < winEndTime)
        winStartInd = find(leg.frameTimes <= winStartTime, 1, 'last');
        winEndInd = find(leg.frameTimes <= (winStartTime + winLen),...
            1, 'last');
        
        swStBit = legSwingStanceNotMove(winStartInd:winEndInd, ...
            zeroXingParams.legInd(i));
        % if fly is moving for the whole window (no zeros in snippet)
        if (isempty(find(swStBit==0, 1)))
            % find fraction of window in stance: stance is 1, swing is -1
            % (mean + 1) /2 is fraction of time in stance
            fracStance = (mean(swStBit) + 1) / 2;
        
            fracStanceLeg = [fracStanceLeg; fracStance];
            winStartTimesAll = [winStartTimesAll; winStartTime];
        end
        
        % update winStartTime by winShift
        winStartTime = winStartTime + winShift;
    end
    % append fracStance for this leg into whole matrix
    fracStanceAll = [fracStanceAll, fracStanceLeg];
end

% compute ratio of fraction of time spent in stance for left and right back
% legs (ind 3 and 6); right:left
otherLegs = [1,4:6];

otherLegsMeanFracStance = mean(fracStanceAll(:,otherLegs),2);

% ratioR3L3FracStance = fracStanceAll(:,3) ./ fracStanceAll(:,6);
ratioR3L3FracStance = fracStanceAll(:,3) ./ otherLegsMeanFracStance;


% for these same time windows, find the mean yaw velocity (FicTrac)
% preallocate vector for mean yaw velocities
meanYawVel = zeros(size(ratioR3L3FracStance));
% loop through all time windows
for i = 1:length(winStartTimesAll)
    winStartInd = find(fictracProc.t <= winStartTimesAll(i), 1, 'last');
    winEndInd = find(fictracProc.t <= (winStartTimesAll(i) + winLen), ...
        1, 'last');
    
    yawVelBit = fictracProc.yawAngVel(winStartInd:winEndInd);
    meanYawVel(i) = mean(yawVelBit);
end


% for these same time windows, find the mean spike rate 
% shift ephys data by ephysOffset prior to behavior
ephysOffset = 0.1; 
% preallocate vector for mean spike rate
meanSpikeRate = zeros(size(ratioR3L3FracStance));
% loop through all time windows
for i = 1:length(winStartTimesAll)
    winStartInd = find(ephysSpikes.t <= ...
        (winStartTimesAll(i)-ephysOffset), 1, 'last');
    winEndInd = find(ephysSpikes.t <= ...
        (winStartTimesAll(i)-ephysOffset + winLen), 1, 'last');
    
    numSpikes = length(find((ephysSpikes.startInd >= winStartInd) & ...
        (ephysSpikes.startInd <= winEndInd)));
    % spike rate is number of spikes in window divided by length of window
    meanSpikeRate(i) = numSpikes / winLen;
end

% plot as scatterplot of ratio vs. yaw velocity, colored by spike rate
figure;
scatter(meanYawVel, ratioR3L3FracStance, 25, meanSpikeRate, 'filled');
% scatter(meanYawVel, fracStanceAll(:,3), 25, meanSpikeRate, 'filled');
xlabel('Yaw Velocity (deg/s)');
ylabel('Ratio of Stance duration R3:rest');
colorbar;


%% plots

f1XLims = [234 240];
figure;
subplot(5,1,1);
plot(leg.frameTimes, srnLegX(:,1:6));
% hold on; 
% plot(leg.frameTimes(zeroVelInd), srnLegX(zeroVelInd,2), '.');
xlim(f1XLims);
legend(zeroXingParams.legNames);

subplot(5,1,2);
plot(leg.frameTimes, srnLegY(:,1:6));
xlim(f1XLims);
legend(zeroXingParams.legNames);

subplot(5,1,3);
imagesc(legSwingStanceNotMove');
swingStanceLims = (f1XLims-leg.frameTimes(1)) .* ...
    (1/median(diff(leg.frameTimes)));
xlim(swingStanceLims);

subplot(5,1,4);
plot(ephysData.t, ephysData.scaledVoltage);
xlim(f1XLims);

% subplot(5,1,4);
% plot(ephysSpikes.t, ephysSpikes.spikeRate);
% xlim(f1XLims);

subplot(5,1,5);
plot(fictracProc.t, fictracProc.yawAngVel);
xlim(f1XLims);