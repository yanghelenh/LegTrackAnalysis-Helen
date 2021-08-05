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

%% looking at step length, spike rate, yaw vel

delay = 0.05; % time in sec for ephys to precede behavior

stepLim = [0 0.6];
spikeLim = [0 200];

colorRes = 255;
colorbarWhite = 0;
zScale = [-200 200];
numTicks = 11; % number of ticks for colorbar

% zScale
colorScale = redblue(colorRes);
minColorInd = 1;

% max difference from value set to white
maxAmp = max(abs(zScale - colorbarWhite));
% use maxAmp to determine scale factor between zScale and colorScale:
% maxAmp corresponds to 1/2 of colorScale; c per z
zcScFctr = ((colorRes - minColorInd)/2) / maxAmp;

% scale factors for colorbar
% exclude 1 from showing on colorbar
colorbarLims = [minColorInd colorRes]; 

% where ticks are, on zScale; will become tick labels
tickLabels = zScale(1):((zScale(2)-zScale(1))/(numTicks - 1)):zScale(2);
% where ticks are, in indicies
tickLocs = zcScFctr .* (tickLabels - colorbarWhite) + ...
    (floor(colorRes/2) + 1);
tickLocs = (tickLocs - minColorInd)/colorRes;
colorbarLims = [0 1];
        
stanceCorrCoeff = zeros(1,length(zeroXingParams.legNames));
% front to back (stance, mostly)
for i = 1:length(zeroXingParams.legNames)
%     allSpikeRates = zeros(size(stepLengths.front2Back.(zeroXingParams.legNames{i})));
%     allYawVel = zeros(size(stepLengths.front2Back.(zeroXingParams.legNames{i})));

    allSpikeRates = [];
    allYawVel = [];
    allStepLengths = [];
    thisStepLength = stepLengths.front2Back.(zeroXingParams.legNames{i});
    
    for j = 1:size(stepInd.front2Back.(zeroXingParams.legNames{i}),1)
        theseInd = stepInd.front2Back.(zeroXingParams.legNames{i});
        startTime = leg.frameTimes(theseInd(j,1));
        endTime = leg.frameTimes(theseInd(j,2));
        
        if (startTime < endTime)
        
            ephysStartInd = find(ephysSpikes.t > (startTime-delay), 1,'first');
            ephysEndInd = find(ephysSpikes.t < (endTime-delay), 1, 'last');

            spikes = find(ephysSpikes.startInd >= ephysStartInd & ...
                ephysSpikes.startInd <= ephysEndInd);

            numSpikes = numel(spikes);

            spikeRate = numSpikes / ...
                (ephysSpikes.t(ephysEndInd) - ephysSpikes.t(ephysStartInd));
            allSpikeRates= [allSpikeRates spikeRate];

            fictracStartInd = find(fictracProc.t > startTime, 1, 'first');
            fictracEndInd = find(fictracProc.t < endTime, 1, 'last');

            avgYawVel = mean(fictracProc.yawAngVel(fictracStartInd:fictracEndInd));
            allYawVel = [allYawVel avgYawVel];
            

            stepLength = thisStepLength(j);
            
            allStepLengths = [allStepLengths stepLength];
        end
    end
    
    yawVelInd = round(zcScFctr .* (allYawVel - colorbarWhite) + ...
        (floor(colorRes/2) + 1));
    
    yawVelInd(yawVelInd < minColorInd) = minColorInd;
    yawVelInd(yawVelInd > colorRes) = colorRes;
    yawVelCol = colorScale(yawVelInd,:);
    
    stanceCorrCoeff(i) = corr(allStepLengths', allSpikeRates');
            
    figure;
%     colormap(colorScale);
%     scatter(allStepLengths, allSpikeRates, 25, yawVelCol, 'filled');
%     colorbarHandle = colorbar('LimitsMode', 'manual', 'Limits',...
%         colorbarLims, 'TicksMode', 'manual', 'Ticks', tickLocs, ...
%         'TickLabelsMode', 'manual', 'TickLabels', tickLabels);
    scatter(allStepLengths, allSpikeRates, 25, 'filled');
    xlim(stepLim);
    ylim(spikeLim);
	xlabel('Step Length (body lengths)');
    ylabel('Spike Rate (spikes/sec)');
%     colorbarHandle.Label.String = 'Yaw Vel (deg/sec)';
    title(zeroXingParams.legNames{i});
%     set(gca, 'Color', [0.3 0.3 0.3]);
end


swingCorrCoeff = zeros(1,length(zeroXingParams.legNames));
% back to front (swing, mostly)
for i = 1:length(zeroXingParams.legNames)
%     allSpikeRates = zeros(size(stepLengths.front2Back.(zeroXingParams.legNames{i})));
%     allYawVel = zeros(size(stepLengths.front2Back.(zeroXingParams.legNames{i})));

    allSpikeRates = [];
    allYawVel = [];
    allStepLengths = [];
    thisStepLength = stepLengths.back2Front.(zeroXingParams.legNames{i});
    
    for j = 1:size(stepInd.back2Front.(zeroXingParams.legNames{i}),1)
        theseInd = stepInd.back2Front.(zeroXingParams.legNames{i});
        startTime = leg.frameTimes(theseInd(j,1));
        endTime = leg.frameTimes(theseInd(j,2));
        
        if (startTime < endTime)
        
            ephysStartInd = find(ephysSpikes.t > (startTime-delay), 1,'first');
            ephysEndInd = find(ephysSpikes.t < (endTime-delay), 1, 'last');

            spikes = find(ephysSpikes.startInd >= ephysStartInd & ...
                ephysSpikes.startInd <= ephysEndInd);

            numSpikes = numel(spikes);

            spikeRate = numSpikes / ...
                (ephysSpikes.t(ephysEndInd) - ephysSpikes.t(ephysStartInd));
            allSpikeRates= [allSpikeRates spikeRate];

            fictracStartInd = find(fictracProc.t > startTime, 1, 'first');
            fictracEndInd = find(fictracProc.t < endTime, 1, 'last');

            avgYawVel = mean(fictracProc.yawAngVel(fictracStartInd:fictracEndInd));
            allYawVel = [allYawVel avgYawVel];
            

            stepLength = thisStepLength(j);
            
            allStepLengths = [allStepLengths stepLength];
        end
    end
    
    yawVelInd = round(zcScFctr .* (allYawVel - colorbarWhite) + ...
        (floor(colorRes/2) + 1));
    
    yawVelInd(yawVelInd < minColorInd) = minColorInd;
    yawVelInd(yawVelInd > colorRes) = colorRes;
    yawVelCol = colorScale(yawVelInd,:);
    
    swingCorrCoeff(i) = corr(allStepLengths', allSpikeRates');
            
    figure;
%     colormap(colorScale);
%     scatter(allStepLengths, allSpikeRates, 25, yawVelCol, 'filled');
%     colorbarHandle = colorbar('LimitsMode', 'manual', 'Limits',...
%         colorbarLims, 'TicksMode', 'manual', 'Ticks', tickLocs, ...
%         'TickLabelsMode', 'manual', 'TickLabels', tickLabels);
    scatter(allStepLengths, allSpikeRates, 25, 'filled');
    xlim(stepLim);
    ylim(spikeLim);
	xlabel('Step Length (body lengths)');
    ylabel('Spike Rate (spikes/sec)');
%     colorbarHandle.Label.String = 'Yaw Vel (deg/sec)';
    title(zeroXingParams.legNames{i});
%     set(gca, 'Color', [0.3 0.3 0.3]);
end

%% looking at step length, spike rate, fwd vel


delay = 0.05; % time in sec for ephys to precede behavior

stepLim = [0 0.6];
spikeLim = [0 200];

colorRes = 255;
colorbarWhite = 0;
zScale = [0 10];
numTicks = 11; % number of ticks for colorbar

% zScale
colorScale = parula(colorRes);
minColorInd = 1; % shift corresponding to 1 vs 0 indexing

zcScFctr = (colorRes - minColorInd) / (zScale(2)-zScale(1));

% scale factors for colorbar
% exclude 1 from showing on colorbar
colorbarLims = [minColorInd colorRes]; 

% where ticks are, on zScale; will become tick labels
colorbarLims = [minColorInd colorRes]; 

% where ticks are, on zScale; will become tick labels
tickLabels = zScale(1):((zScale(2)-zScale(1))/(numTicks - 1)):zScale(2);
% where ticks are, in indicies
tickLocs = zcScFctr .* (tickLabels - zScale(1)) + minColorInd;
tickLocs = (tickLocs - minColorInd)/colorRes;
colorbarLims = [0 1];
        
% front to back (stance, mostly)
for i = 1:length(zeroXingParams.legNames)
%     allSpikeRates = zeros(size(stepLengths.front2Back.(zeroXingParams.legNames{i})));
%     allYawVel = zeros(size(stepLengths.front2Back.(zeroXingParams.legNames{i})));

    allSpikeRates = [];
    allFwdVel = [];
    allStepLengths = [];
    thisStepLength = stepLengths.front2Back.(zeroXingParams.legNames{i});
    
    for j = 1:size(stepInd.front2Back.(zeroXingParams.legNames{i}),1)
        theseInd = stepInd.front2Back.(zeroXingParams.legNames{i});
        startTime = leg.frameTimes(theseInd(j,1));
        endTime = leg.frameTimes(theseInd(j,2));
        
        if (startTime < endTime)
        
            ephysStartInd = find(ephysSpikes.t > startTime, 1,'first');
            ephysEndInd = find(ephysSpikes.t < endTime, 1, 'last');

            spikes = find(ephysSpikes.startInd >= ephysStartInd & ...
                ephysSpikes.startInd <= ephysEndInd);

            numSpikes = numel(spikes);

            spikeRate = numSpikes / ...
                (ephysSpikes.t(ephysEndInd) - ephysSpikes.t(ephysStartInd));
            allSpikeRates= [allSpikeRates spikeRate];

            fictracStartInd = find(fictracProc.t > startTime, 1, 'first');
            fictracEndInd = find(fictracProc.t < endTime, 1, 'last');

            avgFwdVel = mean(fictracProc.fwdVel(fictracStartInd:fictracEndInd));
            allFwdVel = [allFwdVel avgFwdVel];
            

            stepLength = thisStepLength(j);
            
            allStepLengths = [allStepLengths stepLength];
        end
    end
    
    fwdVelInd = round(zcScFctr .* (allFwdVel - zScale(1)) + ...
        minColorInd);
    
    fwdVelInd(fwdVelInd < minColorInd) = minColorInd;
    fwdVelInd(fwdVelInd > colorRes) = colorRes;
    fwdVelCol = colorScale(fwdVelInd,:);
            
    figure;
    colormap(colorScale);
    scatter(allStepLengths, allSpikeRates, 25, fwdVelCol, 'filled');
    colorbarHandle = colorbar('LimitsMode', 'manual', 'Limits',...
        colorbarLims, 'TicksMode', 'manual', 'Ticks', tickLocs, ...
        'TickLabelsMode', 'manual', 'TickLabels', tickLabels);
    xlim(stepLim);
    ylim(spikeLim);
	xlabel('Step Length (body lengths)');
    ylabel('Spike Rate (spikes/sec)');
    colorbarHandle.Label.String = 'Fwd Vel (mm/sec)';
    title(zeroXingParams.legNames{i});
end

        
% back to front (swing, mostly)
for i = 1:length(zeroXingParams.legNames)
%     allSpikeRates = zeros(size(stepLengths.front2Back.(zeroXingParams.legNames{i})));
%     allYawVel = zeros(size(stepLengths.front2Back.(zeroXingParams.legNames{i})));

    allSpikeRates = [];
    allFwdVel = [];
    allStepLengths = [];
    thisStepLength = stepLengths.back2Front.(zeroXingParams.legNames{i});
    
    for j = 1:size(stepInd.back2Front.(zeroXingParams.legNames{i}),1)
        theseInd = stepInd.back2Front.(zeroXingParams.legNames{i});
        startTime = leg.frameTimes(theseInd(j,1));
        endTime = leg.frameTimes(theseInd(j,2));
        
        if (startTime < endTime)
        
            ephysStartInd = find(ephysSpikes.t > startTime, 1,'first');
            ephysEndInd = find(ephysSpikes.t < endTime, 1, 'last');

            spikes = find(ephysSpikes.startInd >= ephysStartInd & ...
                ephysSpikes.startInd <= ephysEndInd);

            numSpikes = numel(spikes);

            spikeRate = numSpikes / ...
                (ephysSpikes.t(ephysEndInd) - ephysSpikes.t(ephysStartInd));
            allSpikeRates= [allSpikeRates spikeRate];

            fictracStartInd = find(fictracProc.t > startTime, 1, 'first');
            fictracEndInd = find(fictracProc.t < endTime, 1, 'last');

            avgFwdVel = mean(fictracProc.fwdVel(fictracStartInd:fictracEndInd));
            allFwdVel = [allFwdVel avgFwdVel];
            

            stepLength = thisStepLength(j);
            
            allStepLengths = [allStepLengths stepLength];
        end
    end
    
    fwdVelInd = round(zcScFctr .* (allFwdVel - zScale(1)) + ...
        minColorInd);
    
    fwdVelInd(fwdVelInd < minColorInd) = minColorInd;
    fwdVelInd(fwdVelInd > colorRes) = colorRes;
    fwdVelCol = colorScale(fwdVelInd,:);
            
    figure;
    colormap(colorScale);
    scatter(allStepLengths, allSpikeRates, 25, fwdVelCol, 'filled');
    colorbarHandle = colorbar('LimitsMode', 'manual', 'Limits',...
        colorbarLims, 'TicksMode', 'manual', 'Ticks', tickLocs, ...
        'TickLabelsMode', 'manual', 'TickLabels', tickLabels);
    xlim(stepLim);
    ylim(spikeLim);
	xlabel('Step Length (body lengths)');
    ylabel('Spike Rate (spikes/sec)');
    colorbarHandle.Label.String = 'Fwd Vel (mm/sec)';
    title(zeroXingParams.legNames{i});
end

%% looking at step direction, spike rate, yaw vel

delay = 0.05; % time in sec for ephys to precede behavior

stepLim = [-180 180];
spikeLim = [0 200];

colorRes = 255;
colorbarWhite = 0;
zScale = [-200 200];
numTicks = 11; % number of ticks for colorbar

% zScale
colorScale = redblue(colorRes);
minColorInd = 1;

% max difference from value set to white
maxAmp = max(abs(zScale - colorbarWhite));
% use maxAmp to determine scale factor between zScale and colorScale:
% maxAmp corresponds to 1/2 of colorScale; c per z
zcScFctr = ((colorRes - minColorInd)/2) / maxAmp;

% scale factors for colorbar
% exclude 1 from showing on colorbar
colorbarLims = [minColorInd colorRes]; 

% where ticks are, on zScale; will become tick labels
tickLabels = zScale(1):((zScale(2)-zScale(1))/(numTicks - 1)):zScale(2);
% where ticks are, in indicies
tickLocs = zcScFctr .* (tickLabels - colorbarWhite) + ...
    (floor(colorRes/2) + 1);
tickLocs = (tickLocs - minColorInd)/colorRes;
colorbarLims = [0 1];
        

stanceDirCorrCoeff = zeros(1, length(zeroXingParams.legNames));
% front to back (stance, mostly)
for i = 1:length(zeroXingParams.legNames)
%     allSpikeRates = zeros(size(stepLengths.front2Back.(zeroXingParams.legNames{i})));
%     allYawVel = zeros(size(stepLengths.front2Back.(zeroXingParams.legNames{i})));

    allSpikeRates = [];
    allYawVel = [];
    allStepDirs = [];
    thisStepDir = wrapTo180(stepDirs.front2Back.(zeroXingParams.legNames{i}) + 90);
    
    for j = 1:size(stepInd.front2Back.(zeroXingParams.legNames{i}),1)
        theseInd = stepInd.front2Back.(zeroXingParams.legNames{i});
        startTime = leg.frameTimes(theseInd(j,1));
        endTime = leg.frameTimes(theseInd(j,2));
        
        if (startTime < endTime)
        
            ephysStartInd = find(ephysSpikes.t > (startTime-delay), 1,'first');
            ephysEndInd = find(ephysSpikes.t < (endTime-delay), 1, 'last');

            spikes = find(ephysSpikes.startInd >= ephysStartInd & ...
                ephysSpikes.startInd <= ephysEndInd);

            numSpikes = numel(spikes);

            spikeRate = numSpikes / ...
                (ephysSpikes.t(ephysEndInd) - ephysSpikes.t(ephysStartInd));
            allSpikeRates= [allSpikeRates spikeRate];

            fictracStartInd = find(fictracProc.t > startTime, 1, 'first');
            fictracEndInd = find(fictracProc.t < endTime, 1, 'last');

            avgYawVel = mean(fictracProc.yawAngVel(fictracStartInd:fictracEndInd));
            allYawVel = [allYawVel avgYawVel];
            

            stepDir = thisStepDir(j);
            
            allStepDirs = [allStepDirs stepDir];
        end
    end
    
    yawVelInd = round(zcScFctr .* (allYawVel - colorbarWhite) + ...
        (floor(colorRes/2) + 1));
    
    yawVelInd(yawVelInd < minColorInd) = minColorInd;
    yawVelInd(yawVelInd > colorRes) = colorRes;
    yawVelCol = colorScale(yawVelInd,:);
    
    stanceDirCorrCoeff(i) = corr(allStepDirs',allSpikeRates');
            
    figure;
    colormap(colorScale);
    scatter(allStepDirs, allSpikeRates, 25, yawVelCol, 'filled');
    colorbarHandle = colorbar('LimitsMode', 'manual', 'Limits',...
        colorbarLims, 'TicksMode', 'manual', 'Ticks', tickLocs, ...
        'TickLabelsMode', 'manual', 'TickLabels', tickLabels);
%     scatter(allStepDirs, allSpikeRates, 25, 'filled');
    xlim(stepLim);
    ylim(spikeLim);
	xlabel('Step Dir (deg)');
    ylabel('Spike Rate (spikes/sec)');
    colorbarHandle.Label.String = 'Yaw Vel (deg/sec)';
    title(zeroXingParams.legNames{i});
    set(gca, 'Color', [0.3 0.3 0.3]);
end

swingDirCorrCoeff = zeros(1,length(zeroXingParams.legNames));
% back to front (swing, mostly)
for i = 1:length(zeroXingParams.legNames)
%     allSpikeRates = zeros(size(stepLengths.front2Back.(zeroXingParams.legNames{i})));
%     allYawVel = zeros(size(stepLengths.front2Back.(zeroXingParams.legNames{i})));

    allSpikeRates = [];
    allYawVel = [];
    allStepDirs = [];
    thisStepDir = wrapTo180(stepDirs.back2Front.(zeroXingParams.legNames{i}) + 90);
    
    for j = 1:size(stepInd.back2Front.(zeroXingParams.legNames{i}),1)
        theseInd = stepInd.back2Front.(zeroXingParams.legNames{i});
        startTime = leg.frameTimes(theseInd(j,1));
        endTime = leg.frameTimes(theseInd(j,2));
        
        if (startTime < endTime)
        
            ephysStartInd = find(ephysSpikes.t > (startTime-delay), 1,'first');
            ephysEndInd = find(ephysSpikes.t < (endTime-delay), 1, 'last');

            spikes = find(ephysSpikes.startInd >= ephysStartInd & ...
                ephysSpikes.startInd <= ephysEndInd);

            numSpikes = numel(spikes);

            spikeRate = numSpikes / ...
                (ephysSpikes.t(ephysEndInd) - ephysSpikes.t(ephysStartInd));
            allSpikeRates= [allSpikeRates spikeRate];

            fictracStartInd = find(fictracProc.t > startTime, 1, 'first');
            fictracEndInd = find(fictracProc.t < endTime, 1, 'last');

            avgYawVel = mean(fictracProc.yawAngVel(fictracStartInd:fictracEndInd));
            allYawVel = [allYawVel avgYawVel];
            

            stepDir = thisStepDir(j);
            
            allStepDirs = [allStepDirs stepDir];
        end
    end
    
    yawVelInd = round(zcScFctr .* (allYawVel - colorbarWhite) + ...
        (floor(colorRes/2) + 1));
    
    yawVelInd(yawVelInd < minColorInd) = minColorInd;
    yawVelInd(yawVelInd > colorRes) = colorRes;
    yawVelCol = colorScale(yawVelInd,:);
    
    swingDirCorrCoeff(i) = corr(allStepDirs',allSpikeRates');
            
    figure;
    colormap(colorScale);
    scatter(allStepDirs, allSpikeRates, 25, yawVelCol, 'filled');
    colorbarHandle = colorbar('LimitsMode', 'manual', 'Limits',...
        colorbarLims, 'TicksMode', 'manual', 'Ticks', tickLocs, ...
        'TickLabelsMode', 'manual', 'TickLabels', tickLabels);
%     scatter(allStepDirs, allSpikeRates, 25, 'filled');
    xlim(stepLim);
    ylim(spikeLim);
	xlabel('Step Length (body lengths)');
    ylabel('Spike Rate (spikes/sec)');
    colorbarHandle.Label.String = 'Yaw Vel (deg/sec)';
    title(zeroXingParams.legNames{i});
    set(gca, 'Color', [0.3 0.3 0.3]);
end

%% plots

f1XLims = [234 240];
figure;
ax1 = subplot(5,1,1);
plot(leg.frameTimes, srnLegX(:,1:6));
% hold on; 
% plot(leg.frameTimes(zeroVelInd), srnLegX(zeroVelInd,2), '.');
legend(zeroXingParams.legNames);
title('Normalized leg position, front-back axis');
ylabel('Body lengths');


ax2 = subplot(5,1,2);
plot(leg.frameTimes, srnLegY(:,1:6));
legend(zeroXingParams.legNames);
title('Normalized leg position, side axis');
ylabel('Body lengths');

% ax3 = subplot(6,1,3);
% % swingStanceLims = (f1XLims-leg.frameTimes(1)) .* ...
% %    (1/median(diff(leg.frameTimes)));
% swingStanceXlims = [leg.frameTimes(1) leg.frameTimes(end)];
% swingStanceYlims = [1 6];
% imagesc(swingStanceXlims, swingStanceYlims, legSwingStanceNotMove');
% title('Swing/stance calls (yellow - stance, blue - swing, teal - not moving');
% ylabel({'R1';'R2';'R3';'L1';'L2';'L3'});
% ylh = get(gca,'ylabel');
% gyl = get(ylh);                                                        
% ylp = get(ylh, 'Position');
% set(ylh, 'Rotation',0, 'Position',ylp, 'VerticalAlignment','middle', ...
%     'HorizontalAlignment','right');

% xlim(swingStanceLims);

ax4 = subplot(5,1,3);
plot(ephysData.t, ephysData.scaledVoltage);
title('Scaled voltage');
ylabel('mV');

% ax4 = subplot(6,1,4);
% plot(ephysSpikes.t, ephysSpikes.spikeRate);

ax5 = subplot(5,1,4);
plot(fictracProc.t, fictracProc.yawAngVel);
title('Yaw velocity');
ylabel('deg/s');

ax6 = subplot(5,1,5);
plot(fictracProc.t, fictracProc.fwdVel);
title('Forward velocity');
ylabel('mm/s');

linkaxes([ax1 ax2 ax4 ax5 ax6], 'x');
xlim(f1XLims);

%% 
figure;
plot(leg.frameTimes, legXVel(:,2));
hold on;
plot(leg.frameTimes(zeroVelInd), legXVel(zeroVelInd,2),'.');