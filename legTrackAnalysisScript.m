% legTracksAnalysisScript.m
%
% Quick and dirty script for starting to look at leg tracking data
%
% various cells are things tried
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
smoParams.padLen = 50; % pad length in samples
smoParams.sigma = 10; % in samples

% parameters for determining whether fly is moving
% notMoveParams.medFiltNumSamps = 10;
% notMoveParams.zeroXVelThreshMed = 0.004;
% notMoveParams.zeroXVelThresh = 0.01;
% notMoveParams.movePosXVelThresh = 0.01;
% notMoveParams.moveNegXVelThresh = -0.02;
% notMoveParams.zeroYVelThreshMed = 0.008;
% notMoveParams.zeroYVelThresh = 0.02;
% notMoveParams.movePosYVelThresh = 0.01;
% notMoveParams.moveNegYVelThresh = -0.02;
% notMoveParams.minBoutLen = 10; % in samples
% notMoveParams.r2LegInd = 2;
% notMoveParams.l2LegInd = 5;


% parameters as of 8/11/21
notMoveParams.medFiltNumSamps = 10;
notMoveParams.zeroXVelThreshMed = 0.0004; %0.004;
notMoveParams.zeroXVelThresh = 0.02;
notMoveParams.movePosXVelThresh = 0.01;
notMoveParams.moveNegXVelThresh = -0.02;
notMoveParams.zeroYVelThreshMed = 0.004; %0.008;
notMoveParams.zeroYVelThresh = 0.04;
notMoveParams.movePosYVelThresh = 0.01;
notMoveParams.moveNegYVelThresh = -0.02;
notMoveParams.minBoutLen = 10; % in samples
notMoveParams.stepNegXVelThresh = - 0.035;
notMoveParams.maxTimeFromStep = 100; % in samples
% merge any not-moving bouts less than this many samples apart
notMoveParams.adjBoutSep = 50;  
notMoveParams.r2LegInd = 2;
notMoveParams.l2LegInd = 5;

% parameters for finding velocity zero-crossings
zeroXingParams.legInd = 1:6;
zeroXingParams.legNames = {'R1', 'R2', 'R3', 'L1', 'L2', 'L3'}; 
zeroXingParams.minStepDur = 2;

%% 200821_fly01_cell01_trial01 - mac
% full path to .trk leg tracking data
trkFilePath = ['D:\WilsonLab\Dropbox (HMS)\APTlegTracking\'...
    'Tracking_v1_mdn_labeled2798_ephys\' ...
    '200821_fly01_cell01_trial01_legVid_crop_caLegTrack_v1_mdn.trk'];
pDataFilePath = ['D:\WilsonLab\Dropbox (HMS)\EphysAnalysis-Helen\pData\'...
    '200821_fly01_cell01_trial01.mat'];

%% 200826_fly01_cell01_trial01 - mac
% full path to .trk leg tracking data
trkFilePath = ['D:\WilsonLab\Dropbox (HMS)\APTlegTracking\'...
    'Tracking_v1_mdn_labeled2798_ephys\' ...
    '200826_fly01_cell01_trial01_legVid_crop_caLegTrack_v1_mdn.trk'];
pDataFilePath = ['D:\WilsonLab\Dropbox (HMS)\EphysAnalysis-Helen\pData\'...
    '200826_fly01_cell01_trial01.mat'];

%% 200821_fly01_cell01_trial01 - windows
% full path to .trk leg tracking data
trkFilePath = ['D:\WilsonLab\Dropbox (HMS)\APTlegTracking\'...
    'Tracking_v1_mdn_labeled2798_ephys\' ...
    '200821_fly01_cell01_trial01_legVid_crop_caLegTrack_v1_mdn.trk'];
pDataFilePath = ['D:\WilsonLab\Dropbox (HMS)\EphysAnalysis-Helen\pData\'...
    '200821_fly01_cell01_trial01.mat'];

%% 200826_fly01_cell01_trial01 - windows
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

% some parameters about .trk file
numFrames = size(legX, 1); % number of frames in this trial
numPoints = size(legX, 2); % number of tracked points

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




%% THIS AND BELOW NOT USED %
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

%% PCA on leg X/Y position to extract 1 value describing position
% DON't USE - sometimes X and Y go in opposite directions for maximizing
%  variance, makes difference smaller

% loop over each leg
for i = 1%:length(zeroXingParams.legInd)
    
    % for conversion to leg indicies
    ind = zeroXingParams.legInd(i);
    
    % PCA function expects each row to be an observation and each column to
    %  be a variable (i.e. each row is a time point, 2 columns for x and y)
    rawMatPCA = [srnLegX(:,ind) srnLegY(:,ind)];
    
    % perform PCA
    [coeff, score, latent] = pca(rawMatPCA);
    
end

%% find reversals of leg position by smoothing aggressively first
% smooth leg position aggressively using sliding window
% find zero crossings in smoothed leg position
% given those zero crossings, take window around that point in original
%  data and find max/min

% CURRENT STATUS: bug in code that computes zero crossings. Fix this!

% srnLegXSmo = smoLegPos(srnLegX, smoParams);
% figure; plot(srnLegX(:,1)); hold on; plot(srnLegXSmo(:,1));

% some parameters
movAvgWinLen = 20; % length of window, in frames, for moving average

% get moving average for legX position
legXMovAvg = movmean(srnLegX,movAvgWinLen);

% compute velocity on movAvg
legXMovAvgVel = findLegVel(legXMovAvg, smoParams);

% get zero crossings on movAvg legX position
movAvgZeroXing = findVelZeroXing(legXMovAvgVel, zeroVelInd, ...
    zeroXingParams);

% plot zero crossings for R1
neg2PosR1Val = legXMovAvgVel(movAvgZeroXing.neg2Pos.R1,1);
pos2NegR1Val = legXMovAvgVel(movAvgZeroXing.pos2Neg.R1,1);


figure; plot(legXMovAvgVel(:,1));
hold on;
plot(movAvgZeroXing.neg2Pos.R1 ,neg2PosR1Val, 'x', 'LineStyle', 'none');
plot(movAvgZeroXing.pos2Neg.R1, pos2NegR1Val, 'o', 'LineStyle', 'none');



figure; plot(srnLegX(:,1)); hold on; plot(movAvg(:,1));

r1MovAvg = movmean(srnLegX(:,1),10);
figure; plot(srnLegX(:,1)); hold on; plot(r1MovAvg);

%% find reversals of leg position, attempt 4; CURRENT METHOD
% smooth leg position aggressively using sliding window
% using overlapping windows, find index of max, min
% if index of max/min isn't the first or last index of the window, flag the
%  window
% find the max/min within flagged windows in raw leg position data

% some parameters
movAvgWinLen = 30; % length of window, in frames, for moving average
maxminWinLen = 30; % length of window, in frames, for finding max/min
adjThresh = 4; % threshold for adjacent indicies

% get moving average for legX position
legXMovAvg = movmean(srnLegX(:,zeroXingParams.legInd),movAvgWinLen);
% legXMovAvg = movmedian(srnLegX(:,zeroXingParams.legInd),movAvgWinLen);
% smoParams.padLen = int32(50); % pad length in samples
% smoParams.sigma = int32(movAvgWinLen); % in samples
% legXMovAvg = smoLegPos(srnLegX, smoParams);

% generate 3D matrix from smoothed leg position, turn each position into
%  window; will be maxminWinLen shorter than original

% preallocate
legXWins = zeros(size(legXMovAvg,1)-maxminWinLen, size(legXMovAvg,2), ...
    maxminWinLen);

% populate legXWins
for i = 1:size(legXWins,1)
    for j = 1:size(legXWins, 2)
        legXWins(i,j,:) = legXMovAvg(i:i+maxminWinLen-1, j);
    end
end

% for each window, find index of max and min
[~,legXWinsMaxInds] = max(legXWins,[],3);
[~,legXWinsMinInds] = min(legXWins,[],3);


% find indicies where max/min is not 1 or window length
[legXWinsMaxRow,legXWinsMaxCol] = ind2sub(size(legXWinsMaxInds), ...
    find(legXWinsMaxInds ~= 1 & legXWinsMaxInds ~= maxminWinLen));

[legXWinsMinRow, legXWinsMinCol] = ind2sub(size(legXWinsMinInds), ...
    find(legXWinsMinInds ~= 1 & legXWinsMinInds ~= maxminWinLen));

% remove not-moving bouts from consideration
% maxNotMovLog = ismember(legXWinsMaxRow, zeroVelInd);
% 
% legXWinsMaxRowMov = legXWinsMaxRow(~maxNotMovLog);
% legXWinsMaxColMov = legXWinsMaxCol(~maxNotMovLog);
% 
% minNotMovLog = ismember(legXWinsMinRow, zeroVelInd);
% 
% legXWinsMinRowMov = legXWinsMinRow(~minNotMovLog);
% legXWinsMinColMov = legXWinsMinCol(~minNotMovLog);


% was removing not moving bouts here, but causes errors in calling some
%  min/max points
maxNotMovLog = ismember(legXWinsMaxRow, zeroVelInd);

legXWinsMaxRowMov = legXWinsMaxRow;
legXWinsMaxColMov = legXWinsMaxCol;

minNotMovLog = ismember(legXWinsMinRow, zeroVelInd);

legXWinsMinRowMov = legXWinsMinRow;
legXWinsMinColMov = legXWinsMinCol;

% merge adjacent windows, return start and end indicies of windows
[maxWinsRowsStartInd,maxWinsStartIndInd,maxWinsRowsEndInd, ...
    maxWinsEndIndInd] = mergeAdjVecVals(...
    legXWinsMaxRowMov,adjThresh);
maxWinsColsEndInd = legXWinsMaxColMov(maxWinsEndIndInd);

[minWinsRowsStartInd,minWinsStartIndInd,minWinsRowsEndInd, ...
    minWinsEndIndInd] = mergeAdjVecVals(...
    legXWinsMinRowMov, adjThresh);
minWinsColsEndInd = legXWinsMinColMov(minWinsEndIndInd);

% convert start and end indicies of windows to window midpoints
% maxWinMidIndInd = floor((maxWinsEndIndInd + maxWinsStartIndInd)./2);
% maxWinsRowsMidInd = legXWinsMaxRowMov(maxWinMidIndInd);
% maxWinsColsMidInd = legXWinsMaxColMov(maxWinMidIndInd);
% 
% minWinMidIndInd = floor((minWinsEndIndInd + minWinsStartIndInd)./2);
% minWinsRowsMidInd = legXWinsMinRowMov(minWinMidIndInd);
% minWinsColsMidInd = legXWinsMinColMov(minWinMidIndInd);


% for each window, find the index in the Gaussian-smoothed data where the
%  max/min actually occurred

% Gaussian process smooth leg position
% srnsLegX = smoLegPos(srnLegX, smoParams);
% srnsLegY = smoLegPos(srnLegY, smoParams);

% actually, don't Gaussian smooth, just use raw leg position
srnsLegX = srnLegX;
srnsLegY = srnLegY;

% generate matrix 

% preallocate
maxInds = zeros(length(maxWinsRowsEndInd), 1);
minInds = zeros(length(minWinsRowsEndInd), 1);

% find max within each window
for i = 1:length(maxWinsRowsEndInd)
    % get start index of window
    startInd = maxWinsRowsEndInd(i) - floor(maxminWinLen/2);
    % if start index is before start of vector, set it as vector start
    if (startInd < 1)
        startInd = 1;
    end
    % get end index of window
    endInd = maxWinsRowsEndInd(i) + floor(maxminWinLen/2);
    % if end index is after the end of the vecotr, set it as vector end
    if (endInd > size(srnsLegX,1))
        endInd = size(srnsLegX,1);
    end
    
    % get window
    window = srnsLegX(startInd:endInd,maxWinsColsEndInd(i));
    
    % get max within window
    [~, maxInd] = max(window);
    
    % save max (convert to index within srnsLegX)
    maxInds(i) = maxInd + startInd - 1;
    
end

% find min within each window
for i = 1:length(minWinsRowsEndInd)
    % get start index of window
    startInd = minWinsRowsEndInd(i) - floor(maxminWinLen/2);
    % if start index is before start of vector, set it as vector start
    if (startInd < 1)
        startInd = 1;
    end
    % get end index of window
    endInd = minWinsRowsEndInd(i) + floor(maxminWinLen/2);
    % if end index is after the end of the vecotr, set it as vector end
    if (endInd > size(srnsLegX,1))
        endInd = size(srnsLegX,1);
    end
    
    % get window
    window = srnsLegX(startInd:endInd,minWinsColsEndInd(i));
    
    % get max within window
    [~, minInd] = min(window);
    
    % save min (convert to index within srnsLegX)
    minInds(i) = minInd + startInd - 1;
end

% variables of interest: maxInds, maxWinsColsEndInd, minInds,
%  minWinsColsEndInd

% % preallocate
% maxInds = zeros(length(maxWinsRowsMidInd), 1);
% minInds = zeros(length(minWinsRowsMidInd), 1);
% 
% % find max within each window
% for i = 1:length(maxWinsRowsMidInd)
%     % get start index of window
%     startInd = maxWinsRowsMidInd(i) - floor(maxminWinLen/2);
%     % if start index is before start of vector, set it as vector start
%     if (startInd < 1)
%         startInd = 1;
%     end
%     % get end index of window
%     endInd = maxWinsRowsMidInd(i) + floor(maxminWinLen/2);
%     % if end index is after the end of the vecotr, set it as vector end
%     if (endInd > size(srnsLegX,1))
%         endInd = size(srnsLegX,1);
%     end
%     
%     % get window
%     window = srnsLegX(startInd:endInd,maxWinsColsMidInd(i));
%     
%     % get max within window
%     [~, maxInd] = max(window);
%     
%     % save max (convert to index within srnsLegX)
%     maxInds(i) = maxInd + startInd - 1;
%     
% end
% 
% % find min within each window
% for i = 1:length(minWinsRowsMidInd)
%     % get start index of window
%     startInd = minWinsRowsMidInd(i) - floor(maxminWinLen/2);
%     % if start index is before start of vector, set it as vector start
%     if (startInd < 1)
%         startInd = 1;
%     end
%     % get end index of window
%     endInd = minWinsRowsMidInd(i) + floor(maxminWinLen/2);
%     % if end index is after the end of the vecotr, set it as vector end
%     if (endInd > size(srnsLegX,1))
%         endInd = size(srnsLegX,1);
%     end
%     
%     % get window
%     window = srnsLegX(startInd:endInd,minWinsColsMidInd(i));
%     
%     % get max within window
%     [~, minInd] = min(window);
%     
%     % save min (convert to index within srnsLegX)
%     minInds(i) = minInd + startInd - 1;
% end
% 
% % variables of interest: maxInds, maxWinsColsMidInd, minInds,
% %  minWinsColsMidInd


% remove indicies within not-moving bouts

maxNotMovLog = ismember(maxInds, zeroVelInd);

maxIndsMov = maxInds(~maxNotMovLog);
maxColsMov = maxWinsColsEndInd(~maxNotMovLog);

minNotMovLog = ismember(minInds, zeroVelInd);

minIndsMov = minInds(~minNotMovLog);
minColsMov = minWinsColsEndInd(~minNotMovLog);

%% plots for identifying leg reversals

% plot
% maxIndsR1 = maxInds(maxWinsColsEndInd == 1);
% minIndsR1 = minInds(minWinsColsEndInd == 1);
% maxValsR1 = srnLegX(maxInds(maxWinsColsEndInd == 1),1);
% minValsR1 = srnLegX(minInds(minWinsColsEndInd == 1),1);
% 
% figure;
% plot(srnLegX(:,1));
% hold on;
% plot(maxIndsR1,maxValsR1, 'x', 'LineStyle', 'none');
% plot(minIndsR1,minValsR1, 'o', 'LineStyle', 'none');

% % plot
% maxIndsR2 = maxInds(maxWinsColsEndInd == 2);
% minIndsR2 = minInds(minWinsColsEndInd == 2);
% maxValsR2 = srnLegX(maxInds(maxWinsColsEndInd == 2),2);
% minValsR2 = srnLegX(minInds(minWinsColsEndInd == 2),2);
% 
% figure;
% plot(srnLegX(:,2));
% hold on;
% plot(maxIndsR2,maxValsR2, 'x', 'LineStyle', 'none');
% plot(minIndsR2,minValsR2, 'o', 'LineStyle', 'none');
% 
% % plot
% maxIndsR3 = maxInds(maxWinsColsEndInd == 3);
% minIndsR3 = minInds(minWinsColsEndInd == 3);
% maxValsR3 = srnLegX(maxInds(maxWinsColsEndInd == 3),3);
% minValsR3 = srnLegX(minInds(minWinsColsEndInd == 3),3);
% 
% figure;
% plot(srnLegX(:,3));
% hold on;
% plot(maxIndsR3,maxValsR3, 'x', 'LineStyle', 'none');
% plot(minIndsR3,minValsR3, 'o', 'LineStyle', 'none');

% plot
maxIndsR1 = maxIndsMov(maxColsMov == 1);
minIndsR1 = minIndsMov(minColsMov == 1);
maxValsR1 = srnLegX(maxIndsMov(maxColsMov == 1),1);
minValsR1 = srnLegX(minIndsMov(minColsMov == 1),1);

figure;
plot(srnLegX(:,1));
hold on;
plot(maxIndsR1,maxValsR1, 'x', 'LineStyle', 'none');
plot(minIndsR1,minValsR1, 'o', 'LineStyle', 'none');

plot(zeroVelInd, srnLegX(zeroVelInd,1), '.', 'LineStyle', 'none');

maxIndsR2 = maxIndsMov(maxColsMov == 2);
minIndsR2 = minIndsMov(minColsMov == 2);
maxValsR2 = srnLegX(maxIndsMov(maxColsMov == 2),2);
minValsR2 = srnLegX(minIndsMov(minColsMov == 2),2);

figure;
plot(srnLegX(:,2));
hold on;
plot(maxIndsR2,maxValsR2, 'x', 'LineStyle', 'none');
plot(minIndsR2,minValsR2, 'o', 'LineStyle', 'none');

plot(zeroVelInd, srnLegX(zeroVelInd,2), '.', 'LineStyle', 'none');

% plot
% maxIndsR1 = maxInds(maxWinsColsMidInd == 1);
% minIndsR1 = minInds(minWinsColsMidInd == 1);
% maxValsR1 = srnLegX(maxInds(maxWinsColsMidInd == 1),1);
% minValsR1 = srnLegX(minInds(minWinsColsMidInd == 1),1);
% 
% figure;
% plot(srnLegX(:,1));
% hold on;
% plot(maxIndsR1,maxValsR1, 'x', 'LineStyle', 'none');
% plot(minIndsR1,minValsR1, 'o', 'LineStyle', 'none');

%% Refine calls of leg reversals using 1st derivative of leg X position

% threshold of what 1st derivative (vel) should exceed in the positive
%  direction for a position maximum
legRevParams.maxPosVelThresh = 0.0001; 
legRevParams.maxNegVelThresh = -0.001; % in negative direction
legRevParams.minPosVelThresh = 0.0001; % in pos dir, for leg position min
legRevParams.minNegVelThresh = -0.001; % in neg dir, for leg position min
% number of frames before and after max/min to check for velocity thresh
legRevParams.numNegVelFrames = 8;
legRevParams.numPosVelFrames = 12;

% indicies (into maxIndsMov) of maxima to remove; initialize
maxRmvInd = [];

% loop through all maxima
for i = 1:length(maxIndsMov)
    % window before maximum
    bfWinStart = maxIndsMov(i) - legRevParams.numPosVelFrames;
    bfWinEnd = maxIndsMov(i) - 1; % index right before 
    % window after maximum
    afWinStart = maxIndsMov(i) + 1; % index right after
    afWinEnd = maxIndsMov(i) + legRevParams.numNegVelFrames;
    
    % if maximum occurs at first or last frame of the trial, only use 1
    %  window to check if maximum is valid
    if (maxIndsMov(i) == 1) % max is first frame of trial
        % only check that after window is valid
        afAvgVel = mean(legXVel(afWinStart:afWinEnd,maxColsMov(i)));
        % if invalid, save index
        if (afAvgVel > legRevParams.maxNegVelThresh)
            maxRmvInd = [maxRmvInd; i];
        end
        % continue to next iteration of loop
        continue;
    elseif (maxIndsMov(i) == numFrames) % max is last frame of trial
        % only check that before window is valid
        bfAvgVel = mean(legXVel(bfWinStart:bfWinEnd, maxColsMov(i)));
        % if invalid, save index
        if (bfAvgVel < legRevParams.maxPosVelThresh)
            maxRmvInd = [maxRmvInd; i];
        end
        % continue to next iteration of loop
        continue;
    end
    
    % check that windows are valid (contained b/w 1 and number of frames)
    %  if max is w/in numVelFrames of start or end of trial, shorten window
    %  accordingly (but can't be length zero, b/c taken care of with check
    %  on whether max is at start or end of trial
    
    % set before window to start with trial start
    if (bfWinStart < 1) 
        bfWinStart = 1;
    end
    % set after window to end with trial end
    if (afWinEnd > numFrames)
        afWinEnd = numFrames;
    end
    
    % get average velocities within window
    bfAvgVel = mean(legXVel(bfWinStart:bfWinEnd, maxColsMov(i)));
    afAvgVel = mean(legXVel(afWinStart:afWinEnd,maxColsMov(i)));
    
    % check that 1st derivative meets threshold for both windows, if not,
    %  add index to those to remove
    if (afAvgVel > legRevParams.maxNegVelThresh) || ...
            (bfAvgVel < legRevParams.maxPosVelThresh)
        maxRmvInd = [maxRmvInd; i];
    end
end

% remove all the indicies flagged for maxima
maxIndsMovCorr = maxIndsMov;
maxIndsMovCorr(maxRmvInd) = [];
maxColsMovCorr = maxColsMov;
maxColsMovCorr(maxRmvInd) = [];



% for minima

% initialize array of indicies to remove
minRmvInd = [];

% loop through all maxima
for i = 1:length(minIndsMov)
    % window before minimum
    bfWinStart = minIndsMov(i) - legRevParams.numNegVelFrames;
    bfWinEnd = minIndsMov(i) - 1; % index right before 
    % window after minimum
    afWinStart = minIndsMov(i) + 1; % index right after
    afWinEnd = minIndsMov(i) + legRevParams.numPosVelFrames;
    
    % if minimum occurs at first or last frame of the trial, only use 1
    %  window to check if minimum is valid
    if (minIndsMov(i) == 1) % min is first frame of trial
        % only check that after window is valid
        afAvgVel = mean(legXVel(afWinStart:afWinEnd, minColsMov(i)));
        % if invalid, save index
        if (afAvgVel < legRevParams.minPosVelThresh)
            minRmvInd = [minRmvInd; i];
        end
        % continue to next iteration of loop
        continue;
    elseif (minIndsMov(i) == numFrames) % min is last frame of trial
        % only check that before window is valid
        bfAvgVel = mean(legXVel(bfWinStart:bfWinEnd, minColsMov(i)));
        % if invalid, save index
        if (bfAvgVel > legRevParams.minNegVelThresh)
            minRmvInd = [minRmvInd; i];
        end
        % continue to next iteration of loop
        continue;
    end
    
    % check that windows are valid (contained b/w 1 and number of frames)
    %  if max is w/in numVelFrames of start or end of trial, shorten window
    %  accordingly (but can't be length zero, b/c taken care of with check
    %  on whether max is at start or end of trial
    
    % set before window to start with trial start
    if (bfWinStart < 1) 
        bfWinStart = 1;
    end
    % set after window to end with trial end
    if (afWinEnd > numFrames)
        afWinEnd = numFrames;
    end
    
    % get average velocities within window
    bfAvgVel = mean(legXVel(bfWinStart:bfWinEnd, minColsMov(i)));
    afAvgVel = mean(legXVel(afWinStart:afWinEnd,minColsMov(i)));
    
    % check that 1st derivative meets threshold for both windows, if not,
    %  add index to those to remove
    if (afAvgVel < legRevParams.minPosVelThresh) || ...
            (bfAvgVel > legRevParams.minNegVelThresh)
        minRmvInd = [minRmvInd; i];
    end
end

% remove all the indicies flagged for minima
minIndsMovCorr = minIndsMov;
minIndsMovCorr(minRmvInd) = [];
minColsMovCorr = minColsMov;
minColsMovCorr(minRmvInd) = [];

%% Refine calls of leg reversals using 1st derivative of leg X position v2
% Instead of taking mean of velocity within window, find mean position,
%  compute velocity as change in distance b/w this point and the max/min
%  over time. Set threshold on this velocity
% This should be less noisy...
% Works, but doesn't seem to be an improvement over v1 (above) 8/31/21
% DON'T USE (for now)

% threshold of what 1st derivative (vel) should exceed in the positive
%  direction for a position maximum
legRevParams.maxPosVelThresh = 0.00015; 
legRevParams.maxNegVelThresh = -0.0015; % in negative direction
legRevParams.minPosVelThresh = 0.00015; % in pos dir, for leg position min
legRevParams.minNegVelThresh = -0.0015; % in neg dir, for leg position min
% number of frames before and after max/min to check for velocity thresh
legRevParams.numNegVelFrames = 8;
legRevParams.numPosVelFrames = 12;

% indicies (into maxIndsMov) of maxima to remove; initialize
maxRmvInd = [];

% loop through all maxima
for i = 1:length(maxIndsMov)
    % window before maximum
    bfWinStart = maxIndsMov(i) - legRevParams.numPosVelFrames;
    bfWinEnd = maxIndsMov(i) - 1; % index right before 
    % window after maximum
    afWinStart = maxIndsMov(i) + 1; % index right after
    afWinEnd = maxIndsMov(i) + legRevParams.numNegVelFrames;
    % this leg postion and index
    thisLegInd = zeroXingParams.legInd(maxColsMov(i));
    thisLegPos = srnLegX(maxIndsMov(i), thisLegInd);
    
    % number of before and after frames
    numBfFrames = legRevParams.numPosVelFrames;
    numAfFrames = legRevParams.numNegVelFrames;
    
    % if maximum occurs at first or last frame of the trial, only use 1
    %  window to check if maximum is valid
    if (maxIndsMov(i) == 1) % max is first frame of trial
        % only check that after window is valid
        afAvgPos = mean(srnLegX(afWinStart:afWinEnd, thisLegInd));
        afAvgVel = (afAvgPos - thisLegPos) / numAfFrames;
        % if invalid, save index
        if (afAvgVel > legRevParams.maxNegVelThresh)
            maxRmvInd = [maxRmvInd; i];
        end
        % continue to next iteration of loop
        continue;
    elseif (maxIndsMov(i) == numFrames) % max is last frame of trial
        % only check that before window is valid
        bfAvgPos = mean(srnLegX(bfWinStart:bfWinEnd, thisLegInd));
        bfAvgVel = (thisLegPos - bfAvgPos) / numBfFrames;
        % if invalid, save index
        if (bfAvgVel < legRevParams.maxPosVelThresh)
            maxRmvInd = [maxRmvInd; i];
        end
        % continue to next iteration of loop
        continue;
    end
    
    % check that windows are valid (contained b/w 1 and number of frames)
    %  if max is w/in numVelFrames of start or end of trial, shorten window
    %  accordingly (but can't be length zero, b/c taken care of with check
    %  on whether max is at start or end of trial
    
    % set before window to start with trial start
    if (bfWinStart < 1) 
        bfWinStart = 1;
        % adjust count on number of frames for appropriate velocity calc
        numBfFrames = bfWinEnd - bfWinStart + 1;
    end
    % set after window to end with trial end
    if (afWinEnd > numFrames)
        afWinEnd = numFrames;
        % adjust count on number of frames for appropriate velocity calc
        numAfFrames = afWinEnd - afWinStart + 1;
    end
    
    % get average velocities within window
    bfAvgPos = mean(srnLegX(bfWinStart:bfWinEnd, thisLegInd));
    bfAvgVel = (thisLegPos - bfAvgPos) / numBfFrames;
    afAvgPos = mean(srnLegX(afWinStart:afWinEnd, thisLegInd));
    afAvgVel = (afAvgPos - thisLegPos) / numAfFrames;

    % check that 1st derivative meets threshold for both windows, if not,
    %  add index to those to remove
    if (afAvgVel > legRevParams.maxNegVelThresh) || ...
            (bfAvgVel < legRevParams.maxPosVelThresh)
        maxRmvInd = [maxRmvInd; i];
    end
end

% remove all the indicies flagged for maxima
maxIndsMovCorr = maxIndsMov;
maxIndsMovCorr(maxRmvInd) = [];
maxColsMovCorr = maxColsMov;
maxColsMovCorr(maxRmvInd) = [];



% for minima

% initialize array of indicies to remove
minRmvInd = [];

% loop through all maxima
for i = 1:length(minIndsMov)
    % window before minimum
    bfWinStart = minIndsMov(i) - legRevParams.numNegVelFrames;
    bfWinEnd = minIndsMov(i) - 1; % index right before 
    % window after minimum
    afWinStart = minIndsMov(i) + 1; % index right after
    afWinEnd = minIndsMov(i) + legRevParams.numPosVelFrames;
    
    % this leg postion and index
    thisLegInd = zeroXingParams.legInd(minColsMov(i));
    thisLegPos = srnLegX(minIndsMov(i), thisLegInd);
    
    % number of before and after frames
    numBfFrames = legRevParams.numNegVelFrames;
    numAfFrames = legRevParams.numPosVelFrames;
    
    % if minimum occurs at first or last frame of the trial, only use 1
    %  window to check if minimum is valid
    if (minIndsMov(i) == 1) % min is first frame of trial
        % only check that after window is valid
        afAvgPos = mean(srnLegX(afWinStart:afWinEnd, thisLegInd));
        afAvgVel = (afAvgPos - thisLegPos) / numAfFrames;
        % if invalid, save index
        if (afAvgVel < legRevParams.minPosVelThresh)
            minRmvInd = [minRmvInd; i];
        end
        % continue to next iteration of loop
        continue;
    elseif (minIndsMov(i) == numFrames) % min is last frame of trial
        % only check that before window is valid
        bfAvgPos = mean(srnLegX(bfWinStart:bfWinEnd, thisLegInd));
        bfAvgVel = (thisLegPos - bfAvgPos) / numBfFrames;
        % if invalid, save index
        if (bfAvgVel > legRevParams.minNegVelThresh)
            minRmvInd = [minRmvInd; i];
        end
        % continue to next iteration of loop
        continue;
    end
    
    % check that windows are valid (contained b/w 1 and number of frames)
    %  if max is w/in numVelFrames of start or end of trial, shorten window
    %  accordingly (but can't be length zero, b/c taken care of with check
    %  on whether max is at start or end of trial
    
    % set before window to start with trial start
    if (bfWinStart < 1) 
        bfWinStart = 1;
        % adjust count on number of frames for appropriate velocity calc
        numBfFrames = bfWinEnd - bfWinStart + 1;
    end
    % set after window to end with trial end
    if (afWinEnd > numFrames)
        afWinEnd = numFrames;
        % adjust count on number of frames for appropriate velocity calc
        numAfFrames = afWinEnd - afWinStart + 1;
    end
    
    % get average velocities within window
    bfAvgPos = mean(srnLegX(bfWinStart:bfWinEnd, thisLegInd));
    bfAvgVel = (thisLegPos - bfAvgPos) / numBfFrames;
    afAvgPos = mean(srnLegX(afWinStart:afWinEnd, thisLegInd));
    afAvgVel = (afAvgPos - thisLegPos) / numAfFrames;
    
    % check that 1st derivative meets threshold for both windows, if not,
    %  add index to those to remove
    if (afAvgVel < legRevParams.minPosVelThresh) || ...
            (bfAvgVel > legRevParams.minNegVelThresh)
        minRmvInd = [minRmvInd; i];
    end
end

% remove all the indicies flagged for minima
minIndsMovCorr = minIndsMov;
minIndsMovCorr(minRmvInd) = [];
minColsMovCorr = minColsMov;
minColsMovCorr(minRmvInd) = [];

%% plots for min/max, after refinement with 1st derivative

for i = 1:length(zeroXingParams.legInd)
    thisLeg = zeroXingParams.legInd(i);
    
    theseMaxInds = maxIndsMovCorr(maxColsMovCorr == thisLeg);
    theseMinInds = minIndsMovCorr(minColsMovCorr == thisLeg);
    theseMaxVals = srnLegX(maxIndsMovCorr(maxColsMovCorr == thisLeg),thisLeg);
    theseMinVals = srnLegX(minIndsMovCorr(minColsMovCorr == thisLeg),thisLeg);
    
    figure;
    plot(srnLegX(:,thisLeg));
    hold on;
    plot(theseMaxInds, theseMaxVals, 'x', 'LineStyle', 'none');
    plot(theseMinInds, theseMinVals, 'o', 'LineStyle', 'none');
    
    % add title; leg and whether X or Y position
    title([zeroXingParams.legNames{i} ' leg X position']);
end

% maxIndsR1 = maxIndsMovCorr(maxColsMovCorr == 1);
% minIndsR1 = minIndsMovCorr(minColsMovCorr == 1);
% maxValsR1 = srnLegX(maxIndsMovCorr(maxColsMovCorr == 1),1);
% minValsR1 = srnLegX(minIndsMovCorr(minColsMovCorr == 1),1);
% 
% figure;
% plot(srnLegX(:,1));
% hold on;
% plot(maxIndsR1,maxValsR1, 'x', 'LineStyle', 'none');
% plot(minIndsR1,minValsR1, 'o', 'LineStyle', 'none');
% 
% plot(zeroVelInd, srnLegX(zeroVelInd,1), '.', 'LineStyle', 'none');
% 
% maxIndsR2 = maxIndsMovCorr(maxColsMovCorr == 2);
% minIndsR2 = minIndsMovCorr(minColsMovCorr == 2);
% maxValsR2 = srnLegX(maxIndsMovCorr(maxColsMovCorr == 2),2);
% minValsR2 = srnLegX(minIndsMovCorr(minColsMovCorr == 2),2);
% 
% figure;
% plot(srnLegX(:,2));
% hold on;
% plot(maxIndsR2,maxValsR2, 'x', 'LineStyle', 'none');
% plot(minIndsR2,minValsR2, 'o', 'LineStyle', 'none');
% 
% plot(zeroVelInd, srnLegX(zeroVelInd,2), '.', 'LineStyle', 'none');

%% convert maxInds, minInds to steps
% each step as triplet [startInd, midInd, endInd], steps defined as
%  starting with leg at front-most position (minInds)
% filter to keep steps within the same moving bout
% CURRENTLY WORKING as of 9/2/21

% get start and end indicies of not moving bouts
% if not moving bouts separated by <= this value, merge them
zeroVelAdjThresh = 50; % number of samples
[notMoveStartInd, ~, notMoveEndInd, ~] = mergeAdjVecVals(...
    zeroVelInd, zeroVelAdjThresh);
% aggressive bout merging here, b/c would otherwise break up steps

% no merging of bouts
% [notMoveStartInd, notMoveEndInd, notMoveBoutDur] = findBouts(zeroVelInd);

% plot shading for not-moving bouts
notMovingX = [notMoveStartInd'; notMoveStartInd'; notMoveEndInd'; notMoveEndInd'];
y0 = ones(size(notMoveStartInd)) * -0.3;
y1 = ones(size(notMoveStartInd)) * 0.6;
notMovingY = [y0'; y1'; y1'; y0'];

figure;
plot(srnLegX(:,2));
hold on;
plot(zeroVelInd, srnLegX(zeroVelInd,2), '.', 'LineStyle', 'none');
patch(notMovingX, notMovingY, 'black', 'FaceAlpha', 0.3');



% initialize step arrays
stepInds = zeros(length(minIndsMovCorr),3); 
% corresponds to Cols arrays when IDing local min/maxes
stepsWhichLeg = zeros(length(minIndsMovCorr),1);

% initialize step counter
counter = 1;

% debugging: for each minIndsMovCorr point that doesn't become a step, keep
%  track of why
disqualPts.notMove = []; % start point in not moving bout
disqualPts.minFMin = []; % start point is followed by another min
disqualPts.noMidpt = []; % no midpoint for this start point in trial
% midpoint not in same moving bout as start point
disqualPts.midptDiffMoveBout = []; 
disqualPts.maxFMax = []; % midpoint is followed by another max
disqualPts.noEndpt = []; % no endpoint for this start point in trial
% endpoint is not in same moving bout as start point
disqualPts.endptDiffMoveBout = [];




% loop through all minInds elements
for i = 1:length(minIndsMovCorr)
    
    % frame index of this potential start point (minimum)
    thisPotStart = minIndsMovCorr(i);
    
    % get leg for potential start point
    thisLeg = minColsMovCorr(i);
    
    % consider minIndsMovCorr and maxIndsMovCorr only for this leg
    thisLegMinInds = minIndsMovCorr(minColsMovCorr == thisLeg);
    thisLegMaxInds = maxIndsMovCorr(maxColsMovCorr == thisLeg);
    
    % index of this potential start point into
    %  thisLegMinInds/thisLegMaxInds
    thisPotStartInd = find(thisLegMinInds == thisPotStart);
    
    
    % check if this minInd is within not moving bout
    
    % index into notMoveStartInd of not move bout start point immediately
    % preceding this minInd
    thisNotMoveStartInd = find(thisPotStart > notMoveStartInd, 1, ...
        'last');
    % if there could be a not moving bout that this belongs to (empty if
    %  this minInd falls earlier than start of first not moving bout; i.e.
    %  if fly is walking at very beginning of trial)
    if ~(isempty(thisNotMoveStartInd))
        % get frame index value for start of not moving bout
        thisNotMoveStart = notMoveStartInd(thisNotMoveStartInd);
        % get frame index value for end of not moving bout (starts and ends
        %  are paired)
        thisNotMoveEnd = notMoveEndInd(thisNotMoveStartInd);
        
        % check if this potential start point falls within this not moving 
        %  bout (is its value less than thisNotMoveEnd?)
        % if yes, continue this loop without executing rest
        if (thisPotStart < thisNotMoveEnd) 
            disqualPts.notMove = [disqualPts.notMove; i];
            continue;
        end
    end
    
    % convert start and end indicies of not moving bouts to start and end
    %  indicies for moving bout that this step start point belongs to
    
    % index of notMoveEndInd for last time this minInd exceeds, to get
    %  start index of moving bout
    thisMoveStart = find(thisPotStart > notMoveEndInd, 1, 'last');
    
    % if this minInd never greater than notMoveEndInd, fly must be moving
    %  at the beginning
    if ~isempty(thisMoveStart)
        thisMoveStartInd = notMoveEndInd(thisMoveStart) + 1;
    else
        thisMoveStartInd = 1; % start index is start of trial
    end
    
    % index of notMoveStartInd for first time notMoveStartInd is greater
    % than thisMoveStartInd (so when this moving bout ends)
    thisMoveEnd = find(notMoveStartInd > thisMoveStartInd, 1, 'first');
    
    % if end of move period can't be found b/c no subsequent start index,
    %  fly must be moving at the end of the trial
    if ~isempty(thisMoveEnd)
        thisMoveEndInd = notMoveStartInd(thisMoveEnd) - 1;
    else
        thisMoveEndInd = size(srnLegX,1); % end index is end of trial
    end
    
    
    % using this potential start point, find indicies
    %  that define step mid point (from maxInd) and end point (from minInd)
    stepStartPt = thisPotStart;
    
    % check whether this step start point (as a min), is followed by
    %  another min or by a max
    % should be a max, but if not, don't use this start point
    % check if there is a next min point
    % this start index is the last min point for this leg
    if (thisPotStartInd == length(thisLegMinInds))
        % if there's no next min point, there must not be end point for
        % this step
        disqualPts.noEndpt = [disqualPts.noEndpt; i];
        continue;
    else % there are additional points
        nextMinPt = thisLegMinInds(thisPotStartInd + 1); 
        nextMaxPt = thisLegMaxInds(find(thisLegMaxInds > stepStartPt, 1,...
            'first'));
        % min follows start point, not max; don't use this start point
        if (nextMinPt < nextMaxPt)
            disqualPts.minFMin = [disqualPts.minFMin; i];
            continue;
        end
    end
    
    % find mid point - maxInds value that immediately follows step start pt
    % index into thisLegMaxInds of this point
    midPtInd = find(thisLegMaxInds > stepStartPt, 1, 'first');
    % if there is no such point, this full step is not found in the trial,
    %  continue this loop without executing rest
    if isempty(midPtInd)
        disqualPts.noMidpt = [disqualPts.noMidpt; i];
        continue;
    end
    
    % convert this index into maxInd into frame index of step mid point
    stepMidPt = thisLegMaxInds(midPtInd);
    
    % check that this midpoint is within the same moving bout as the step
    %  start point; if not, continue loop without executing the rest
    if (stepMidPt < thisMoveStartInd) || (stepMidPt > thisMoveEndInd)
        disqualPts.midptDiffMoveBout = [disqualPts.midptDiffMoveBout; i];
        continue;
    end
    
    % check whether this step mid point (as a max), is followed by
    %  another max or by a min
    % should be a min, but if not, don't use this start point
    % check whether there is a subsequent max point (if not, don't use this
    %  check)
    if (midPtInd ~= length(thisLegMaxInds))
        nextMaxPt = thisLegMaxInds(midPtInd + 1);
        nextMinPt = thisLegMinInds(thisPotStartInd + 1); % same as above
        if (nextMinPt > nextMaxPt)
            disqualPts.maxFMax = [disqualPts.maxFMax; i];
            continue;
        end
    end
    
    
    % find end point - minInd value that immediately follows step mid pt
    % index into thisLegMinInds of this point
    endPtInd = find(thisLegMinInds > stepMidPt, 1, 'first');
    % if there is no such point, this full step is not found in the trial,
    %  continue this loop without executing rest
    if isempty(endPtInd)
        disqualPts.noEndpt = [disqualPts.noEndPt; i];
        continue;
    end

    % convert this index into minInd into frame index of step end point
    stepEndPt = thisLegMinInds(endPtInd);
    
    % check that this endpoint is within the same moving bout as the step
    %  start point; if not, continue loop without exectuing the rest
    if (stepEndPt < thisMoveStartInd) || (stepEndPt > thisMoveEndInd)
        disqualPts.endptDiffMoveBout = [disqualPts.endptDiffMoveBout; i];
        continue;
    end
    
    % if this part of the code is reached, the step start, mid, and end
    %  points are all valid; add this step to the step arrays
    % step indicies
    stepInds(counter,:) = [stepStartPt, stepMidPt, stepEndPt];
    % which leg
    stepsWhichLeg(counter) = thisLeg;
    
    % update counter
    counter = counter + 1;
end

% stepInds and stepsWhichLeg are shorter than they were intialized, because
% not all minInds points became steps; clear the excess zeros from the end
stepInds = stepInds(1:(counter-1), :);
stepsWhichLeg = stepsWhichLeg(1:(counter-1));

%% plot steps, overlayed on leg positions in X or Y  
% one plot for each leg

% leg X positions
for i = 1:length(zeroXingParams.legInd)
    thisLeg = zeroXingParams.legInd(i);
    
    % stepInds for this leg
    stepIndsThisLeg = stepInds(stepsWhichLeg==thisLeg,:);
    
    % subtract 1 from end index
    stepIndsThisLeg(:,3) = stepIndsThisLeg(:,3) - 1;
    
    figure;
    
    plot(srnLegX(:,thisLeg), 'k'); % plot leg position
    hold on;
    
    % get leg position value for each step
    % preallocate
    stepLegPos = zeros(size(stepIndsThisLeg));

    for j = 1:size(stepIndsThisLeg,1)
        stepLegPos(j,:) = srnLegX(stepIndsThisLeg(j,:),thisLeg)';
    end
    
    % plot steps
    plot(stepIndsThisLeg',stepLegPos', 'x', 'MarkerSize', 10);
    
    % add title
    title([zeroXingParams.legNames{i} ' leg X position']);
    
end

% leg Y positions
for i = 1:length(zeroXingParams.legInd)
    thisLeg = zeroXingParams.legInd(i);
    
    % stepInds for this leg
    stepIndsThisLeg = stepInds(stepsWhichLeg==thisLeg,:);
    
    % subtract 1 from end index
    stepIndsThisLeg(:,3) = stepIndsThisLeg(:,3) - 1;
    
    figure;
    
    plot(srnLegY(:,thisLeg), 'k'); % plot leg position
    hold on;
    
    % get leg position value for each step
    % preallocate
    stepLegPos = zeros(size(stepIndsThisLeg));

    for j = 1:size(stepIndsThisLeg,1)
        stepLegPos(j,:) = srnLegY(stepIndsThisLeg(j,:),thisLeg)';
    end
    
    % plot steps
    plot(stepIndsThisLeg',stepLegPos', 'x', 'MarkerSize', 10);
    
    % add title; leg and whether X or Y position
    title([zeroXingParams.legNames{i} ' leg Y position']);
    
end

%% Improve ID of walking/not-walking: parameter tweaking
% TESTING ONLY - currently (8/12/21) not in use

% notMoveParams.medFiltNumSamps = 10;
% notMoveParams.zeroVelThreshMed = 0.0001; %0.004;
% notMoveParams.zeroVelThresh = 0.01;
% notMoveParams.movePosVelThresh = 0.01;
% notMoveParams.moveNegVelThresh = -0.02;
% notMoveParams.minBoutLen = 10; % in samples
% notMoveParams.r2LegInd = 2;
% notMoveParams.l2LegInd = 5;

notMoveParams.medFiltNumSamps = 10;
notMoveParams.zeroXVelThreshMed = 0.0004; %0.004;
notMoveParams.zeroXVelThresh = 0.02;
notMoveParams.movePosXVelThresh = 0.01;
notMoveParams.moveNegXVelThresh = -0.02;
notMoveParams.zeroYVelThreshMed = 0.004; %0.008;
notMoveParams.zeroYVelThresh = 0.04;
notMoveParams.movePosYVelThresh = 0.01;
notMoveParams.moveNegYVelThresh = -0.02;
notMoveParams.minBoutLen = 10; % in samples
notMoveParams.stepNegXVelThresh = - 0.035;
notMoveParams.maxTimeFromStep = 100; % in samples
% merge any not-moving bouts less than this many samples apart
notMoveParams.adjBoutSep = 20;  
notMoveParams.r2LegInd = 2;
notMoveParams.l2LegInd = 5;



% find indicies of bouts when the fly isn't moving, based on 2 midlegs not
%  moving in direction parallel to long axis of fly
% zeroVelInd = findFlyNotMovingMidlegs(legXVel, notMoveParams);


% find indicies of bouts when fly isn't moving, based on 2 midlegs not
%  moving in x and y direction (updated from just x direction motion)
zeroVelInd = findFlyNotMovingMidlegsXY(legXVel, legYVel, notMoveParams);


% plot leg position and not moving bouts
figure;
plot(srnLegX(:,2));
hold on;
plot(zeroVelInd, srnLegX(zeroVelInd,2), '.', 'LineStyle', 'none');

figure;
plot(legXVel(:,2));
hold on;
plot(zeroVelInd, legXVel(zeroVelInd,2), '.', 'LineStyle', 'none');

%% Improve ID of walking/not-walking: PCA, ICA, k-means
% DON'T USE!! %
% doesn't separate walking from not-walking cleanly enough

% generate matrix for PCA: X, Y for 2 midlegs
velMatPCA = [legXVel(:,notMoveParams.r2LegInd), ...
    legXVel(:,notMoveParams.l2LegInd), ...
    legYVel(:,notMoveParams.r2LegInd),...
    legYVel(:,notMoveParams.l2LegInd)];

velMatPCA = [legXVel(:,notMoveParams.r2LegInd), ...
    legYVel(:,notMoveParams.r2LegInd)];

% moving average for velocities
legVelXMovAvg = movmean(legXVel(:,zeroXingParams.legInd),10);
legVelYMovAvg = movmean(legYVel(:,zeroXingParams.legInd),10);

velMatPCA = [legVelXMovAvg(:,notMoveParams.r2LegInd), ...
    legVelXMovAvg(:,notMoveParams.l2LegInd), ...
    legVelYMovAvg(:,notMoveParams.r2LegInd),...
    legVelYMovAvg(:,notMoveParams.l2LegInd)];

[coeff, score, ~,~,explained,~] = pca(velMatPCA);

% colormap
cm = repmat([0 0 1], size(legXVel,1),1);
for i = 1:length(zeroVelInd)
    cm(zeroVelInd(i),:) = [1 0 0];
end


figure;
s = scatter(score(:,1), score(:,2), 36, cm, 'filled', 'MarkerFaceAlpha', ...
    0.03, 'MarkerEdgeAlpha', 0.03);

figure; 
s = scatter(score(:,1), score(:,2), 'filled', 'MarkerFaceAlpha', ...
    0.03, 'MarkerEdgeAlpha', 0.03);


% k-means clustering
legVelXMovAvg = movmean(legXVel(:,zeroXingParams.legInd),10);
legVelYMovAvg = movmean(legYVel(:,zeroXingParams.legInd),10);

velMatPCA = [legVelXMovAvg(:,notMoveParams.r2LegInd), ...
    legVelXMovAvg(:,notMoveParams.l2LegInd)];


clusters = kmeans(velMatPCA,3);

% colormap
cm = zeros(size(legXVel,1),3);
for i = 1:length(clusters)
    switch clusters(i)
        case 1
            cm(i,:) = [0 0 1];
        case 2
            cm(i,:) = [1 0 0];
        case 3
            cm(i,:) = [0 1 0];
    end
end

figure;
s = scatter(legXVel(:,notMoveParams.r2LegInd), ...
    legXVel(:, notMoveParams.l2LegInd), ...
    36, cm, 'filled', 'MarkerFaceAlpha', 0.03, 'MarkerEdgeAlpha', 0.03);

% ICA
velMatPCA = [legVelXMovAvg(:,notMoveParams.r2LegInd), ...
    legVelXMovAvg(:,notMoveParams.l2LegInd), ...
    legVelYMovAvg(:,notMoveParams.r2LegInd),...
    legVelYMovAvg(:,notMoveParams.l2LegInd)];

icaMdl = rica(velMatPCA, 2);
icaFeatures = transform(icaMdl, velMatPCA);

cm = repmat([0 0 1], size(legXVel,1),1);
for i = 1:length(zeroVelInd)
    cm(zeroVelInd(i),:) = [1 0 0];
end

figure;
scatter(icaFeatures(:,1),icaFeatures(:,2), 36, cm, 'filled', 'MarkerFaceAlpha', ...
    0.03, 'MarkerEdgeAlpha', 0.03);

%% Compute parameters for steps %%

%% Step length (in xy plane)
%% Step direction (in xy plane)
%% Step duration (for each half step as well as whole step)
%% Velocity of each half step (in xy plane, x and y components)

%% Swing/stance call for each half step using comparison to FicTrac

% interpolate fictrac forward velocity to legVidFrametimes
interpFtFwd = interp1(fictracProc.t, fictracProc.fwdVel, leg.frameTimes);




%% Swing/stance call for each half step 

