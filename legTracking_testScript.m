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

%% load in data
% file paths
trialPath = ...
    ['/Users/hyang/Documents/LegTrackAnalysis-Helen/TestFiles/' ...
    '190507_fly01_fov02_trial01'];

legTrackPath = ...
    '/Users/hyang/Documents/LegTrackAnalysis-Helen/TestFiles';
legTrackFile = '190507_fly01_fov02_trial01_legVid_caLegTrack_v1_mdn.trk';

legVidPath = [trialPath filesep 'rawLegVid'];



curDir = pwd; % remember current directory

% load in leg tracking data from .trk file; need to load as .mat
cd(legTrackPath);
load(legTrackFile, '-mat', 'pTrk');

% load in leg video data, FicTrac data
cd(trialPath);
load('legVidDat.mat', 'legVidFrameTimes'); % leg vid frame times
load('pData.mat', 'fictrac'); % filtered & smoothed FicTrac data
load('fictracDat.mat'); % only preprocessed FicTrac data
ficTracTimes = t;

%% shift, rotate, normalize leg position data

% leg tracking
% get x and y coordinates separately
legX = squeeze(pTrk(:,1,:))';
legY = squeeze(pTrk(:,2,:))';

% get leg image size
legTrackInd = 76; % which frame (doesn't really matter)
legImgFullPath = sprintf('%s%slegvid-%d.tiff', legVidPath, ...
    filesep, legTrackInd-1);
legVidImg =  imread(legImgFullPath);
legVidMinInt = 50;
legVidMaxInt = 150;

% first dimension when read in with imread is vertical size, second is
%  horizontal
legImgSize = size(legVidImg);
legImgVertPx = legImgSize(1);
legImgHorizPx = legImgSize(2);

% fit line to 4 points defining head and 3 thorax-coxa joints; [7-10]
lineFitInd = 7:10;
% preallocate vector for saving coefficients of line fit
bodyLineCoeffs = zeros(size(legX,1), 2);

% get coefficients of line for each video frame; first is slope, second is
%  y-intercept
for i = 1:size(legX, 1)
    xPts = legX(i, lineFitInd);
    yPts = legY(i, lineFitInd);
    
    bodyLineCoeffs(i,:) = polyfit(xPts, yPts, 1);
end

% get projection of midleg thorax-coxa point onto body line (define this as
%  midpoint of fly)
% midpoints
midPtInd = 9;
[midXProj, midYProj] = projPt2Line(bodyLineCoeffs(:,1), ...
    bodyLineCoeffs(:,2), legX(:,midPtInd), legY(:,midPtInd));

% project head point and abdomen point onto body line, get distance between
%  them to get body length
headPtInd = 7;
abdPtInd = 11;
[headXProj, headYProj] = projPt2Line(bodyLineCoeffs(:,1), ...
    bodyLineCoeffs(:,2), legX(:,headPtInd), legY(:,headPtInd));
[abdXProj, abdYProj] = projPt2Line(bodyLineCoeffs(:,1), ...
    bodyLineCoeffs(:,2), legX(:,abdPtInd), legY(:,abdPtInd));

bodyLen = distBtw2Pts(headXProj, headYProj, abdXProj, abdYProj);

% shift all points so that they're relative to the midpoint, which is now
%  defined as (0,0)
% preallocate
legXShift = zeros(size(legX));
legYShift = zeros(size(legY));
% shift each point
for i = 1:size(legX,2)
    legXShift(:,i) = legX(:,i) - midXProj;
    legYShift(:,i) = legY(:,i) - midYProj;
end

% get body line coefficients, after shift
bodyLineShiftCoeffs = zeros(size(legXShift,1),2);

for i = 1:size(legXShift, 1)
    xPts = legXShift(i, lineFitInd);
    yPts = legYShift(i, lineFitInd);
    
    bodyLineShiftCoeffs(i,:) = polyfit(xPts, yPts, 1);
end

% for whole time series, define the body length and the body angle
% for each frame, rotate around the midpoint by the body angle so x and y
%  axes align with forward and lateral axes of fly
% normalize units to body length (from pixels)

% body length as median of body lengths throughout trial
flyBodyLength = median(bodyLen);

% tangent of body angle = median slope of body line, 
%  post shift (y-int now ~0)
flyBodyAng = atand(median(bodyLineShiftCoeffs(:,1)));

% rotate all points about midpt, now (0,0), by inverse of fly body angle
% preallocate
legXShiftRot = zeros(size(legXShift));
legYShiftRot = zeros(size(legYShift));
% loop through each marked point and rotate
for i = 1:size(legXShift,2)
    [legXShiftRot(:,i), legYShiftRot(:,i)] = ...
        rotatePts2D(legXShift(:,i), legYShift(:,i), -1*flyBodyAng);
end

% normalize to body length
legXShiftRotNorm = legXShiftRot / flyBodyLength;
legYShiftRotNorm = legYShiftRot / flyBodyLength;

% multiply legYShiftRotNorm by -1 so that fly's right side has positive
%  coordinates and left side has negative coordinates
legYShiftRotNorm = -1 * legYShiftRotNorm;

%% determine when fly isn't moving
% when both right and left midleg velocity close to zero (by threshold)

% smooth leg position traces using Gaussian process smoothing
padLen = int32(50);
sigma = int32(10);

% preallocate
legXShiftRotNormSmo = zeros(size(legXShiftRotNorm));

for i = 1:size(legXShiftRotNormSmo,2)
    smoPosPy = py.proc_utils.safe_interp_conv_smooth_ball_data(...
        legXShiftRotNorm(:,i)', padLen, sigma);
    % convert from python to matlab data format
    smoPos = cell2mat(cell(smoPosPy.tolist()));
    legXShiftRotNormSmo(:,i) = smoPos';
end

% get leg x velocities by taking gradient
legXVel = zeros(size(legXShiftRotNormSmo));
for i = 1:size(legXVel,2)
    legXVel(:,i) = gradient(legXShiftRotNormSmo(:,i));
end

% median filter leg velocities, mid legs
medFiltR2 = medfilt1(legXVel(:,2),10);
medFiltL2 = medfilt1(legXVel(:,5),10);

% threshold for moving in median filtered velocity, non-median filtered
% velocity
zeroVelThreshMed = 0.004;
zeroVelThresh = 0.01;
% not moving as speed less than threshold
r2ZeroVelInd = find((abs(medFiltR2) < zeroVelThreshMed) & ...
    (abs(legXVel(:,2)) < zeroVelThresh));
% r2ZeroVelIndDiff = diff(r2ZeroVelInd);
% r2ZeroVelIndDiff = [1; r2ZeroVelIndDiff];

% also has to be within 2 frames of subsequent not-moving frame (so zero
%  crossings aren't counted, only more sustained stops)
% r2ZeroVelIndLong = r2ZeroVelInd(r2ZeroVelIndDiff<=2);

% find inverse of R2 zeroVelInd -> moving index
allInd = (1:length(legXVel(:,2)))';
r2MoveInd = setdiff(allInd, r2ZeroVelInd); 

% find moving bouts
[r2MoveBoutStarts, r2MoveBoutEnds, r2MoveBoutDur] = findBouts(r2MoveInd);

% find bouts where movement doesn't exceed threshold
movePosVelThresh = 0.01; % threshold for positive velocity
moveNegVelThresh = -0.02; % threshold for negative velocity

% loop through all moving bouts
for i = 1:length(r2MoveBoutStarts)
    boutInd = (r2MoveBoutStarts(i):r2MoveBoutEnds(i))';
    boutVels = legXVel(boutInd,2);
    % if movement bout doesn't exceed threshold, append to zeroVelInd
    if ~(sum(boutVels >= movePosVelThresh) || ...
            sum(boutVels <= moveNegVelThresh))
        r2ZeroVelInd = [r2ZeroVelInd; boutInd];
    end
end

% resort r2ZeroVelInd in ascending order
r2ZeroVelInd = sort(r2ZeroVelInd);

% same for left mid leg
l2ZeroVelInd = find((abs(medFiltL2) < zeroVelThreshMed) & ...
    (abs(legXVel(:,5)) < zeroVelThresh));
% l2ZeroVelIndDiff = diff(l2ZeroVelInd);
% l2ZeroVelIndDiff = [1; l2ZeroVelIndDiff];

% l2ZeroVelIndLong = l2ZeroVelInd(l2ZeroVelIndDiff<=2);

l2MoveInd = setdiff(allInd, l2ZeroVelInd);

[l2MoveBoutStarts, l2MoveBoutEnds, l2MoveBoutDur] = findBouts(l2MoveInd);

for i = 1:length(l2MoveBoutStarts)
    boutInd = (l2MoveBoutStarts(i):l2MoveBoutEnds(i))';
    boutVels = legXVel(boutInd,5);
    % if movement bout doesn't exceed threshold, append to zeroVelInd
    if ~(sum(boutVels >= movePosVelThresh) || ...
            sum(boutVels <= moveNegVelThresh))
        l2ZeroVelInd = [l2ZeroVelInd; boutInd];
    end
end

l2ZeroVelInd = sort(l2ZeroVelInd);

% fly not moving when both right and left mid legs not moving (fly can
%  groom with front or back legs and be stationary, but they don't groom
%  with midlegs)
% zeroVelInd = intersect(r2ZeroVelIndLong, l2ZeroVelIndLong);
zeroVelInd = intersect(r2ZeroVelInd, l2ZeroVelInd);

% each not moving bout must be of minimum length; if not, remove
minZeroVelBoutLen = 5; % minimum length of 5 samples

% find start, end, and duration of zero velocity bouts
[zeroVelBoutStarts, zeroVelBoutEnds, zeroVelBoutDur] = ...
    findBouts(zeroVelInd);

% find bouts that are less than minimum long
zeroVelShortBouts = find(zeroVelBoutDur < minZeroVelBoutLen);

% remove these bouts from zeroVelInd - set all values within zeroVelInd
%  that fall into these short bouts to NaN
% loop through all short bouts, check if zeroVelInd falls within start and
%  end
for i = 1:length(zeroVelShortBouts)
    startInd = zeroVelBoutStarts(zeroVelShortBouts(i));
    endInd = zeroVelBoutEnds(zeroVelShortBouts(i));
    
    zeroVelInd((zeroVelInd >= startInd)&(zeroVelInd <= endInd)) = NaN;
end
% remove NaNs from zeroVelInd
zeroVelInd(isnan(zeroVelInd)) = [];


% % fill in gaps
% allInd = (1:length(legXVel(:,2)))';
% moveInd = setdiff(allInd, zeroVelInd); 
% 
% % find moving bouts
% [moveBoutStarts, moveBoutEnds, moveBoutDur] = findBouts(moveInd);
% 
% % loop through all moving bouts
% for i = 1:length(moveBoutStarts)
%     boutInd = (moveBoutStarts(i):moveBoutEnds(i))';
%     boutR2Vels = legXVel(boutInd,2);
%     boutL2Vels = legXVel(boutInd,5);
%     % if movement bout doesn't exceed threshold, append to zeroVelInd
%     if ~(sum(boutR2Vels >= movePosVelThresh) || ...
%             sum(boutR2Vels <= moveNegVelThresh)) && ...
%             ~(sum(boutL2Vels >= movePosVelThresh) || ...
%             sum(boutL2Vels <= moveNegVelThresh))
%         zeroVelInd = [zeroVelInd; boutInd];
%     end
% end
% 
% % resort r2ZeroVelInd in ascending order
% zeroVelInd = sort(zeroVelInd);

%% find zero crossings in leg X velocity

legInd = [1:6]; % indicies for legs
% leg names corresponding to leg indicies
legNames = {'R1', 'R2', 'R3', 'L1', 'L2', 'L3'}; 
allInd = (1:length(legXVel(:,2)))';
minStepDur = 2; % minimum step duration, in samples
moveInd = setdiff(allInd, zeroVelInd);

% find moving bouts
[moveBoutStarts, moveBoutEnds, moveBoutDur] = findBouts(moveInd);

% preallocate struct for saving zero crossing info (needs to be struct b/c
%  number of zero crossings not the same across legs)
for i = 1:length(legInd)
    zeroXing.neg2Pos.(legNames{legInd(i)}) = [];
    zeroXing.pos2Neg.(legNames{legInd(i)}) = [];
end
    
% loop through all moving bouts
for i = 1:length(moveBoutStarts)
    boutInd = (moveBoutStarts(i):moveBoutEnds(i))';
    
    % loop through all legs
    for j = 1:length(legInd)
        % find indicies of all positive velocities
        posInd = find(legXVel(boutInd,legInd(j)) >= 0) + boutInd(1);
        
        % find indicies of negative velocities
        negInd = find(legXVel(boutInd,legInd(j)) < 0) + boutInd(1);
        
        % moving bout must have positive and negative velocities to count
        if (~isempty(posInd) && ~isempty(negInd))
            % find bout starts and ends
            [posBoutStarts, posBoutEnds, posBoutDur] = findBouts(posInd);
            [negBoutStarts, negBoutEnds, negBoutDur] = findBouts(negInd);

            % find indicies that bracket zero crossing
            % negative to postive -> negBoutEnd + 1 = posBoutStart; positive
            %  bout duration defined by boutStart must be > 2
            % postive to negatvie -> posBoutEnd + 1 = negBoutStart; negative
            %  bout duration defined by boutStart must be > 2
            neg2PosZeroXingLog = ismember(posBoutStarts, (negBoutEnds+1));
            neg2PosZeroXingEnds = ...
                posBoutStarts(posBoutDur(neg2PosZeroXingLog) > minStepDur);
            neg2PosZeroXingStarts = neg2PosZeroXingEnds - 1;

            pos2NegZeroXingLog = ismember(negBoutStarts, (posBoutEnds+1));
            pos2NegZeroXingEnds = ...
                negBoutStarts(negBoutDur(pos2NegZeroXingLog) > minStepDur);
            pos2NegZeroXingStarts = pos2NegZeroXingEnds - 1;

            % append to existing
            zeroXing.neg2Pos.(legNames{legInd(j)}) = [...
                zeroXing.neg2Pos.(legNames{legInd(j)}); ...
                neg2PosZeroXingStarts];
            zeroXing.pos2Neg.(legNames{legInd(j)}) = [...
                zeroXing.pos2Neg.(legNames{legInd(j)}); ...
                pos2NegZeroXingStarts];
        end
    end
end

%% determine stance vs. swing
% for each frame, compare leg movement direction to FicTrac movement
% direction

% determine stance vs. swing: only when fly is moving, compare to direction
% of fictrac forward velocity
% interpolate fictrac forward velocity to legVidFrametimes
interpFtFwd = interp1(fictrac.t, fictrac.fwdVel, legVidFrameTimes);

% binarize fictrac forward velocity
ftFwdInd = interpFtFwd >= 0; % logical
% convert to +1 for moving forward, -1 for moving backwards
ftFwd = zeros(size(legXVel,2),1);
ftFwd(ftFwdInd) = 1; % fwd = 1
ftFwd(~ftFwdInd) = -1; % not fwd = -1

% binarize leg velocity
binLegXVel = legXVel >=0;
% covert to +1 for moving forward, -1 for moving backwards
legXVelDir = zeros(size(legXVel));
legXVelDir(binLegXVel) = 1;
legXVelDir(~binLegXVel) = -1;

legInd = [1:6]; % indicies for legs
% preallocate stance and swing matricies
legStance = false(size(legXShiftRotNormSmo,1), length(legInd));
legSwing = false(size(legXShiftRotNormSmo,1), length(legInd));
% non-logical that indicates swing, stance, or not moving in 1 vector: 1 is
%  stance, -1 is swing, 0 is not moving
legSwingStanceNotMove = zeros(size(legXShiftRotNormSmo,1), length(legInd));

% stance is leg and fictrac moving in same direction (multiply to 1), swing
%  is leg and fictrac moving in opposite directions (multiply to -1)
% assign not moving bouts to be 0
% loop through all legs
for i = 1:length(legInd)
    legFtDir = legXVelDir(:,legInd(i)) .* ftFwd;
    legFtDir(zeroVelInd) = 0;
    
    legStance(legFtDir==1,legInd(i))=true;
    legSwing(legFtDir==-1,legInd(i))=true;
    
    legSwingStanceNotMove(:,legInd(i)) = legFtDir;
    
end

%%
% correct for noise in swing, stance assignment by setting minimum bout
%  length

% % cutoff for min bout length, in samples
% minBoutLenSamp = 5; % 5 samples, ~22.5 ms
%     
% for i = 1:length(legInd)
%     % indicies when transitions between stance, swing, not moving happened
%     transInd = find(diff(legSwingStanceNotMove(:,legInd(i))) ~=0);
%     % duration in samples ofeach swing, stance, not moving bout
%     boutDur = diff(transInd);
% 
%     % deal with bouts that are too short (less than minBoutLen) - merge
%     % with other, adjacent short bouts or merge into longer sequence
%     whichBout = 1;
%     while whichBout<=length(boutDur)
%         % bout is too short
%         if (boutDur(whichBout) < minBoutLenSamp)
%             % index into transInd for start of bout
%             boutStartInd = whichBout;
% 
%             % continue until bouts stop being too short
%             for k = (whichBout+1):length(boutDur)
%                 if (boutDur(k) < minBoutLenSamp)
%                     whichBout = k;
%                 else
%                     break;
%                 end
%             end
% 
%             % index into moveLogical for bout transitions
%             boutStartMLInd = transInd(boutStartInd);
%             % index into tranInd for end of bout
%             boutEndInd = whichBout +1;
%             % equivalent for index into moveLogical
%             boutEndMLInd = transInd(boutEndInd) - 1;
% 
%             % which kind of bout (swing, stance, not moving) for first short
%             % bout
%             boutType = legSwingStanceNotMove(boutStartMLInd,legInd(i));
% 
%             % multiple short bouts, to be merged
%             if (whichBout ~= boutStartInd)
%                 % assign all of these short bouts to 1 longer bout, of the
%                 % same type as the first bout
%                 legSwingStanceNotMove(...
%                     boutStartMLInd:boutEndMLInd,legInd(i)) = boutType;
%             % one short bout, type of prior bout
%             else
%                 priorBoutType = legSwingStanceNotMove(boutStartMLInd - 1);
%                 legSwingStanceNotMove(...
%                     boutStartMLInd:boutEndMLInd,legInd(i)) = ...
%                     priorBoutType;
%             end
%         end
%         whichBout = whichBout + 1;
%     end
% 
% end

%% determine start and end indicies of half step
for i = 1:length(legInd)
    stepInd.front2Back.(legNames{legInd(i)}) = [];
    stepInd.back2Front.(legNames{legInd(i)}) = [];
end

% neg2pos is front-most position; pos2neg is back-most position
for i = 1:length(legInd)
    for j = 1:length(zeroXing.pos2Neg.(legNames{legInd(i)}))
        curInd = zeroXing.pos2Neg.(legNames{legInd(i)})(j);
        midInd = zeroXing.neg2Pos.(legNames{legInd(i)})(find(...
            zeroXing.neg2Pos.(legNames{legInd(i)}) > curInd, 1, 'first'));
        endInd = zeroXing.pos2Neg.(legNames{legInd(i)})(find(...
            zeroXing.pos2Neg.(legNames{legInd(i)}) > curInd, 1, 'first'))...
            - 1;
        
        % save into struct
        if (~isempty(midInd) && ~isempty(endInd))
            stepInd.back2Front.(legNames{legInd(i)}) = ...
                [stepInd.back2Front.(legNames{legInd(i)}); ...
                curInd (midInd-1)];
            stepInd.front2Back.(legNames{legInd(i)}) = ...
                [stepInd.front2Back.(legNames{legInd(i)}); ...
                midInd endInd];
        end
    end
end

for i = 1:length(legInd)
    meanStepVel.front2Back.(legNames{legInd(i)}) = zeros(...
        size(stepInd.front2Back.(legNames{legInd(i)}),1),1);
    meanStepVel.back2Front.(legNames{legInd(i)}) = zeros(...
        size(stepInd.back2Front.(legNames{legInd(i)}),1),1);
    for j = 1:size(stepInd.front2Back.(legNames{legInd(i)}),1)
        meanStepVel.front2Back.(legNames{legInd(i)})(j) = ...
            mean(legXVel(stepInd.front2Back.(...
            legNames{legInd(i)})(j,1):...
            stepInd.front2Back.(legNames{legInd(i)})(j,2),legInd(i)));
        meanStepVel.back2Front.(legNames{legInd(i)})(j) = ...
            mean(legXVel(stepInd.back2Front.(...
            legNames{legInd(i)})(j,1):...
            stepInd.back2Front.(legNames{legInd(i)})(j,2),legInd(i)));
    end
    fprintf('Leg %s F2B: %d\n', legNames{legInd(i)}, ...
        mean(meanStepVel.front2Back.(legNames{legInd(i)})));
    fprintf('Leg %s B2F: %d\n', legNames{legInd(i)}, ...
        mean(meanStepVel.back2Front.(legNames{legInd(i)})));
end

%% determine stance and swing
% for each step (defined by zero crossings), compare leg movement direction
%  to FicTrac direction (i.e. average over frames)


%% phase of each leg (in X dimension), using Hilbert transform
% preallocate
legXPhase = zeros(size(legXShiftRotNorm,1), length(legInd));
legXPhaseSmo = zeros(size(legXShiftRotNormSmo,1), length(legInd));
% get leg phase on both smoothed and not smoothed leg position
for i = 1:length(legInd)
    legXPhase(:,legInd(i)) = angle(hilbert(legXShiftRotNorm(:,legInd(i))));
    legXPhaseSmo(:,legInd(i)) = ...
        angle(hilbert(legXShiftRotNormSmo(:,legInd(i))));
end

%% direction each leg moves in xy plane

% direction at each time point
% preallocate; 1 less after difference
moveAngFrames = zeros(size(legXShiftRotNorm, 1) - 1, length(legInd));

% compute direction moved for each leg
for i = 1:length(legInd)
    moveAngFrames(:, legInd(i)) = findAngle2Pts(...
        legXShiftRotNorm(1:(end-1),legInd(i)), ...
        legYShiftRotNorm(1:(end-1),legInd(i)), ...
        legXShiftRotNorm(2:end, legInd(i)), ...
        legYShiftRotNorm(2:end, legInd(i)));
end

% direction for each half step, as defined by start and end pts
% preallocate
for i = 1:length(legInd)
    moveAngSteps.front2Back.(legNames{legInd(i)}) = ...
        zeros(size(stepInd.front2Back.(legNames{legInd(i)}),1),1);
    moveAngSteps.back2Front.(legNames{legInd(i)}) = ...
        zeros(size(stepInd.back2Front.(legNames{legInd(i)}),1),1);
end
% compute direction moved for each leg
for i = 1:length(legInd)
    moveAngSteps.front2Back.(legNames{legInd(i)}) = findAngle2Pts(...
        legXShiftRotNorm(stepInd.front2Back.(legNames{legInd(i)})(:,1),...
        legInd(i)), ...
        legYShiftRotNorm(stepInd.front2Back.(legNames{legInd(i)})(:,1),...
        legInd(i)), ...
        legXShiftRotNorm(stepInd.front2Back.(legNames{legInd(i)})(:,2),...
        legInd(i)), ...
        legYShiftRotNorm(stepInd.front2Back.(legNames{legInd(i)})(:,2),...
        legInd(i)));
    moveAngSteps.back2Front.(legNames{legInd(i)}) = findAngle2Pts(...
        legXShiftRotNorm(stepInd.back2Front.(legNames{legInd(i)})(:,1),...
        legInd(i)), ...
        legYShiftRotNorm(stepInd.back2Front.(legNames{legInd(i)})(:,1),...
        legInd(i)), ...
        legXShiftRotNorm(stepInd.back2Front.(legNames{legInd(i)})(:,2),...
        legInd(i)), ...
        legYShiftRotNorm(stepInd.back2Front.(legNames{legInd(i)})(:,2),...
        legInd(i)));
end

%% random plots

% scatterplot of midpoints
midX = legX(:,9);
midY = legY(:,9);

figure;
scatter(midX, midY, '.', 'MarkerFaceAlpha', 0.2, 'MarkerEdgeAlpha', 0.2);
xlim([0 legImgHorizPx]);
ylim([0 legImgVertPx]);

% colormap for points
cmap = lines(size(legX,2));
% swap out last 4 as same color for body
bodyCol = cmap(1,:);
cmap(1:6,:) = cmap(2:7,:);
cmap(7:end,:) = repmat(bodyCol,5,1);

figure;
imagesc(legVidImg,[legVidMinInt legVidMaxInt]);
colormap('gray');
hold on;

figure;
hold on;
scatter(legX(legTrackInd,:),legY(legTrackInd,:),50,cmap,'x', ...
    'LineWidth', 1.5);
plot(legX(legTrackInd,lineFitInd), ...
    polyval(bodyLineCoeffs(legTrackInd,:),...
    legX(legTrackInd,lineFitInd)), 'r');
scatter(midXProj(legTrackInd), midYProj(legTrackInd),'g');
scatter(headXProj(legTrackInd), headYProj(legTrackInd), 'r');
scatter(abdXProj(legTrackInd), abdYProj(legTrackInd), 'b');
axis square;

figure;
hold on;
scatter(legXShift(legTrackInd,:),legYShift(legTrackInd,:),50,cmap,...
    'x', 'LineWidth', 1.5);
scatter(legXShiftRot(legTrackInd,:),legYShiftRot(legTrackInd,:),30,...
    cmap,'x', 'LineWidth', 1.5);
axis square;

% plot red dots for when fly stopped, on top of R2 leg velocity
figure; 
plot(legVidFrameTimes, legXVel(:,2)); 
hold on; 
plot(legVidFrameTimes(zeroVelInd), legXVel(zeroVelInd,2), '.');

% plot leg position trajectories, with red dots for stance
figure; 
plot(legVidFrameTimes, legXShiftRotNorm(:,legInd)); 
hold on; 
for i = 1:length(legInd)
    plot(legVidFrameTimes(legStance(:,legInd(i))), ...
        legXShiftRotNorm(legStance(:,legInd(i)),legInd(i)),'.r');
end

% plot swing/stance as heatmap
figure;
imagesc(legSwingStanceNotMove');
xlim([10000 11000]);

% plot velocity zero crossings on leg positions
figure;
plot(legVidFrameTimes, legXShiftRotNorm(:,2));
hold on;
plot(legVidFrameTimes(zeroXing.neg2Pos.R2),...
    legXShiftRotNorm(zeroXing.neg2Pos.R2,2), '.');
plot(legVidFrameTimes(zeroXing.pos2Neg.R2),...
    legXShiftRotNorm(zeroXing.pos2Neg.R2,2), '.');

figure;
plot(legVidFrameTimes, legXShiftRotNorm(:,1));
hold on;
plot(legVidFrameTimes(zeroXing.neg2Pos.R1),...
    legXShiftRotNorm(zeroXing.neg2Pos.R1,1), '.');
plot(legVidFrameTimes(zeroXing.pos2Neg.R1),...
    legXShiftRotNorm(zeroXing.pos2Neg.R1,1), '.');

figure;
plot(legVidFrameTimes, legXShiftRotNorm(:,3));
hold on;
plot(legVidFrameTimes(zeroXing.neg2Pos.R3),...
    legXShiftRotNorm(zeroXing.neg2Pos.R3,3), '.');
plot(legVidFrameTimes(zeroXing.pos2Neg.R3),...
    legXShiftRotNorm(zeroXing.pos2Neg.R3,3), '.');

% leg phase with zero crossings
figure;
plot(legVidFrameTimes, legXPhase(:,2));
hold on;
plot(legVidFrameTimes(zeroXing.neg2Pos.R2),...
    legXPhase(zeroXing.neg2Pos.R2,2), '.');
plot(legVidFrameTimes(zeroXing.pos2Neg.R2),...
    legXPhase(zeroXing.pos2Neg.R2,2), '.');

% plot leg phase
figure;
plot(legVidFrameTimes, legXPhase(:,2));
hold on;
plot(legVidFrameTimes, legXPhaseSmo(:,2));


% plot leg position scatters
figure;
scatter(legYShiftRotNorm(stepInd.front2Back.R1(:,1),1),...
    legXShiftRotNorm(stepInd.front2Back.R1(:,1),1),'.');
hold on;
scatter(legYShiftRotNorm(stepInd.front2Back.R1(:,2),1),...
    legXShiftRotNorm(stepInd.front2Back.R1(:,2),1),'.');

figure;
scatter(legYShiftRotNorm(stepInd.front2Back.L1(:,1),4),...
    legXShiftRotNorm(stepInd.front2Back.L1(:,1),4),'.');
hold on;
scatter(legYShiftRotNorm(stepInd.front2Back.L1(:,2),4),...
    legXShiftRotNorm(stepInd.front2Back.L1(:,2),4),'.');
