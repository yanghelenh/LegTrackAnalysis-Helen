%% plot swing/stance, leg positions, ephys voltage, FicTrac velocities
f1XLims = [234 240];
figure;
ax1 = subplot(6,1,1);
plot(leg.frameTimes, legTrack.srnLegX(:,1:6));
% hold on; 
% plot(leg.frameTimes(zeroVelInd), srnLegX(zeroVelInd,2), '.');
legend(legIDs.names);
title('Normalized leg position, front-back axis');
ylabel('Body lengths');


ax2 = subplot(6,1,2);
plot(frameTimes, legTrack.srnLegY(:,1:6));
legend(legIDs.names);
title('Normalized leg position, side axis');
ylabel('Body lengths');

ax3 = subplot(6,1,3);
% swingStanceLims = (f1XLims-leg.frameTimes(1)) .* ...
%    (1/median(diff(leg.frameTimes)));
swingStanceXlims = [frameTimes(1) frameTimes(end)];
swingStanceYlims = [1 6];
imagesc(swingStanceXlims, swingStanceYlims, framesSwingStance');
title('Swing/stance calls (yellow - stance, blue - swing, teal - not moving');
ylabel({'R1';'R2';'R3';'L1';'L2';'L3'});
ylh = get(gca,'ylabel');
gyl = get(ylh);                                                        
ylp = get(ylh, 'Position');
set(ylh, 'Rotation',0, 'Position',ylp, 'VerticalAlignment','middle', ...
    'HorizontalAlignment','right');

xlim(swingStanceXlims);

ax4 = subplot(6,1,4);
plot(ephysData.t, ephysData.scaledVoltage);
title('Scaled voltage');
ylabel('mV');

% ax4 = subplot(6,1,4);
% plot(ephysSpikes.t, ephysSpikes.spikeRate);

ax5 = subplot(6,1,5);
plot(fictracProc.t, fictracProc.yawAngVel);
title('Yaw velocity');
ylabel('deg/s');

ax6 = subplot(6,1,6);
plot(fictracProc.t, fictracProc.fwdVel);
title('Forward velocity');
ylabel('mm/s');

linkaxes([ax1 ax2 ax3 ax4 ax5 ax6], 'x');
xlim(f1XLims);

%%

phaseFieldnames = fieldnames(phaseDiffs);

whichInds = 4:6;

figure;

for i = whichInds
    % for plotting, repeat 1st point at end to connect
    theta = [phaseBinMids phaseBinMids(1)];
    rho = [phasePDFs.(phaseFieldnames{i}) ...
        phasePDFs.(phaseFieldnames{i})(1)];
    
    polarplot(deg2rad(theta), rho);
    hold on;
end

legend(phaseFieldnames(whichInds));
title('Leg phase relationships');

%% scatter plot phase difference vs. yaw/fwd velocity

% interpolate FicTrac yaw velocity to leg track time scale
ftYawVelLegT = interp1(fictracProc.t, fictracProc.yawAngVel, legTrack.t);
ftFwdVelLegT = interp1(fictracProc.t, fictracProc.fwdVel, legTrack.t);

% scatterplot against each phase variable

phaseFieldnames = fieldnames(phaseDiffs);

% yaw, loop through all
for i = 1:length(phaseFieldnames)
    % remove NaNs from phase diff
    thisPhaseDiff = phaseDiffs.(phaseFieldnames{i});
    notNanLog = ~isnan(thisPhaseDiff);
    thisPhaseDiffNoNaN = thisPhaseDiff(notNanLog);
    thisYawVelNoNaN = ftYawVelLegT(notNanLog);
    
    figure;
    s = scatter(thisPhaseDiffNoNaN,thisYawVelNoNaN, 10, 'filled');
    s.MarkerFaceAlpha = 0.1;
    s.MarkerEdgeAlpha = 0;
    
    % correlation coefficient
    corrVal = corrcoef(thisPhaseDiffNoNaN,thisYawVelNoNaN);
    
    title(sprintf('%s phase diff vs. yaw vel, r=%0.3f', ...
        phaseFieldnames{i}, corrVal(1,2)));
    
end

%% fwd, loop through all
for i = 1:length(phaseFieldnames)
    % remove NaNs from phase diff
    thisPhaseDiff = phaseDiffs.(phaseFieldnames{i});
    notNanLog = ~isnan(thisPhaseDiff);
    thisPhaseDiffNoNaN = thisPhaseDiff(notNanLog);
    thisFwdVelNoNaN = ftFwdVelLegT(notNanLog);
    
    figure;
    s = scatter(thisFwdVelNoNaN, thisPhaseDiffNoNaN,10, 'filled');
        s.MarkerFaceAlpha = 0.1;
    s.MarkerEdgeAlpha = 0;
    
    % correlation coefficient
    corrVal = corrcoef(thisPhaseDiffNoNaN,thisFwdVelNoNaN);
    
    title(sprintf('%s phase diff vs. fwd vel, r=%0.3f', ...
        phaseFieldnames{i}, corrVal(1,2)));
    
end

%% plot phase diff against spike rate

delay = 0;
spikeRate = getSpikeRateSpecTimeWithDelay(...
    ephysSpikes.t(ephysSpikes.startInd), legTrack.t, delay);

phaseFieldnames = fieldnames(phaseDiffs);

for i = 1:length(phaseFieldnames)
    % remove NaNs from phase diff
    thisPhaseDiff = phaseDiffs.(phaseFieldnames{i});
    notNanLog = ~isnan(thisPhaseDiff);
    thisPhaseDiffNoNaN = thisPhaseDiff(notNanLog);
    thisSpikerateNoNaN = spikeRate(notNanLog);
    
    figure; 
    s = scatter(thisPhaseDiffNoNaN, thisSpikerateNoNaN, 10, 'filled');
%     s.Marker = '.';
    s.MarkerFaceAlpha = 0.1;
    s.MarkerEdgeAlpha = 0;
    
    % correlation coefficient
    corrVal = corrcoef(thisPhaseDiffNoNaN, thisSpikerateNoNaN);
    
    title(sprintf('%s phase diff vs. spike rate, delay=%d ms, r=%0.3f',...
        phaseFieldnames{i}, delay * 1000, corrVal(1,2)));
end

%% plot phase against spike rate
delay = 0;
spikeRate = getSpikeRateSpecTimeWithDelay(...
    ephysSpikes.t(ephysSpikes.startInd), legTrack.t, delay);

% loop through all legs
for i = 1:size(legPhase,2)
    % remove NaNs from phase
    notNanLog = ~isnan(legPhase(:,i));
    thisLegNotNaNPhase = legPhase(notNanLog,i);
    thisSpikerateNoNaN = spikeRate(notNanLog);
    
    figure;
    s = scatter(thisLegNotNaNPhase, thisSpikerateNoNaN, 10, 'filled');
%     s.Marker = '.';
    s.MarkerFaceAlpha = 0.1;
    s.MarkerEdgeAlpha = 0;
    
    % correlation coefficient
    corrVal = corrcoef(thisLegNotNaNPhase, thisSpikerateNoNaN);
    
    title(sprintf('Leg %s phase vs spike rate, delay=%d ms, r=%0.3f',...
        legIDs.names{i}, delay * 1000, corrVal(1,2)));
end

%% get mean spike rate for each phase
delay = 0;
spikeRate = getSpikeRateSpecTimeWithDelay(...
    ephysSpikes.t(ephysSpikes.startInd), legTrack.t, delay);

numPhaseBins = 120;
binMinVal = 0;
binMaxVal = 360;
phaseBinEdges = linspace(binMinVal,binMaxVal, numPhaseBins + 1);
phaseBinMids = zeros(1,numPhaseBins);
for i = 1:numPhaseBins
    phaseBinMids(i) = (phaseBinEdges(i+1) - phaseBinEdges(i))/2 + ...
        phaseBinEdges(i);
end
phaseBinStarts = phaseBinEdges(1:(end-1));

binCounts = zeros(numPhaseBins,6);
binSums = zeros(numPhaseBins,6);
    
% loop through all legs
for i = 1:6
    % remove NaNs from phase
    notNanLog = ~isnan(legPhase(:,i));
    thisLegNotNaNPhase = legPhase(notNanLog,i);
    thisSpikerateNoNaN = spikeRate(notNanLog);
    
    % loop through all phase vals
    for j = 1:length(thisLegNotNaNPhase)
        % find bin this phase belongs to
        binInd = find(thisLegNotNaNPhase(j) >= phaseBinStarts, 1, 'last');
        
        
        % add spike rate to that bin
        binSums(binInd,i) = binSums(binInd,i) + thisSpikerateNoNaN(j);
        % increment counter
        binCounts(binInd,i) = binCounts(binInd,i) + 1;
    end
end

% get bin mean spike rates
binMeanSpikeRate = binSums ./ binCounts;

% plot mean spike rate vs. phase
for i = 1:6
    figure;
    plot(phaseBinMids,binMeanSpikeRate(:,i));
    legend(legIDs.names{i});
end

    figure;
    plot(phaseBinMids,binMeanSpikeRate);
    legend(legIDs.names);

    
%% get spike counts binned by phase, normalized by time spent at that phase
%  alternative approach to mean spike rate per bin

delay = 0; % offset b/w phase and spike time, in sec

numPhaseBins = 120;
binMinVal = 0;
binMaxVal = 360;
phaseBinEdges = linspace(binMinVal,binMaxVal, numPhaseBins + 1);
phaseBinMids = zeros(1,numPhaseBins);
for i = 1:numPhaseBins
    phaseBinMids(i) = (phaseBinEdges(i+1) - phaseBinEdges(i))/2 + ...
        phaseBinEdges(i);
end
phaseBinStarts = phaseBinEdges(1:(end-1));
phaseBinEnds = phaseBinEdges(2:end);

% preallocate
binSpikeCounts = zeros(numPhaseBins, 6);
binPhaseCounts = zeros(numPhaseBins, 6);

spikeTimes = ephysSpikes.t(ephysSpikes.startInd);

% loop through all legs
for i = 1:6
    % this leg phases
    thisLegPhases = legPhase(:,i);
    
    % loop through all spikes, assign to phase bin
    for j = 1:length(spikeTimes)
        % find frame during which this spike occured
        thisSpikeFrameInd = find(spikeTimes(j)-delay >= legTrack.t, 1, 'last');
        % get phase at this frame
        thisFramePhase = thisLegPhases(thisSpikeFrameInd);
        
        % find phase bin index this phase belongs to (if not NaN)
        if (~isnan(thisFramePhase))
            thisBinInd = find(thisFramePhase >= phaseBinStarts, 1, 'last');
            % increment count in appropriate bin
            binSpikeCounts(thisBinInd,i) = binSpikeCounts(thisBinInd,i) + 1;
        end
    end
    

    % loop through all phases, assign to phase bin
    for j = 1:length(thisLegPhases)
        % current pahse
        thisFramePhase = thisLegPhases(j);
        
        % find phase bin index this phase belongs to (if not NaN)
        if (~isnan(thisFramePhase))
            thisBinInd = find(thisFramePhase >= phaseBinStarts, 1, 'last');
            % increment count in appropriate bin
            binPhaseCounts(thisBinInd,i) = binPhaseCounts(thisBinInd,i) + 1;
        end
    end
end

% convert phase counts into total time
ifi = median(diff(legTrack.t));
binPhaseTimes = binPhaseCounts * ifi;

% convert spike counts and phase times into spike rate
binSpikeRate = binSpikeCounts ./ binPhaseTimes;

% plot spike rate vs. phase
for i = 1:6
    figure;
    plot(phaseBinMids,binSpikeRate(:,i));
    legend(legIDs.names{i});
end

figure;
plot(phaseBinMids,binSpikeRate);
legend(legIDs.names);
    
% plot phase bin counts
figure;
plot(phaseBinMids,binPhaseCounts);
legend(legIDs.names);



%% for each step parameter, plot scatter against spike rate
% subplots for each leg

% change this up to change what is being plotted
% thisStepParam = stepLengths;
% whichParamStr = 'Step Lengths (Body Lengths)';

% stepDirsShift = wrapTo180(stepDirections+90);
% thisStepParam = stepDirsShift;
% whichParamStr = 'Step Directions (deg)';

% thisStepParam = stepDurations;
% whichParamStr = 'Step Durations (sec)';

% thisStepParam = stepVelocities;
% whichParamStr = 'Step Speeds (Body Lengths/sec)';

% thisStepParam = stepVelX;
% whichParamStr = 'Step X Velocities (Body Lengths/sec)';

% thisStepParam = stepVelY;
% whichParamStr = 'Step Y Velocities (Body Lengths/sec)';

% thisStepParam = stepFtFwd;
% whichParamStr = 'FicTrac Fwd Vel (mm/sec)';

% thisStepParam = stepFtLat;
% whichParamStr = 'FicTrac Lateral Vel (mm/sec)';

% thisStepParam = stepFtYaw;
% whichParamStr = 'FicTrac Yaw Vel (deg/sec)';

thisStepParam = strideSpeed;
whichParamStr = 'Stride Speed (Body Lengths/s)';

% loop through all delay times
for j = 1:length(tDelay)
    % this delay stepSpikeRate
    thisDelayStepSpikes = stepSpikeRate(:,:,j);
    
    figure('Position', [10 10 1200 900]);
    suptitle(sprintf('%s, delay %d ms',whichParamStr,tDelay(j) * 1000));
    % 
    for i = 1:length(zeroXingParams.legInd)
        thisLeg = zeroXingParams.legInd(i);
        % steps for this leg
        thisLegStepParams = thisStepParam(stepsWhichLeg == thisLeg, :);
        % ephys for steps for this leg
        thisLegStepSpikes = thisDelayStepSpikes(stepsWhichLeg == thisLeg, :);

        % get correlation coefficients
        corrVal = corrcoef(thisLegStepParams(:,1),thisLegStepSpikes(:,1));
%         step2Coeff = corrcoef(thisLegStepParams(:,2),thisLegStepSpikes(:,2));

        subplot(2,3,i)
        scatter(thisLegStepParams(:,1),thisLegStepSpikes(:,1), 50, '.');
    %     hold on;
    %     scatter(thisLegStepParams(:,2),thisLegStepSpikes(:,2), 50,'.');

%         title(sprintf('Leg %s, r=%0.3f, %0.3f', ...
%             zeroXingParams.legNames{i}, corrVal(1,2), ...
%             step2Coeff(1,2)));
        title(sprintf('Leg %s, r=%0.3f', ...
            zeroXingParams.legNames{i}, corrVal(1,2)));

        ylabel('Spike Rate (spikes/sec)');
    end
end