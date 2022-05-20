%% for each step, get spike rate; variable time offset
% for each half step

% delays in sec, neg delay is ephys b/f behavior
tDelay = [-0.5 -0.2 -0.1 -0.05 -0.025 0 0.025 0.05 0.1 0.2 0.5]; 

% preallocate for spike rate matrix for (number of steps x 2 (for 2 half
%  steps) x num delays
stepSpikeRate = zeros(size(legSteps.stepInds,1),2, length(tDelay));

% times when all spikes occured
spikeTimes = ephysSpikes.t(ephysSpikes.startInd);


% loop through all time delays
for j = 1:length(tDelay)
    % incorporate time offset b/w ephys and behavior
    spikeTimesDelay = spikeTimes - tDelay(j);

    % loop through all steps
    for i = 1:size(legSteps.stepInds,1)
        % get step start, mid, end times
        stepStartT = legTrack.t(legSteps.stepInds(i,1));
        stepMidT = legTrack.t(legSteps.stepInds(i,2));
        stepEndT = legTrack.t(legSteps.stepInds(i,3));

        % for first half step, figure out how many spikes 
        numSpikesStartMid = sum((spikeTimesDelay >= stepStartT) &...
            (spikeTimesDelay < stepMidT));
        % convert to step rate
        stepSpikeRate(i,1,j) = numSpikesStartMid / (stepMidT-stepStartT);

        % for second half step, figure out how many spikes 
        numSpikesMidEnd = sum((spikeTimesDelay >= stepMidT) &...
            (spikeTimesDelay < stepEndT));
        % convert to step rate
        stepSpikeRate(i,2,j) = numSpikesMidEnd / (stepEndT-stepMidT); 
    end
end

%% for each step, get average Vm; variable time offset
%  for each half step

% delays in sec, neg delay is ephys b/f behavior
tDelay = [-0.5 -0.2 -0.1 -0.05 -0.025 0 0.025 0.05 0.1 0.2 0.5]; 

% preallocate for spike rate matrix for (number of steps x 2 (for 2 half
%  steps) x num delays
stepVmRate = zeros(size(legSteps.stepInds,1),2, length(tDelay));

for j = 1:length(tDelay)
    
    % incorporate time offset b/w ephys and behavior
    vmTDelay = ephysSpikes.t - tDelay(j);
    
    % loop through all steps
    for i = 1:size(legSteps.stepInds,1)
        % get step start, mid, end times
        stepStartT = legTrack.t(legSteps.stepInds(i,1));
        stepMidT = legTrack.t(legSteps.stepInds(i,2));
        stepEndT = legTrack.t(legSteps.stepInds(i,3));
        
        % for first half step, get mean Vm
        % find ind into vmDelay of step start, mid
        startInd = find(vmTDelay >= stepStartT,1,'first');
        midInd = find(vmTDelay >= stepMidT,1,'first');
        
        thisMeanVm = mean(ephysSpikes.medFiltV(startInd:midInd));
        
        % convert
        stepVmRate(i,1,j) = thisMeanVm;
        
        % for second half step
        midInd = find(vmTDelay >= stepMidT,1,'first');
        endInd = find(vmTDelay >= stepEndT, 1, 'first');
        
        thisMeanVm = mean(ephysSpikes.medFiltV(midInd:endInd));
        
        stepVmRate(i,2,j) = thisMeanVm;
        
    end
end

%% stance step speed plotted on same plot for legs on each side

% which parameter, comment out as needed

% thisStepParam = legSteps.stepSpeeds(:,2);
thisStepParam = legSteps.stepLengths(:,2);
whichParamStr = 'Stance Step Speed (body lengths/s)';
thisEphysSpikes = stepSpikeRate;
thisWhichLeg = legSteps.stepWhichLeg;

xLimits = [0 0.75];
yLimits = [0 150];
% colormap
colMap = lines(6);

% loop through all tDelay
for j = 4
    
    % spikes just for this delay
    thisDelayEphysSpikes = thisEphysSpikes(:,j);
    
    figure('Position', [10 10 1000 400]);
    
    % ipsi legs (right)
    subplot(1,2,1);
    
    for i = 1:3
        thisLeg = legIDs.ind(i);
        % steps for this leg
        thisLegStepParams = thisStepParam(thisWhichLeg == thisLeg);
        % ephys for steps for this leg
        thisLegStepSpikes = thisDelayEphysSpikes(thisWhichLeg == thisLeg);
        
        % remove outliers
        [~,stepParamOutInd] = rmoutliers(thisLegStepParams);
        stepParamNoOut = thisLegStepParams(~stepParamOutInd);
        spikeNoOut = thisLegStepSpikes(~stepParamOutInd);
        
        % fit line
        p = polyfit(stepParamNoOut, spikeNoOut, 1);
        f = polyval(p,linspace(xLimits(1),xLimits(2)));

        scatter(stepParamNoOut,spikeNoOut, 80, '.', ...
            'MarkerEdgeColor', colMap(i,:), 'MarkerFaceColor', colMap(i,:));
        hold on;
        plot(linspace(xLimits(1),xLimits(2)),f, 'Color', colMap(i,:), ...
            'LineWidth', 2);
        legend;
    end
    ylim(yLimits);
    xlim(xLimits);
    
    % contra legs (left)
    subplot(1,2,2);
    for i = 4:6
        thisLeg = legIDs.ind(i);
        % steps for this leg
        thisLegStepParams = thisStepParam(thisWhichLeg == thisLeg);
        % ephys for steps for this leg
        thisLegStepSpikes = thisDelayEphysSpikes(thisWhichLeg == thisLeg);
        
        % remove outliers
        [~,stepParamOutInd] = rmoutliers(thisLegStepParams);
        stepParamNoOut = thisLegStepParams(~stepParamOutInd);
        spikeNoOut = thisLegStepSpikes(~stepParamOutInd);
        
        % fit line
        p = polyfit(stepParamNoOut, spikeNoOut, 1);
        f = polyval(p,linspace(xLimits(1),xLimits(2)));

        scatter(thisLegStepParams,thisLegStepSpikes, 80, '.',...
            'MarkerEdgeColor', colMap(i,:), 'MarkerFaceColor', colMap(i,:));
        hold on;
        plot(linspace(xLimits(1),xLimits(2)),f, 'Color', colMap(i,:), ...
            'LineWidth', 2);  
        legend;
    end
    ylim(yLimits);
    xlim(xLimits);

end
