%% get data from multiple pData files

% delays in sec, neg delay is ephys b/f behavior
tDelay = [-0.5 -0.2 -0.1 -0.05 -0.025 0 0.025 0.05 0.1 0.2 0.5]; 

% user selects pData files

% prompt user to select pData files
[pDataFNames, pDataPath] = uigetfile('*.mat', 'Select pData files', ...
    pDataDir(), 'MultiSelect', 'on');

% if only 1 pData file selected, not cell array; make sure loop still
%  works 
if (iscell(pDataFNames))
    numPDataFiles = length(pDataFNames);
else
    numPDataFiles = 1;
end


for k = 1:numPDataFiles

    % handle whether it's a cell array or not
    if (iscell(pDataFNames))
        pDataName = pDataFNames{k};
    else
        pDataName = pDataFNames;
    end

    pDataFullPath = [pDataPath pDataName];
    
    % load data
    load(pDataFullPath, 'legSteps','ephysSpikes','legTrack');
    
    legStepsAll{k} = legSteps;
    ephysSpikesAll{k} = ephysSpikes;
    legTrackAll{k} = legTrack;
end

% loop through all trials

for k = 1:numPDataFiles

    % preallocate for spike rate matrix for (number of steps x 2 (for 2 half
    %  steps) x num delays
    stepSpikeRate = zeros(size(legStepsAll{k}.stepInds,1),2, length(tDelay));

    % times when all spikes occured
    spikeTimes = ephysSpikesAll{k}.t(ephysSpikesAll{k}.startInd);


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
    
    stepSpikeRateAll{k} = stepSpikeRate;

end

%% for each step parameter, plot scatter against spike rate
% subplots for each leg

thisStepParam = 'stepLengths';
whichParamStr = 'Step Lengths (Body Lengths)';

% loop through all delay times
for j = 1:length(tDelay)
    
    figure('Position', [10 10 1200 900]);
    suptitle(sprintf('%s, delay %d ms',whichParamStr,tDelay(j) * 1000));
    % 
    for i = 1:length(legIDs.ind)
        % initialize variables
        thisLegStepParamsAll = [];
        thisLegStepSpikesAll = [];
        
        for k = 1:numPDataFiles
            thisLeg = legIDs.ind(i);
            
            % this delay stepSpikeRate
            thisDelayStepSpikes = stepSpikeRateAll{k}(:,:,j);
            
            thisLegSteps = legStepsAll{k};
            thisParam = thisLegSteps.(thisStepParam);
            thisWhichLeg = thisLegSteps.stepWhichLeg;
    
            % steps for this leg
            thisLegStepParams = thisParam(thisWhichLeg == thisLeg, :);
            % ephys for steps for this leg
            thisLegStepSpikes = thisDelayStepSpikes(thisWhichLeg == thisLeg, :);


    %         step2Coeff = corrcoef(thisLegStepParams(:,2),thisLegStepSpikes(:,2));
        
        
            thisLegStepParamsAll = [thisLegStepParamsAll; thisLegStepParams];
            thisLegStepSpikesAll = [thisLegStepSpikesAll; thisLegStepSpikes];
        end
        
            % get correlation coefficients
            corrVal = corrcoef(thisLegStepParamsAll(:,2),thisLegStepSpikesAll(:,2));

            subplot(2,3,i)
            scatter(thisLegStepParamsAll(:,2),thisLegStepSpikesAll(:,2), 50, '.');
        
   
    %     hold on;
    %     scatter(thisLegStepParams(:,2),thisLegStepSpikes(:,2), 50,'.');

%         title(sprintf('Leg %s, r=%0.3f, %0.3f', ...
%             zeroXingParams.legNames{i}, corrVal(1,2), ...
%             step2Coeff(1,2)));
        title(sprintf('Leg %s, r=%0.3f', ...
            legIDs.names{i}, corrVal(1,2)));

        ylabel('Spike Rate (spikes/sec)');
    end
end



%% for each step parameter, plot scatter against spike rate
% subplots for each leg

% change this up to change what is being plotted
% thisStepParam = legSteps.stepLengths;
% whichParamStr = 'Step Lengths (Body Lengths)';

thisStepParam = legSteps.stepXLengths;
whichParamStr = 'Step X Lengths (Body Lengths)';

% stepDirsShift = wrapTo180(legSteps.stepDirections+90);
% thisStepParam = stepDirsShift;
% whichParamStr = 'Step Directions (deg)';

% thisStepParam = legSteps.stepDurations;
% whichParamStr = 'Step Durations (sec)';
% stepsWhichLeg = legSteps.stepWhichLeg;

% thisStepParam = legSteps.stepSpeeds;
% whichParamStr = 'Step Speeds (Body Lengths/sec)';
stepsWhichLeg = legSteps.stepWhichLeg;

% thisStepParam = legSteps.stepVelX;
% whichParamStr = 'Step X Velocities (Body Lengths/sec)';

% thisStepParam = legSteps.stepVelY;
% whichParamStr = 'Step Y Velocities (Body Lengths/sec)';

% thisStepParam = legSteps.stepFtFwd;
% whichParamStr = 'FicTrac Fwd Vel (mm/sec)';

% thisStepParam = stepFtLat;
% whichParamStr = 'FicTrac Lateral Vel (mm/sec)';

% thisStepParam = legSteps.stepFtYaw;
% whichParamStr = 'FicTrac Yaw Vel (deg/sec)';

% thisStepParam = strideSpeed;
% whichParamStr = 'Stride Speed (Body Lengths/s)';

% loop through all delay times
for j = 1:length(tDelay)
    % this delay stepSpikeRate
    thisDelayStepSpikes = stepSpikeRate(:,:,j);
    
    figure('Position', [10 10 1200 900]);
    suptitle(sprintf('%s, delay %d ms',whichParamStr,tDelay(j) * 1000));
    % 
    for i = 1:length(legIDs.ind)
        thisLeg = legIDs.ind(i);
        % steps for this leg
        thisLegStepParams = thisStepParam(stepsWhichLeg == thisLeg, :);
        % ephys for steps for this leg
        thisLegStepSpikes = thisDelayStepSpikes(stepsWhichLeg == thisLeg, :);

        % get correlation coefficients
        corrVal = corrcoef(thisLegStepParams(:,2),thisLegStepSpikes(:,2));
%         step2Coeff = corrcoef(thisLegStepParams(:,2),thisLegStepSpikes(:,2));

        subplot(2,3,i)
        scatter(thisLegStepParams(:,2),thisLegStepSpikes(:,2), 50, '.');
    %     hold on;
    %     scatter(thisLegStepParams(:,2),thisLegStepSpikes(:,2), 50,'.');

%         title(sprintf('Leg %s, r=%0.3f, %0.3f', ...
%             zeroXingParams.legNames{i}, corrVal(1,2), ...
%             step2Coeff(1,2)));
        title(sprintf('Leg %s, r=%0.3f', ...
            legIDs.names{i}, corrVal(1,2)));

        ylabel('Spike Rate (spikes/sec)');
    end
end

%% Scatterplot for two behavioral variables

% param1 = legSteps.stepSpeeds;
% param1Str = 'Step Speed (body lengths/s)';

% param1 = legSteps.stepLengths;
% param1Str = 'Step Lengths (body lengths/s)';

param1 = legSteps.stepXLengths;
param1Str = 'Step X Lengths (body lengths/s)';

% param1 = legSteps.stepVelX;
% param1Str = 'Step X Velocity (body lengths/s)';

% param1 = legSteps.stepVelY;
% param1Str = 'Step Y Velocity (body lengths/s)';

param2 = legSteps.stepFtYaw;
param2Str = 'FicTrac Yaw Vel (deg/s)';

% param2 = legSteps.stepFtFwd;
% param2Str = 'FicTrac Fwd Vel (mm/s)';

figure;
suptitle(sprintf('%s vs %s', param2Str, param1Str));
% 
for i = 1:length(legIDs.ind)
    thisLeg = legIDs.ind(i);
    % param1 for this leg
    thisLegParam1 = param1(stepsWhichLeg == thisLeg, :);
    % param2 for this leg
    thisLegParam2 = param2(stepsWhichLeg == thisLeg, :);
    
    % get correlation coefficients
    corrVal = corrcoef(thisLegParam1(:,2),thisLegParam2(:,2));
%     step2Coeff = corrcoef(thisLegParam1(:,2),thisLegParam2(:,2));
    
    subplot(2,3,i)
    scatter(thisLegParam1(:,2),thisLegParam2(:,2), 50, '.');
%     hold on;
%     scatter(thisLegParam1(:,2),thisLegParam2(:,2), 50,'.');
    
%     title(sprintf('Leg %s, r=%0.3f, %0.3f', ...
%         zeroXingParams.legNames{i}, corrVal(1,2), ...
%         step2Coeff(1,2)));
    
    title(sprintf('Leg %s, r=%0.3f', ...
        legIDs.names{i}, corrVal(1,2)));
    
    ylabel(param2Str);
    xlabel(param1Str);
end