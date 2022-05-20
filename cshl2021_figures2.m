%% for each step parameter, plot scatter against spike rate
% subplots for each leg

% change this up to change what is being plotted
% thisStepParam = legSteps.stepLengths;
% whichParamStr = 'Step Lengths (Body Lengths)';

% thisStepParam = legSteps.stepXLengths;
% whichParamStr = 'Step X Lengths (Body Lengths)';

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

thisStepParam = legSteps.stepFtFwd;
whichParamStr = 'FicTrac Fwd Vel (mm/sec)';

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

%% for each step parameter, plot scatter against mean Vm
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
    thisDelayStepSpikes = stepVmRate(:,:,j);
    
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

        ylabel('Mean Vm (mV)');
    end
end

%% Scatterplot for two behavioral variables

% param1 = legSteps.stepSpeeds;
% param1Str = 'Step Speed (body lengths/s)';

param1 = legSteps.stepLengths;
param1Str = 'Step Lengths (body lengths/s)';

% param1 = legSteps.stepXLengths;
% param1Str = 'Step X Lengths (body lengths/s)';

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