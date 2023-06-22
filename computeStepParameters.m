% computeStepParameters.m
%
% Function that takes step indices [start, mid, end] and computes step
%  parameters: step lengths, step directions, step durations, step speeds,
%  step X velocities, step Y velocities
%
% INPUTS:
%   legSteps - struct of leg step data, from processLegTrack() with fields:
%       maxIndsAll
%       minIndsAll
%       maxWhichLeg
%       minWhichLeg
%       userSelVal
%       stepInds
%       stepWhichLeg
%   legTrack - struct of leg tracking data, output of preprocessLegTrack
%   fictrac - struct of FicTrac data, output of filtFicTrac_all
%   ephysSpikes -struct of ephys data processed for spiking, output of
%       interactEphysGetSpikes
%   
% OUTPUTS:
%   legSteps - struct of leg step data, updated with additional fields.
%     Each new field is n x 2 matrix of steps x (val start to mid), (val
%     mid to end), except stepT, which is n x 3 (t at start, mid, end)
%       stepLengths - step length in xy plane
%       stepXLengths - step length in x direction
%       stepYLengths - step length in y direction
%       stepDirections - step direction, degrees, in xy plane
%       stepDurations - step duration (in sec)
%       stepSpeeds - step speed in xy plane (body lengths/sec)
%       stepVelX - step velocity in x direction (body lengths/sec)
%       stepVelY - step velocity in y direction (body lengths/sec)
%       stepAEPX - step anterior extreme position in X
%       stepAEPY - step Y position for AEP as defined in X
%       stepPEPX - step posterior extreme position in X
%       stepPEPY - step Y position for PEP as defined in X
%       stepT - time at start, mid, and end points of step (sec)
%       stepFtFwd - FicTrac forward velocity during step (mm/sec)
%       stepFtLat - FicTrac lateral velocity during step (mm/sec)
%       stepFtYaw - FicTrac lateral velocity during step (deg/sec)
%       stepSpikeRate - spike rate during step (Hz), calculated as number
%           of spikes/step duration
%       stepSpikeRateInt - spike rate during step (Hz), calculated as
%           average spike rate during step (avg of ephysSpikes.spikeRate)
%       stepMedFiltV - median filtered voltage during step (mV), calculated
%           as average of ephysSpikes.medFiltV
%
% CREATED: 10/1/21 - HHY
%
% UPDATED:
%   10/2/21 - HHY
%   7/6/22 - HHY - add stepAEPX, stepPEPX, stepAEPY, stepPEPY and 
%       stepYLengths to legSteps output struct
%   2/12/23 - HHY - add stepSpikeRate, stepSpikeRateInt, stepMedFiltV to
%       legSteps output struct
%   3/4/23 - HHY - update so that it doesn't crash if ephysSpikes doesn't
%       exist
%   6/21/23 - HHY - update to use srnf instead of srn leg pos (outliers
%       removed)
%
function legSteps = computeStepParameters(legSteps, legTrack, fictrac, ...
    ephysSpikes)
    
    % interpolate FicTrac position to leg frame times (spline)
    interpFtFwdPos = interp1(fictrac.t, fictrac.fwdCumPos, ...
        legTrack.t, 'spline');
    interpFtYawPos = interp1(fictrac.t, fictrac.yawAngCumPos, ...
        legTrack.t, 'spline');
    interpFtLatPos = interp1(fictrac.t, fictrac.slideCumPos, ...
        legTrack.t, 'spline');

    % preallocate
    stepLengths = zeros(size(legSteps.stepInds,1), 2); 
    stepXLengths = zeros(size(legSteps.stepInds,1), 2);  
    stepYLengths = zeros(size(legSteps.stepInds,1), 2);
    stepDirections = zeros(size(legSteps.stepInds,1), 2);
    stepDurations = zeros(size(legSteps.stepInds,1), 2);
    stepSpeeds = zeros(size(legSteps.stepInds,1), 2);
    stepVelX = zeros(size(legSteps.stepInds,1), 2);
    stepVelY = zeros(size(legSteps.stepInds,1), 2);
    stepAEPX = zeros(size(legSteps.stepInds,1),2);
    stepAEPY = zeros(size(legSteps.stepInds,1),2);
    stepPEPX = zeros(size(legSteps.stepInds,1),2);
    stepPEPY = zeros(size(legSteps.stepInds,1),2);
    stepT = zeros(size(legSteps.stepInds,1),3);
    stepFtFwd = zeros(size(legSteps.stepInds,1),2);
    stepFtYaw = zeros(size(legSteps.stepInds,1),2);
    stepFtLat = zeros(size(legSteps.stepInds,1),2);
    stepSpikeRate = zeros(size(legSteps.stepInds,1),2);
    stepSpikeRateInt = zeros(size(legSteps.stepInds,1),2);
    stepMedFiltV = zeros(size(legSteps.stepInds,1),2);
    
    % loop through all steps, compute parameters for each half step
    for i = 1:size(legSteps.stepInds, 1)
        % get X position for step start
        stepStartX = legTrack.srnfLegX(legSteps.stepInds(i,1), ...
            legSteps.stepWhichLeg(i));
        % get Y position for step start
        stepStartY = legTrack.srnfLegY(legSteps.stepInds(i,1), ...
            legSteps.stepWhichLeg(i));
        % get X position for step mid
        stepMidX = legTrack.srnfLegX(legSteps.stepInds(i,2), ...
            legSteps.stepWhichLeg(i));
        % get Y position for step mid
        stepMidY = legTrack.srnfLegY(legSteps.stepInds(i,2), ...
            legSteps.stepWhichLeg(i));
        % get X position for step end
        stepEndX = legTrack.srnfLegX(legSteps.stepInds(i,3), ...
            legSteps.stepWhichLeg(i));
        % get Y position for step end
        stepEndY = legTrack.srnfLegY(legSteps.stepInds(i,3), ...
            legSteps.stepWhichLeg(i));
        
        

        % get length of first half step (start to mid)
        stepLengths(i,1) = distBtw2Pts(stepStartX, stepStartY, ...
            stepMidX, stepMidY);
        % get length of second half step (mid to end)
        stepLengths(i,2) = distBtw2Pts(stepMidX, stepMidY, stepEndX, ...
            stepEndY);
        
        % step X lengths
        stepXLengths(i,1) = stepMidX - stepStartX;
        stepXLengths(i,2) = stepEndX - stepMidX;

        % step Y lengths
        stepYLengths(i,1) = stepMidY - stepStartY;
        stepYLengths(i,2) = stepEndY - stepMidY;

        % step AEP
        % get whether start or mid position is most anterior (min position)
        [~, aepInd] = min([stepStartX, stepMidX]);
        % match AEP X to AEP Y; PEP is one that's not AEP
        if (aepInd == 1)
            stepAEPX(i,1) = stepStartX;
            stepAEPY(i,1) = stepStartY;
            stepPEPX(i,1) = stepMidX;
            stepPEPY(i,1) = stepMidY;
        else
            stepAEPX(i,1) = stepMidX;
            stepAEPY(i,1) = stepMidY;
            stepPEPX(i,1) = stepStartX;
            stepPEPY(i,1) = stepStartY;
        end
        % same for second half step
        [~, aepInd] = min([stepMidX, stepEndX]);
        if (aepInd == 1)
            stepAEPX(i,2) = stepMidX;
            stepAEPY(i,2) = stepMidY;
            stepPEPX(i,2) = stepEndX;
            stepPEPY(i,2) = stepEndY;
        else
            stepAEPX(i,2) = stepEndY;
            stepAEPY(i,2) = stepEndY;
            stepPEPX(i,2) = stepMidX;
            stepPEPY(i,2) = stepEndY;
        end

        % get direction of first half step (start to mid)
        stepDirections(i,1) = findAngle2Pts(stepStartX, stepStartY, ...
            stepMidX, stepMidY);
        % get direction of second half step (mid to end)
        stepDirections(i,2) = findAngle2Pts(stepMidX, stepMidY, ...
            stepEndX, stepEndY);

        % get time for step start
        stepStartT = legTrack.t(legSteps.stepInds(i,1));
        % get time for step mid
        stepMidT = legTrack.t(legSteps.stepInds(i,2));
        % get time for step end
        stepEndT = legTrack.t(legSteps.stepInds(i,3));

        % step time
        stepT(i,:) = [stepStartT, stepMidT, stepEndT];

        % get duration of first half step (start to mid)
        stepDurations(i,1) = stepMidT - stepStartT;
        % get duration of second half step (mid to end)
        stepDurations(i,2) = stepEndT - stepMidT;

        % get velocity in xy plane of first half step (start to mid)
        stepSpeeds(i,1) = stepLengths(i,1) / stepDurations(i,1);
        % get velocity in xy plane of second half step (mid to end)
        stepSpeeds(i,2) = stepLengths(i,2) / stepDurations(i,2);

        % get velocity in x for first half of step (start to mid)
        stepVelX(i,1) = (stepMidX - stepStartX) / stepDurations(i,1);
        % get velocity in x for second half step (mid to end)
        stepVelX(i,2) = (stepEndX - stepMidX) / stepDurations(i,2);

        % get velocity in y for first half of step (start to mid)
        stepVelY(i,1) = (stepMidY - stepStartY) / stepDurations(i,1);
        % get velocity in y for second half step (mid to end)
        stepVelY(i,2) = (stepEndY - stepMidY) / stepDurations(i,2); 
        
        % FicTrac fwd position for step start point
        stepStartFtFwdPos = interpFtFwdPos(legSteps.stepInds(i,1));
        % FicTrac fwd position for step mid point
        stepMidFtFwdPos = interpFtFwdPos(legSteps.stepInds(i,2));
        % FicTrac fwd position for step end point
        stepEndFtFwdPos = interpFtFwdPos(legSteps.stepInds(i,3));

        % FicTrac fwd velocity, first half step
        stepFtFwd(i,1) = (stepMidFtFwdPos-stepStartFtFwdPos) / ...
            stepDurations(i,1);
        % second half step
        stepFtFwd(i,2) = (stepEndFtFwdPos - stepMidFtFwdPos) / ...
            stepDurations(i,2);

        % FicTrac yaw position for step start point
        stepStartFtYawPos = interpFtYawPos(legSteps.stepInds(i,1));
        % FicTrac yaw position for step mid point
        stepMidFtYawPos = interpFtYawPos(legSteps.stepInds(i,2));
        % FicTrac yaw position for step end point
        stepEndFtYawPos = interpFtYawPos(legSteps.stepInds(i,3));

        % FicTrac yaw velocity, first half step
        stepFtYaw(i,1) = (stepMidFtYawPos-stepStartFtYawPos) / ...
            stepDurations(i,1);
        % second half step
        stepFtYaw(i,2) = (stepEndFtYawPos - stepMidFtYawPos) / ...
            stepDurations(i,2);

        % FicTrac lat position for step start point
        stepStartFtLatPos = interpFtLatPos(legSteps.stepInds(i,1));
        % FicTrac lat position for step mid point
        stepMidFtLatPos = interpFtLatPos(legSteps.stepInds(i,2));
        % FicTrac lat position for step end point
        stepEndFtLatPos = interpFtLatPos(legSteps.stepInds(i,3));

        % FicTrac lat velocity, first half step
        stepFtLat(i,1) = (stepMidFtLatPos-stepStartFtLatPos) / ...
            stepDurations(i,1);
        % second half step
        stepFtLat(i,2) = (stepEndFtLatPos - stepMidFtLatPos) / ...
            stepDurations(i,2);

        % process ephys info for leg steps if present
        if ~isempty(ephysSpikes)
            % indices into ephysSpikes vectors for time points matching start,
            %  mid, end points of step
            stepStartEphysInd = find(ephysSpikes.t >= stepStartT, 1, 'first');
            stepMidEphysInd = find(ephysSpikes.t >= stepMidT, 1, 'first');
            stepEndEphysInd = find(ephysSpikes.t >= stepEndT, 1, 'first');
    
            % number of spikes during first half step
            numSpikes1 = length(find((ephysSpikes.startInd > stepStartEphysInd) & ...
                (ephysSpikes.startInd <= stepMidEphysInd)));
            % second half step
            numSpikes2 = length(find((ephysSpikes.startInd > stepMidEphysInd) & ...
                (ephysSpikes.startInd <= stepEndEphysInd)));
    
            % convert number of spikes to spike rate (divide by step duration)
            stepSpikeRate(i,1) = numSpikes1 / stepDurations(i,1);
            stepSpikeRate(i,2) = numSpikes2 / stepDurations(i,2);
    
            % mean spike rate during first half step
            stepSpikeRateInt(i,1) = ...
                mean(ephysSpikes.spikeRate(stepStartEphysInd:stepMidEphysInd));
            % mean spike rate during second half step
            stepSpikeRateInt(i,2) = ...
                mean(ephysSpikes.spikeRate(stepMidEphysInd:stepEndEphysInd));
    
            % mean filtered voltage during first half step
            stepMedFiltV(i,1) = ...
                mean(ephysSpikes.medFiltV(stepStartEphysInd:stepMidEphysInd));
            % mean spike rate during second half step
            stepMedFiltV(i,2) = ...
                mean(ephysSpikes.medFiltV(stepMidEphysInd:stepEndEphysInd));

        end

    end
    
    % add these parameters to legSteps struct
    legSteps.stepLengths = stepLengths;
    legSteps.stepXLengths = stepXLengths;
    legSteps.stepYLengths = stepYLengths;
    legSteps.stepDirections = stepDirections;
    legSteps.stepDurations = stepDurations;
    legSteps.stepSpeeds = stepSpeeds;
    legSteps.stepVelX = stepVelX;
    legSteps.stepVelY = stepVelY;
    legSteps.stepAEPX = stepAEPX;
    legSteps.stepAEPY = stepAEPY;
    legSteps.stepPEPX = stepPEPX;
    legSteps.stepPEPY = stepPEPY;
    legSteps.stepT = stepT;
    legSteps.stepFtFwd = stepFtFwd;
    legSteps.stepFtLat = stepFtLat;
    legSteps.stepFtYaw = stepFtYaw;
    % only when there's ephys data
    if ~isempty(ephysSpikes)
        legSteps.stepSpikeRate = stepSpikeRate;
        legSteps.stepSpikeRateInt = stepSpikeRateInt;
        legSteps.stepMedFiltV = stepMedFiltV;
    end
end