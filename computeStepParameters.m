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
%   
% OUTPUTS:
%   legSteps - struct of leg step data, updated with additional fields.
%     Each new field is n x 2 matrix of steps x (val start to mid), (val
%     mid to end)
%       stepLengths - step length in xy plane
%       stepDirections - step direction, degrees, in xy plane
%       stepDurations - step duration (in sec)
%       stepSpeeds - step speed in xy plane (body lengths/sec)
%       stepVelX - step velocity in x direction (body lengths/sec)
%       stepVelY - step velocity in y direction (body lengths/sec)
%
% CREATED: 10/1/21 - HHY
%
% UPDATED:
%   10/1/21 - HHY
%
function legSteps = computeStepParameters(legSteps, legTrack)
    
    % preallocate
    stepLengths = zeros(size(legSteps.stepInds,1), 2); 
    stepDirections = zeros(size(legSteps.stepInds,1), 2);
    stepDurations = zeros(size(legSteps.stepInds,1), 2);
    stepSpeeds = zeros(size(legSteps.stepInds,1), 2);
    stepVelX = zeros(size(legSteps.stepInds,1), 2);
    stepVelY = zeros(size(legSteps.stepInds,1), 2);
    
    % loop through all steps, compute parameters for each half step
    for i = 1:size(legSteps.stepInds, 1)
        % get X position for step start
        stepStartX = legTrack.srnLegX(legSteps.stepInds(i,1), ...
            legSteps.stepWhichLeg(i));
        % get Y position for step start
        stepStartY = legTrack.srnLegY(legSteps.stepInds(i,1), ...
            legSteps.stepWhichLeg(i));
        % get X position for step mid
        stepMidX = legTrack.srnLegX(legSteps.stepInds(i,2), ...
            legSteps.stepWhichLeg(i));
        % get Y position for step mid
        stepMidY = legTrack.srnLegY(legSteps.stepInds(i,2), ...
            legSteps.stepWhichLeg(i));
        % get X position for step end
        stepEndX = legTrack.srnLegX(legSteps.stepInds(i,3), ...
            legSteps.stepWhichLeg(i));
        % get Y position for step end
        stepEndY = legTrack.srnLegY(legSteps.stepInds(i,3), ...
            legSteps.stepWhichLeg(i));

        % get length of first half step (start to mid)
        stepLengths(i,1) = distBtw2Pts(stepStartX, stepStartY, ...
            stepMidX, stepMidY);
        % get length of second half step (mid to end)
        stepLengths(i,2) = distBtw2Pts(stepMidX, stepMidY, stepEndX, ...
            stepEndY);

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
    end
    
    % add these parameters to legSteps struct
    legSteps.stepLengths = stepLengths;
    legSteps.stepDirections = stepDirections;
    legSteps.stepDurations = stepDurations;
    legSteps.stepSpeeds = stepSpeeds;
    legSteps.stepVelX = stepVelX;
    legSteps.stepVelY = stepVelY;
end