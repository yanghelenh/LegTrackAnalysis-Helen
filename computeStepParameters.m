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
%     mid to end), except stepPosX and stepPosY, which are n x 3 (val
%     start, mid, end)
%       stepLengths - step length in xy plane
%       stepXLengths - step length in x direction
%       stepYLengths - step length in y direction
%       stepDirections - step direction, degrees, in xy plane
%       stepDurations - step duration (in sec)
%       stepSpeeds - step speed in xy plane (body lengths/sec)
%       stepVelX - step velocity in x direction (body lengths/sec)
%       stepVelY - step velocity in y direction (body lengths/sec)
%       stepPosX - step position in x axis
%       stepPosY - step position in y axis
%       stepFtFwd - FicTrac forward velocity during step (mm/sec)
%       stepFtLat - FicTrac lateral velocity during step (mm/sec)
%       stepFtYaw - FicTrac lateral velocity during step (deg/sec)
%
% CREATED: 10/1/21 - HHY
%
% UPDATED:
%   10/2/21 - HHY
%   7/6/22 - HHY - add stepPosX, stepPosY, and stepYLengths to legSteps 
%       output struct
%
function legSteps = computeStepParameters(legSteps, legTrack, fictrac)
    
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
    stepPosX = zeros(size(legSteps.stepInds,1),3);
    stepPosY = zeros(size(legSteps.stepInds,1),3);
    stepFtFwd = zeros(size(legSteps.stepInds,1),2);
    stepFtYaw = zeros(size(legSteps.stepInds,1),2);
    stepFtLat = zeros(size(legSteps.stepInds,1),2);
    
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
        
        % step X lengths
        stepXLengths(i,1) = stepMidX - stepStartX;
        stepXLengths(i,2) = stepEndX - stepMidX;

        % step Y lengths
        stepYLengths(i,1) = stepMidY - stepStartY;
        stepYLengths(i,2) = stepEndY - stepMidY;

        % step X and Y positions
        stepPosX(i,:) = [stepStartX, stepMidX, stepEndX];
        stepPosY(i,:) = [stepStartY, stepMidY, stepEndY];

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
    legSteps.stepPosX = stepPosX;
    legSteps.stepPosY = stepPosY;
    legSteps.stepFtFwd = stepFtFwd;
    legSteps.stepFtLat = stepFtLat;
    legSteps.stepFtYaw = stepFtYaw;
end