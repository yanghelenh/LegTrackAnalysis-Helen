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

end