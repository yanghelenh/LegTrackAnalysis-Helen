% callSwingStanceSteps.m
%
% Function that takes in steps (defined by start, mid, end indices),
%  step parameters (calculated in computeStepParameters()), and FicTrac
%  info to call each step's half step as swing or stance
% Option to call swing/stance based on half step duration: shorter duration
%  of two is swing, other is stance
% Option to call swing/stance based on FicTrac velocity: same direction is
%  stance, opposite direction is swing
% Select one of two options
%
% INPUTS:
%   legSteps - struct of leg step data, output from processLegTrack(), with
%       fields:
%       maxIndsAll
%       minIndsAll
%       maxWhichLeg
%       minWhichLeg
%       userSelVal
%       stepInds
%       stepWhichLeg
%       stepLengths
%       stepDirections
%       stepDurations
%       stepSpeeds
%       stepVelX
%       stepVelY
%   legTrack - struct of leg tracking data, output of preprocessLegTrack
%   whichMethod - string for which method to determine swing/stance:
%       'duration' or 'fictrac'
%   fictrac - struct of fictrac data
%
% OUTPUTS:
%   legSteps - struct of leg step data, updated with additional fields:
%       stepSwingStance - n x 2 matrix for swing/stance calls for each half
%           step; 1 for stance, -1 for swing
%       swingStanceMethod - whichMethod from inputs
%
% CREATED: 10/1/21 - HHY
%
% UPDATED:
%   10/1/21 - HHY
%
function legSteps = callSwingStanceSteps(legSteps, legTrack, ...
    whichMethod, fictrac)

end