% interactGetNotMovingInd.m
%
% Function for selecting not-moving times. 
% Brings up GUI with plot of mid-leg positions, shaded regions of not 
%  moving & red dots on not-moving, and sliders to adjust all parameters 
%  that determine not-moving portions.
%
% INPUTS:
%   legTrack - struct of leg tracking data, output of preprocessLegTrack
%   notMoveParams - struct of initial values of parameters
%       medFiltNumSamps - number of samples for median filtering leg
%           velocity
%       zeroXVelThreshMed - zero velocity threshold (on speed, actually),
%           for median filtered leg X velocity
%       zeroYVelThreshMed - zero velocity threshold (on speed, actually),
%           for median filtered leg Y velocity
%       zeroXVelThresh - zero velocity threshold on leg X velocity, w/o
%           filtering
%       zeroYVelThresh - zero velocity threshold on leg Y velocity, w/o
%           filtering
%       movePosXVelThresh - threshold for movement, positive X velocity
%       moveNegXVelThresh - threshold for movement, negative X velocity
%       movePosYVelThresh - threshold for movement, positive Y velocity
%       moveNegYVelThresh - threshold for movement, negative Y velocity
%       minBoutLen - minimum length of not-moving bout, in samples
%       stepNegXVelThresh - threshold for step, negative X velocity
%       maxTimeFromStep - maximum number of samples from step, for midpoint
%           of not-moving bout
%       adjBoutSepInit - if not-moving bouts are less than this many 
%           samples apart initially, merge them
%       adjBoutSepEnd - if not-moving bouts are less than this many 
%           samples apart at the end of the analysis, merge them
%   r2LegInd - index of right mid-leg
%   l2LegInd - index of left mid-leg
%
% OUTPUTS:
%   notMoveInd - indicies for when fly not moving
%   notMoveBout - n x 2 matrix of not move bout start (col 1) and end
%       (col 2) indicies
%   moveInd - indicies for when fly moving (inverse of notMoveInd)
%   moveBout - n x 2 matrix of move bout start (col 1) and end (col 2)
%       indicies
%   notMoveParams - struct of final values of parameters
%
% CREATED: 9/24/21 - HHY
%
% UPDATED:
%   9/24/21 - HHY
%
function [notMoveInd, notMoveBout, moveInd, moveBout, notMoveParams] = ...
    interactGetNotMovingInd(legTrack, notMoveParams, r2LegInd, l2LegInd)

end