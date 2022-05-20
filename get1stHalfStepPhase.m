% get1stHalfStepPhase.m
%
% Helper function for getLegPhaseFromSteps() that takes in step start max
%  position and step mid min position to convert a vector of leg positions
%  between the two to leg phases between 0 and 180 degrees, with a linear
%  mapping.
% For a linear mapping: phase = 180(pos - startPos)/(midPos - startPos)
%
% INPUTS:
%   stepStartPos - scalar value for position of step start
%   stepMidPos - scalar value for position of step mid
%   legPos - vector of leg positions to convert to leg phases
%
% OUTPUTS:
%   legPhase - vector of leg phases (between 0 and 180), same size as legPos
%
% CREATED: 10/21/21 - HHY
%
% UPDATED:
%   10/21/21 - HHY
%
function legPhase = get1stHalfStepPhase(stepStartPos, stepMidPos, legPos)

    % phase val for midpoint of step; here, 180 for degrees, but could be
    %  radians
    MIDSTEP_PHASE_VAL = 180; 

    legPhase = (MIDSTEP_PHASE_VAL * (legPos - stepStartPos)) / ...
        (stepMidPos - stepStartPos);
end