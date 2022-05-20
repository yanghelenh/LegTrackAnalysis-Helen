% get2ndHalfStepPhase.m
%
% Helper function for getLegPhaseFromSteps() that takes in step mid min
%  position and step end max position to convert a vector of leg positions
%  between the two to leg phases between 180 and 360 degrees, with a linear
%  mapping.
% For a linear mapping: 
%    phase = 180(pos + endPos - 2(midPos))/(endPos - midPos)
%
% INPUTS:
%   stepMidPos - scalar value for position of step mid
%   stepEndPos - scalar value for position of step end
%   legPos - vector of leg positions to convert to leg phases
%
% OUTPUTS:
%   legPhase - vector of leg phases (between 180 and 360), same size as 
%     legPos
%
% CREATED: 10/21/21 - HHY
%
% UPDATED:
%   10/21/21 - HHY
%
function legPhase = get2ndHalfStepPhase(stepMidPos, stepEndPos, legPos)

    % phase val for midpoint of step; here, 180 for degrees, but could be
    %  radians
    MIDSTEP_PHASE_VAL = 180; 

    legPhase = (MIDSTEP_PHASE_VAL*(legPos+stepEndPos-(2*stepMidPos))) / ...
        (stepEndPos - stepMidPos);
    
end