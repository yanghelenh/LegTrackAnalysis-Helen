% determineLegPhase.m
%
% Function to determine phase of each leg, in degrees from -180 to 180.
% Uses Hilbert transform.
%
% INPUTS:
%   legPos - matrix of leg positions, num frames x num tracked pts
%   legInd - indicies of leg tracked pts
%
% OUTPUTS:
%   legPhase - matrix of leg phases, from negative -180 to 180; num frames x
%       num tracked pts
%
% CREATED: 11/19/20 - HHY
%
% UPDATED:
%   11/19/20 - HHY
%
function legPhase = determineLegPhase(legPos, legInd)
    % preallocate
    legPhase = zeros(size(legPos,1), length(legInd));
    
    % get leg phase, in degrees
    for i = 1:length(legInd)
        legPhase(:,legInd(i)) = angle(hilbert(legPos(:,legInd(i)))) ...
            .* (180/pi);
    end
end