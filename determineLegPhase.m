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
%   9/14/21 - HHY - update how radian to degree conversion works
%   9/20/21 - HHY - z-score leg position for Hilbert transform
%
function legPhase = determineLegPhase(legPos, legInd)
    % preallocate
    legPhase = zeros(size(legPos,1), length(legInd));
    
    
    % loop through all legs
    for i = 1:length(legInd)
        
        % z-score leg position (otherwise, Hilbert transform doesn't work)
        zscLegPos = zscore(legPos(:,legInd(i)));
         
        % get leg phase, in degrees
        legPhase(:,legInd(i)) = rad2deg(angle(hilbert(zscLegPos)));
    end
end