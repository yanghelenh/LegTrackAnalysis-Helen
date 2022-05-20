% getLegPhaseDiffs.m
%
% Function that takes in leg phase (computed in getLegPhaseFromSteps()) and
%  returns pairwise differences in phase between legs
% Leg pairs: 
%   R1-L1, R2-L2, R3-L3, R1-R2, R2-R3, R1-R3, L1-L2, L2-L3, L1-L3, L1-R2,
%     R2-L3, R1-L2, L2-R3
% Leg phase can be NaN so phase difference will be NaN if at least one of
%   the pair is NaN
% 
% INPUTS:
%   legPhase - numFrames x 6 matrix of leg phases at each frame for all 6
%       legs
%   legIDs - struct of leg info
%       ind - indices of legs
%       names - names of legs
%   units - 'radians' or 'degrees' for units of leg phase
%
% OUTPUTS:
%   phaseDiffs - struct of phase differences between legs, with fields
%       left/right pair, same leg: r1l1, r2l2, r3l3
%       pairs of legs on same side: r1r2, r2r3, r1r3, l1l2, l2l3, l1l3
%       pairs of legs in same tripod: l1r2, r2l3, r1l2, l2r3
%
% CREATED: 10/26/21 - HHY
%
% UPDATED:
%   10/26/21 - HHY
%
function phaseDiffs = getLegPhaseDiffs(legPhase, legIDs, units)

    % phase difference wrapped to 360 or 2 pi depending on units
    if (strcmpi(units, 'radians'))
        wrapFn = @wrapTo2Pi;
    elseif (strcmpi(units, 'degrees'))
        wrapFn = @wrapTo360;
    else
        disp('units must be radians or degrees');
        phaseDiffs = [];
        return;
    end
    
    % get indices for each leg
    r1Ind = legIDs.ind(strcmpi(legIDs.names, 'r1'));
    r2Ind = legIDs.ind(strcmpi(legIDs.names, 'r2'));
    r3Ind = legIDs.ind(strcmpi(legIDs.names, 'r3'));
    l1Ind = legIDs.ind(strcmpi(legIDs.names, 'l1'));
    l2Ind = legIDs.ind(strcmpi(legIDs.names, 'l2'));
    l3Ind = legIDs.ind(strcmpi(legIDs.names, 'l3'));
    
    % get phase differences
    phaseDiffs.r1l1 = wrapFn(legPhase(:,r1Ind) - legPhase(:,l1Ind));
    phaseDiffs.r2l2 = wrapFn(legPhase(:,r2Ind) - legPhase(:,l2Ind));
    phaseDiffs.r3l3 = wrapFn(legPhase(:,r3Ind) - legPhase(:,l3Ind));
    
    phaseDiffs.r1r2 = wrapFn(legPhase(:,r1Ind) - legPhase(:,r2Ind));
    phaseDiffs.r2r3 = wrapFn(legPhase(:,r2Ind) - legPhase(:,r3Ind));
    phaseDiffs.r1r3 = wrapFn(legPhase(:,r1Ind) - legPhase(:,r3Ind));
    
    phaseDiffs.l1l2 = wrapFn(legPhase(:,l1Ind) - legPhase(:,l2Ind));
    phaseDiffs.l2l3 = wrapFn(legPhase(:,l2Ind) - legPhase(:,l3Ind));
    phaseDiffs.l1l3 = wrapFn(legPhase(:,l1Ind) - legPhase(:,l3Ind));
    
    phaseDiffs.l1r2 = wrapFn(legPhase(:,l1Ind) - legPhase(:,r2Ind));
    phaseDiffs.r2l3 = wrapFn(legPhase(:,r2Ind) - legPhase(:,l3Ind));
    phaseDiffs.r1l2 = wrapFn(legPhase(:,r1Ind) - legPhase(:,l2Ind));
    phaseDiffs.l2r3 = wrapFn(legPhase(:,l2Ind) - legPhase(:,r3Ind));

    
end