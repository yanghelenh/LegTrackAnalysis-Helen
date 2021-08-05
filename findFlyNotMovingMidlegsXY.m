% findFlyNotMovingMidlegsXY.m
%
% Function to find when the fly isn't moving, based on the criteria that
%  the 2 midlegs are not moving in X and Y (velocity below threshold).
% Finds based on smoothed velocity below threshold and non-smoothed
%  velocity also below (separate) threshold
% Returns indicies of not moving times
%
% INPUTS:
%   xVel - matrix of leg velocities in X direction
%   yVel - matrix of leg velocities in Y direction
%   notMoveParams - struct of parameters
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
%       r2LegInd - index of right mid-leg
%       l2LegInd - index of left mid-leg
%
% OUTPUTS:
%   zeroVelInd - indices for when fly isn't moving
%
% CREATED: 12/7/20 - HHY
%
% UPDATED:
%   12/7/20 - HHY
%
function zeroVelInd = findFlyNotMovingMidlegsXY(xVel, yVel, notMoveParams)

    % median filter leg velocities, mid legs
    medFiltXR2 = medfilt1(xVel(:,notMoveParams.r2LegInd),...
        notMoveParams.medFiltNumSamps);
    medFiltXL2 = medfilt1(xVel(:,notMoveParams.l2LegInd),...
        notMoveParams.medFiltNumSamps);
    medFiltYR2 = medfilt1(yVel(:,notMoveParams.r2LegInd),...
        notMoveParams.medFiltNumSamps);
    medFiltYL2 = medfilt1(yVel(:,notMoveParams.l2LegInd),...
        notMoveParams.medFiltNumSamps);

    % not moving as speed less than threshold, for x and y velocities
    r2ZeroXVelInd = find((abs(medFiltXR2) < ...
        notMoveParams.zeroXVelThreshMed) & ...
        (abs(xVel(:,notMoveParams.r2LegInd)) < ...
        notMoveParams.zeroXVelThresh));
    r2ZeroYVelInd = find((abs(medFiltYR2) < ...
        notMoveParams.zeroYVelThreshMed) & ...
        (abs(yVel(:,notMoveParams.r2LegInd)) < ...
        notMoveParams.zeroYVelThresh));
    % not moving as intersection between both x and y
    r2ZeroVelInd = intersect(r2ZeroXVelInd, r2ZeroYVelInd);

    % find inverse of R2 zeroVelInd -> moving index
    allInd = (1:length(xVel(:,notMoveParams.r2LegInd)))';
    r2MoveInd = setdiff(allInd, r2ZeroVelInd); 

    % find moving bouts
    [r2MoveBoutStarts, r2MoveBoutEnds, ~] = findBouts(r2MoveInd);

    % loop through all moving bouts
    for i = 1:length(r2MoveBoutStarts)
        boutInd = (r2MoveBoutStarts(i):r2MoveBoutEnds(i))';
        boutXVels = xVel(boutInd,notMoveParams.r2LegInd);
        boutYVels = yVel(boutInd,notMoveParams.r2LegInd);
        % if movement bout doesn't exceed threshold, append to zeroVelInd;
        %  has to fail to exceed thresh in both x and y
        if ~(sum(boutXVels >= notMoveParams.movePosXVelThresh) || ...
                sum(boutXVels <= notMoveParams.moveNegXVelThresh)) && ...
            ~(sum(boutYVels >= notMoveParams.movePosYVelThresh) || ...
                sum(boutYVels <= notMoveParams.moveNegYVelThresh))
            r2ZeroVelInd = [r2ZeroVelInd; boutInd];
        end
    end

    % resort r2ZeroVelInd in ascending order
    r2ZeroVelInd = sort(r2ZeroVelInd);

    % same for left mid leg
    l2ZeroXVelInd = find((abs(medFiltXL2) < ...
        notMoveParams.zeroXVelThreshMed) & ...
        (abs(xVel(:,notMoveParams.l2LegInd)) < ...
        notMoveParams.zeroXVelThresh));
    l2ZeroYVelInd = find((abs(medFiltYL2) < ...
        notMoveParams.zeroYVelThreshMed) & ...
        (abs(yVel(:,notMoveParams.l2LegInd)) < ...
        notMoveParams.zeroYVelThresh));
    l2ZeroVelInd = intersect(l2ZeroXVelInd, l2ZeroYVelInd);

    l2MoveInd = setdiff(allInd, l2ZeroVelInd);

    [l2MoveBoutStarts, l2MoveBoutEnds, ~] = findBouts(l2MoveInd);

    for i = 1:length(l2MoveBoutStarts)
        boutInd = (l2MoveBoutStarts(i):l2MoveBoutEnds(i))';
        boutXVels = xVel(boutInd,notMoveParams.l2LegInd);
        boutYVels = yVel(boutInd,notMoveParams.l2LegInd);
        % if movement bout doesn't exceed threshold, append to zeroVelInd
        if ~(sum(boutXVels >= notMoveParams.movePosXVelThresh) || ...
                sum(boutXVels <= notMoveParams.moveNegXVelThresh)) && ...
           ~(sum(boutYVels >= notMoveParams.movePosYVelThresh) || ...
                sum(boutYVels <= notMoveParams.moveNegYVelThresh))
            l2ZeroVelInd = [l2ZeroVelInd; boutInd];
        end
    end

    l2ZeroVelInd = sort(l2ZeroVelInd);

    % fly not moving when both right and left mid legs not moving (fly can
    %  groom with front or back legs and be stationary, but they don't groom
    %  with midlegs)
    zeroVelInd = intersect(r2ZeroVelInd, l2ZeroVelInd);

    % each not moving bout must be of minimum length; if not, remove
    % find start, end, and duration of zero velocity bouts
    [zeroVelBoutStarts, zeroVelBoutEnds, zeroVelBoutDur] = ...
        findBouts(zeroVelInd);

    % find bouts that are less than minimum long
    zeroVelShortBouts = find(zeroVelBoutDur < notMoveParams.minBoutLen);

    % remove these bouts from zeroVelInd - set all values within zeroVelInd
    %  that fall into these short bouts to NaN
    % loop through all short bouts, check if zeroVelInd falls within start and
    %  end
    for i = 1:length(zeroVelShortBouts)
        startInd = zeroVelBoutStarts(zeroVelShortBouts(i));
        endInd = zeroVelBoutEnds(zeroVelShortBouts(i));

        zeroVelInd((zeroVelInd >= startInd)&(zeroVelInd <= endInd)) = NaN;
    end
    % remove NaNs from zeroVelInd
    zeroVelInd(isnan(zeroVelInd)) = [];
end