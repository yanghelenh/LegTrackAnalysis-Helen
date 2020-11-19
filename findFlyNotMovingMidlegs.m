% findFlyNotMovingMidlegs.m
%
% Function to find when the fly isn't moving, based on the criteria that
%  the 2 midlegs are not moving in X (velocity below threshold).
% Finds based on smoothed velocity below threshold and non-smoothed
%  velocity also below (separate) threshold
% Returns indicies of not moving times
%
% INPUTS:
%   srnLegX - matrix of leg X coordinates after smoothing, rotating, and
%       normalization (output of shiftRotateNormalizeLegPos)
%   srnLegY - matrix of leg Y coordinates after smoothing, rotating, and
%       normalization (output of shiftRotateNormalizeLegPos)
%   notMoveParams - struct of parameters
%       padLen - pad length, for Gaussian process smoothing of leg position
%       sigma - sigma for Gaussian process smoothing of leg position
%       medFiltNumSamps - number of samples for median filtering leg
%           velocity
%       zeroVelThreshMed - zero velocity threshold (on speed, actually),
%           for median filtered leg velocity
%       zeroVelThresh - zero velocity threshold on leg velocity, w/o
%           filtering
%       movePosVelThresh - threshold for movement, positive velocity
%       moveNegVelThresh - threshold for movement, negative velocity
%       minBoutLen - minimum length of not-moving bout, in samples
%       r2LegInd - index of right mid-leg
%       l2LegInd - index of left mid-leg
%
% OUTPUTS:
%   zeroVelInd - indices for when fly isn't moving
%
% CREATED: 11/16/20 - HHY
%
% UPDATED:
%   11/16/20 - HHY
%
function zeroVelInd = findFlyNotMovingMidlegs(srnLegX, srnLegY, ...
    notMoveParams)

    % preallocate
    srnLegXSmo = zeros(size(srnLegX));

    % loop through all tracked points
    for i = 1:size(srnLegXSmo,2)
        smoPosPy = py.proc_utils.safe_interp_conv_smooth_ball_data(...
            srnLegX(:,i)', notMoveParams.padLen, notMoveParams.sigma);
        % convert from python to matlab data format
        smoPos = cell2mat(cell(smoPosPy.tolist()));
        srnLegXSmo(:,i) = smoPos';
    end

    % get leg x velocities by taking gradient
    legXVel = zeros(size(srnLegXSmo));
    for i = 1:size(legXVel,2)
        legXVel(:,i) = gradient(srnLegXSmo(:,i));
    end

    % median filter leg velocities, mid legs
    medFiltR2 = medfilt1(legXVel(:,notMoveParams.r2LegInd),...
        notMoveParams.medFiltNumSamps);
    medFiltL2 = medfilt1(legXVel(:,notMoveParams.l2LegInd),...
        notMoveParams.medFiltNumSamps);

    % not moving as speed less than threshold
    r2ZeroVelInd = find((abs(medFiltR2) < ...
        notMoveParams.zeroVelThreshMed) & ...
        (abs(legXVel(:,notMoveParams.r2LegInd)) < ...
        notMoveParams.zeroVelThresh));

    % find inverse of R2 zeroVelInd -> moving index
    allInd = (1:length(legXVel(:,notMoveParams.r2LegInd)))';
    r2MoveInd = setdiff(allInd, r2ZeroVelInd); 

    % find moving bouts
    [r2MoveBoutStarts, r2MoveBoutEnds, ~] = findBouts(r2MoveInd);

    % loop through all moving bouts
    for i = 1:length(r2MoveBoutStarts)
        boutInd = (r2MoveBoutStarts(i):r2MoveBoutEnds(i))';
        boutVels = legXVel(boutInd,2);
        % if movement bout doesn't exceed threshold, append to zeroVelInd
        if ~(sum(boutVels >= notMoveParams.movePosVelThresh) || ...
                sum(boutVels <= notMoveParams.moveNegVelThresh))
            r2ZeroVelInd = [r2ZeroVelInd; boutInd];
        end
    end

    % resort r2ZeroVelInd in ascending order
    r2ZeroVelInd = sort(r2ZeroVelInd);

    % same for left mid leg
    l2ZeroVelInd = find((abs(medFiltL2) < ...
        notMoveParams.zeroVelThreshMed) & ...
        (abs(legXVel(:,notMoveParams.l2LegInd)) < zeroVelThresh));

    l2MoveInd = setdiff(allInd, l2ZeroVelInd);

    [l2MoveBoutStarts, l2MoveBoutEnds, ~] = findBouts(l2MoveInd);

    for i = 1:length(l2MoveBoutStarts)
        boutInd = (l2MoveBoutStarts(i):l2MoveBoutEnds(i))';
        boutVels = legXVel(boutInd,notMoveParams.l2LegInd);
        % if movement bout doesn't exceed threshold, append to zeroVelInd
        if ~(sum(boutVels >= notMoveParams.movePosVelThresh) || ...
                sum(boutVels <= notMoveParams.moveNegVelThresh))
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