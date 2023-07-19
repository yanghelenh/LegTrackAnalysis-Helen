% findYawVelPeaksFT.m
%
% Helper function for saveBallLegStepParamCond_bouts() that takes in
%  fictracSmo, cond, and moveNotMove and returns times for valid yaw
%  velocity peaks as well as start/end times for corresponding bouts.
% Boolean for whether to extract right or left turns
% 
% Adapted from findCondYawVelPeaks(), which operates on free-walking data.
% Use times instead of indices b/c legStep indices and FicTrac indices
%  aren't the same, unlike in free-walking data. Also, indices just in case
%
% INPUTS:
%   fictracSmo - pData output struct, from computeSmoFictrac()
%   cond - pass from saveBallLegStepParamCond_bouts()
%   moveNotMove - pData output struct
%   rightTurn - boolean for whether to extract right turns (false = left)
%
% OUTPUTS:
%   yawVelPeakTimes - vector of times corresponding to yaw
%       velocity peaks (that meet cond criteria), length n
%   boutStartTimes - vector of times corresponding to start of bouts,
%       length n
%   boutEndTimes - vector of times corresponding to end of bouts, length n
%   yawVelPeakInd - vector of indices of fictracSmo/Proc corresponding to
%       yaw velocity peaks (that meet cond criteria), length n
%   boutStartInd - vector of indices corresponding to start of bouts,
%       length n
%   boutEndInd - vector of indices corresponding to end of bouts, length n
%
% CREATED: 7/17/23 - HHY
%
% UPDATED:
%   7/17/23- HHY
%
function [yawVelPeakTimes, boutStartTimes, boutEndTimes, ...
    yawVelPeakInd, boutStartInd, boutEndInd] = ...
    findYawVelPeaksFT(fictracSmo, minYawThresh, turnDur, moveNotMove, ...
    rightTurn)

    % compensate for bias in fly's yaw and slide velocities
%     fictracSmo.yawAngVel = fictracSmo.yawAngVel - fictracSmo.angVelBias;
%     fictracSmo.slideVel = fictracSmo.slideVel - fictracSmo.slideVelBias;

    % check if we're extracting right or left turns
    if (rightTurn)
        angVelSmoS = fictracSmo.yawAngVel;
    else
        angVelSmoS = -1 * fictracSmo.yawAngVel;
    end
        
    % find all yaw velocity peaks, no conditioning yet
    % always operates on fictracSmo.yawAngVel
    [~, pkInds] = findpeaks(angVelSmoS);

    % remove peaks during not moving times
    pkInds = setdiff(pkInds, moveNotMove.ftNotMoveInd,'stable');

    % remove peaks during dropInd
    pkInds = setdiff(pkInds, fictracSmo.dropInd, 'stable');

    % convert pkInds to logical (for fictrac indices)
    pkLog = false(size(fictracSmo.t));
    pkLog(pkInds) = true;

    % peaks have to be greater than minYawThresh
    minCondLog = angVelSmoS > minYawThresh;
    pkLog = pkLog & minCondLog;

    % convert logical back to indices
    pkInds = find(pkLog);

    % all indices where yaw velocity is less than min
    yawMinInd = find(angVelSmoS < minYawThresh);

    % preallocate start and end index vectors
    pkStartInd = zeros(size(pkInds));
    pkEndInd = zeros(size(pkInds));

    % find start and end for each peak
    % if start and end can't be found (too close to edge, remove peaks)
    rmvInd = [];
    for i = 1:length(pkInds)
        thisStart = find(yawMinInd < pkInds(i), 1, 'last');
        % if no index found, set peak start to beginning of trial
        if isempty(thisStart)
            rmvInd = [rmvInd i];
        else
            % +1 because presumably next index is greater than min
            pkStartInd(i) = yawMinInd(thisStart) + 1;
        end

        thisEnd = find(yawMinInd > pkInds(i), 1, 'first');
        % if no index found, set peak end to end of trial
        if isempty(thisEnd)
            rmvInd = [rmvInd i];
        else
            % -1 because this is first that is less than min
            pkEndInd(i) = yawMinInd(thisEnd) - 1;
        end
    end
    % remove peaks too close to edge
    pkInds(rmvInd) = [];
    pkStartInd(rmvInd) = [];
    pkEndInd(rmvInd) = [];

    initSizeStart = length(pkStartInd);

    % remove any peaks whose edges fall within dropInd
    [pkStartInd, startNotRmvInd] = setdiff(pkStartInd, fictracSmo.dropInd, ...
        'stable');
    startRmvInd = setxor(1:initSizeStart, startNotRmvInd);
    pkInds(startRmvInd) = [];
    pkEndInd(startRmvInd) = [];

    initSizeEnd = length(pkEndInd);
    [pkEndInd, endNotRmvInd] = setdiff(pkEndInd, fictracSmo.dropInd, ...
        'stable');
    endRmvInd = setxor(1:initSizeEnd, endNotRmvInd);
    pkInds(endRmvInd) = [];
    pkStartInd(endRmvInd) = [];




    % check if any of the peaks share edges - look only at peak starts
    %  (will be the same as peak ends)

    % get indices into bodytraj, for not unique starts 
    [~, uniInd]  = unique(pkStartInd, 'stable');
    nonUniInd = setdiff(1:numel(pkStartInd), uniInd);
    startNonUniInd = pkStartInd(nonUniInd);

    % indices of peaks themselves, for non-unique peaks
    peakNonUniInd = pkInds(nonUniInd); 


    % if peaks share edges, keep only peak with greatest yaw velocity
    if (~isempty(startNonUniInd))
        notKeepInd = []; % initialize tracker of peaks to discard
        for i = 1:length(startNonUniInd)
            thisPeakStart = startNonUniInd(i);

            % find all peaks that share this start
            sharedInd = find(pkStartInd == thisPeakStart);

            % find peak with greatest yaw velocity
            bInd = pkInds(sharedInd);
            yawVals = angVelSmoS(bInd);

            [~, keepPkInd] = max(yawVals);
            keepPkBInd = bInd(keepPkInd);
            notKeepPkBInd = setdiff(bInd, keepPkBInd, 'stable');

            if (iscolumn(notKeepPkBInd))
                notKeepPkBInd = notKeepPkBInd';
            end

            % track all peaks to discard
            notKeepInd = [notKeepInd notKeepPkBInd];
        end

        % discard not keep peaks
        % discard for peaks, get indices of ones to keep
        [pkInds, keepInd] = setdiff(pkInds,notKeepInd,'stable');
        pkStartInd = pkStartInd(keepInd);
        pkEndInd = pkEndInd(keepInd);
        pkInds = pkInds(keepInd);
    end

    % loop through all peaks and check that the bout meets the duration
    %  requirements
    % initialize to keep track of peaks to remove
    rmInd = []; 

    for i = 1:length(pkInds)

        thisBoutStartT = fictracSmo.t(pkStartInd(i));
        thisBoutEndT = fictracSmo.t(pkEndInd(i));

        thisBoutDur = thisBoutEndT - thisBoutStartT;

        % if this bout duration is too short or too long, flag this index
        %  for deletion
        if((thisBoutDur < turnDur(1)) || ...
                (thisBoutDur > turnDur(2)))
            rmInd = [rmInd i];
        end
    end
    % remove any bouts that don't meet the duration criteria
    pkInds(rmInd) = [];
    pkStartInd(rmInd) = [];
    pkEndInd(rmInd) = [];

    % outputs, indices
    yawVelPeakInd = pkInds;
    boutStartInd = pkStartInd;
    boutEndInd = pkEndInd;

    % outputs, times
    yawVelPeakTimes = fictracSmo.t(pkInds);
    boutStartTimes = fictracSmo.t(pkInds);
    boutEndTimes = fictracSmo.t(pkInds);
end