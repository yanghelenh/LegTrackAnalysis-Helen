% extractTrialsOpto.m
%
% Helper function for computeLegFictracOpto_1Fly() that takes in one
%  trial's worth of data for one behavior parameter and divides it into
%  reps based on the opto stimulation
% Output appends trials into cell array given as input 
% Has option to align position data by subtracting value at stimulation
%  start
%
% REPLACES extractFictracTrialsOpto()
%
% INPUTS:
%   reps - cell array of size # NDs x # durations, to append individual 
%       trials into
%   repsOptoTimes - cell array of size # NDs x # durations, to keep
%       track of opto stim start times for each trial, to append into
%   behVar - behavior variable values
%   t - time for behVar
%   opto - struct of opto stimulation data
%   NDs - vector of all NDs 
%   durs - vector of all durations of stimulation
%   bwStimDur - scalar value of time between stimulations to consider, in
%       seconds, rep will be stimulation plus this time before and after
%   norm2StimStart - logical for whether to normalize to value at
%       stimulation start
%
% OUTPUTS:
%   reps - same cell array as inputs, with additional reps added
%   repsOptoTimes - same cell array as inputs, with additional opto
%       stim start times added
%   durTs - cell array of same size as durs, with time at each sample point
%       for each duration; stimulation starts at 0 (in sec)
%
% CREATED: 3/29/22 - HHY
%
% UPDATED:
%   3/29/22 - HHY
%   5/13/22 - HHY - made this function to replace
%       extractFictracTrialsOpto(), by adding option to subtract value at
%       stimulation start
%   5/18/22 - HHY - remove portion that removes trials with NaNs to
%       preserve rep indexing b/w FicTrac and leg tracking
%   5/19/22 - HHY add variable repsOptoStartTime to keep track of
%       optogenetic stimulation start times for each rep
%
function [reps, repsOptoTimes, durTs] = extractTrialsOpto(reps, ...
    repsOptoTimes, behVar, t, opto, NDs, durs, bwStimDur, ...
    norm2StimStart)

    % check that behVar is column vector; if not, make it so
    if (isrow(behVar))
        behVar = behVar';
    end

    % get ND for this trial
    thisND = opto.stimParams.ndFilter;

    % get index corresponding to this ND
    thisNDInd = find(NDs == thisND);

    % if this ND isn't one of the specified ones, return
    if (isempty(thisNDInd))
        return;
    end

    % get ifi for behavior var
    ifi = mean(diff(t));

    % convert stimulation durations and between stimulation duration into 
    %  number of behavior samples
    dursNumSamps = floor(durs / ifi);
    bwStimNumSamps = floor(bwStimDur / ifi);

    % get times for each duration
    durTs = cell(size(durs));

    for i = 1:length(durs)
        thisDurSamps = 1:(dursNumSamps(i) + bwStimNumSamps * 2);

        durTs{i} = ((thisDurSamps - 1) * ifi) - (bwStimNumSamps * ifi);
    end

    % loop through all stimulations
    for i = 1:length(opto.stimCmdDurs)
        % get duration for this trial
        thisDur = opto.stimCmdDurs(i);

        % get index corresponding to this duration
        thisDurInd = find(durs == thisDur);

        % if this duration isn't one of the specified ones, go to next
        %  stimulation
        if (isempty(thisDurInd))
            continue;
        end

        % get behavior index corresponding to stimulation start
        stimStartInd = find(t >= opto.stimStartTimes(i),1,'first');

        % index for rep start (stim start with b/w stim period appended)
        repStartInd = stimStartInd - bwStimNumSamps;
        % index for rep end (stim start + stim duration + b/w stim period)
        repEndInd = stimStartInd + dursNumSamps(thisDurInd) + bwStimNumSamps - 1;

        % check that start index is valid (1 or greater), otherwise, skip
        %  this rep
        if (repStartInd < 1)
            continue;
        % check that end index is valid (not greater than length of whole
        %  trial; otherwise, skip this rep
        elseif (repEndInd > length(t))
            continue;
        end

        % get this rep behavior data
        thisRepBeh = behVar(repStartInd:repEndInd);
        thisRepStartVal = behVar(stimStartInd);

%         Removed this 5/18/22 - trials will match b/w FicTrac and leg
%          tracking, by index 
%         % check for NaNs in behavior data; if present, skip this rep
%         if (any(isnan(thisRepBeh)))
%             continue;
%         end

        % check if behavior data should be normalized to value at start of
        % stimulation (subtract out that value)
        if (norm2StimStart)
            thisRepBeh = thisRepBeh - thisRepStartVal;
        end

        % append this rep behavior data to data cell array
        % each rep is a row
        reps{thisNDInd,thisDurInd} = [reps{thisNDInd,thisDurInd}; ...
            thisRepBeh'];

        % append this rep opto stim start time to appropriate cell array
        repsOptoTimes{thisNDInd,thisDurInd} = [...
            repsOptoTimes{thisNDInd,thisDurInd}; ...
            opto.stimStartTimes(i)];
    end
end