% findBoutsDurStim2Rmv.m
%
% Helper function for saveBallLegStepParamCond_bouts() that takes in start
%  and end times of turning bouts as well as start and end times of
%  stimulation bouts and returns the indices of turning bouts that overlap
%  with the stimulation bouts. These turning bouts are to be removed from
%  consideration.
%
% INPUTS:
%   boutStartTimes - start times of turning bouts, vector length n
%   boutEndTimes - end times of turning bouts, matched to start times,
%       length n
%   stimStartTimes - start times of stimulation, vector length m
%   stimEndTimes - end times of stimulation bouts, matched to start times,
%       length m
%   postStimExclDur - additional time, in sec, after stimulation to also
%       exclude turning bouts from overlapping with
%
% OUTPUTS:
%   boutInd2Rmv - indices into boutStartTimes/boutEndTimes of bouts to
%       remove b/c they overlap with the given stimulation
%
% CREATED: 6/22/23 - HHY
%
% UPDATED:
%   6/22/23 - HHY
%
function boutInd2Rmv = findBoutsDurStim2Rmv(boutStartTimes, ...
    boutEndTimes, stimStartTimes, stimEndTimes, postStimExclDur)

    % get end times to consider, with additional buffer
    stimEndTimesBuff = stimEndTimes + postStimExclDur;

    % initialize vector tracking bouts to remove
    boutInd2Rmv = [];

    % loop through all turning bouts, remove ones that overlap with
    %  stimulation bouts
    for i = 1:length(boutStartTimes)
        thisBoutStart = boutStartTimes(i);
        thisBoutEnd = boutEndTimes(i);

        % check if this bout start is greater than any stim start and also
        %  less than the corresponding stim end
        if(any((stimStartTimes <= thisBoutStart) & ...
            (stimEndTimesBuff >= thisBoutStart)))
            % if yes, flag this bout for removal
            boutInd2Rmv = [boutInd2Rmv; i];
        
        % check if this bout end is less than any stim end and also greater 
        %  than the correspond stim start
        elseif(any((stimEndTimesBuff >= thisBoutEnd) & ...
            (stimStartTimes <= thisBoutEnd)))
            % if yes, flag this bout for removal
            boutInd2Rmv = [boutInd2Rmv; i];
        end
    end
    
end