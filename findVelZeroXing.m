% findVelZeroXing.m
%
% Function to find zero crossings in velocity, while specifying a minimum 
%  number of samples spent between zero crossings.
% Finds zero crossings after removing not-moving times, which would
%  otherwise confuse zero-crossing calls
% Specifies whether zero crossing is with velocity going negative to
%  positive or positive to negative
%
% INPUTS:
%   vel - matrix of velocities, num samples x num labeled pts
%   zeroVelInd - indicies of zero velocity; removed from consideration in
%       zero-crossings call
%   zeroXingParams - struct of parameters, with fields
%       legInd - indices into vel labeled pts, for columns to find zero
%           crossings on
%       legNames - names for each vel column corresponding to legInd
%       minStepDur - minimum duration between zero crossings, in samples
%
% OUTPUTS:
%   zeroXing - struct of zero crossings with fields
%       neg2Pos - indicies of negative to positive zero crossings, as a
%           separate field for each legInd, named legNames
%       pos2Neg - indicies of positive to negative zero crossings, with
%           fields as neg2Pos
%
% CREATED: 11/18/20 - HHY
%
% UPDATED:
%   11/18/20 - HHY
%   8/5/21 - HHY - fix bug in identifying bouts of (+) and (-) velocities
%       (shifted by 1)
%   8/6/21 - HHY - DON'T USE! there is some bug in calling zero crossings,
%       maybe having to do with bout edges, but unclear yet
%
function zeroXing = findVelZeroXing(vel, zeroVelInd, zeroXingParams)

    % get indicies when fly is moving (inverse of zeroVelInd), when to find
    %  zero crossings
    allInd = (1:length(vel(:,2)))';
    moveInd = setdiff(allInd, zeroVelInd);

    % find moving bouts
    [moveBoutStarts, moveBoutEnds, ~] = findBouts(moveInd);

    % preallocate struct for saving zero crossing info (needs to be struct b/c
    %  number of zero crossings not the same across legs)
    for i = 1:length(zeroXingParams.legInd)
        zeroXing.neg2Pos.(zeroXingParams.legNames{zeroXingParams.legInd(i)}) = [];
        zeroXing.pos2Neg.(zeroXingParams.legNames{zeroXingParams.legInd(i)}) = [];
    end

    % loop through all moving bouts
    for i = 1:length(moveBoutStarts)
        boutInd = (moveBoutStarts(i):moveBoutEnds(i))';

        % loop through all legs
        for j = 1:length(zeroXingParams.legInd)
            % find indicies of all positive velocities
            posInd = find(...
                vel(boutInd,zeroXingParams.legInd(j)) >= 0) ...
                + boutInd(1) - 1;

            % find indicies of negative velocities
            negInd = find(...
                vel(boutInd,zeroXingParams.legInd(j)) < 0) ...
                + boutInd(1) - 1;

            % moving bout must have positive and negative velocities to 
            %  count
            if (~isempty(posInd) && ~isempty(negInd))
                % find bout starts and ends
                [posBoutStarts, posBoutEnds, posBoutDur] = ...
                    findBouts(posInd);
                [negBoutStarts, negBoutEnds, negBoutDur] = ...
                    findBouts(negInd);

                % find indicies that bracket zero crossing
                % negative to postive -> negBoutEnd + 1 = posBoutStart;
                %  positive bout duration defined by boutStart must be >
                %  minimum
                % postive to negatvie -> posBoutEnd + 1 = negBoutStart; 
                %  negative bout duration defined by boutStart must be >
                %  minimum
                neg2PosZeroXingLog = ismember(posBoutStarts, (negBoutEnds+1));
                neg2PosZeroXingEnds = ...
                    posBoutStarts(posBoutDur(neg2PosZeroXingLog) > ...
                    zeroXingParams.minStepDur);
                neg2PosZeroXingStarts = neg2PosZeroXingEnds - 1;

                pos2NegZeroXingLog = ismember(negBoutStarts, (posBoutEnds+1));
                pos2NegZeroXingEnds = ...
                    negBoutStarts(negBoutDur(pos2NegZeroXingLog) > ...
                    zeroXingParams.minStepDur);
                pos2NegZeroXingStarts = pos2NegZeroXingEnds - 1;

                % append to existing
                zeroXing.neg2Pos.(zeroXingParams.legNames{zeroXingParams.legInd(j)})...
                    = ...
                    [zeroXing.neg2Pos.(zeroXingParams.legNames{zeroXingParams.legInd(j)}); ...
                    neg2PosZeroXingStarts];
                zeroXing.pos2Neg.(zeroXingParams.legNames{zeroXingParams.legInd(j)})...
                    = ...
                    [zeroXing.pos2Neg.(zeroXingParams.legNames{zeroXingParams.legInd(j)}); ...
                    pos2NegZeroXingStarts];
            end
        end
    end

end