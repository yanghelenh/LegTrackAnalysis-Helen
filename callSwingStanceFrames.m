% callSwingStanceFrames.m
%
% Function to call swing phase and stance phase for each leg on each video 
%  frame.
% If the leg and ball move in the same direction, leg is in stance. If the
%  leg and the ball move in opposite directions, leg is in swing. So,
%  comparing leg velocity to FicTrac fwdVel, same sign is stance, opposite
%  signs is swing.
% Does not call stance/swing when fly is not moving.
%
% INPUTS:
%   fictrac - struct of FicTrac data, output from preprocessFictrac()
%   legVidFrameTimes - timing for leg video frames
%   legXVel - velocity of legs in X direction, as matrix, output of
%       findLegVel()
%   zeroVelInd - indicies when fly is not moving
%   legInd - indicies of leg tracked pts
%
% OUTPUTS:
%   legStance - logical for whether leg is in stance, for all legs; 
%       frames x legs matrix
%   legSwing - logical for whether leg is in swing, for all legs; 
%       frames x legs matrix
%   legSwingStanceNotMove - matrix for swing (-1), stance (1), not moving
%       (0)
%
% CREATED: 11/19/20 - HHY
%
% UPDATED:
%   11/19/20 - HHY
%
function [legStance, legSwing, legSwingStanceNotMove] = ...
    callSwingStanceFrames(fictrac, legVidFrameTimes, legXVel, ...
    zeroVelInd, legInd)

    % determine stance vs. swing: only when fly is moving, compare to 
    %  direction of fictrac forward velocity
    % interpolate fictrac forward velocity to legVidFrametimes
    interpFtFwd = interp1(fictrac.t, fictrac.fwdVel, legVidFrameTimes);

    % binarize fictrac forward velocity
    ftFwdInd = interpFtFwd >= 0; % logical
    % convert to +1 for moving forward, -1 for moving backwards
    ftFwd = zeros(size(legXVel,2),1);
    ftFwd(ftFwdInd) = 1; % fwd = 1
    ftFwd(~ftFwdInd) = -1; % not fwd = -1

    % binarize leg velocity
    binLegXVel = legXVel >=0;
    % covert to +1 for moving forward, -1 for moving backwards
    legXVelDir = zeros(size(legXVel));
    legXVelDir(binLegXVel) = 1;
    legXVelDir(~binLegXVel) = -1;

    % preallocate stance and swing matricies
    legStance = false(size(legXVel,1), length(legInd));
    legSwing = false(size(legXVel,1), length(legInd));
    % non-logical that indicates swing, stance, or not moving in 1 vector:
    %  1 is stance, -1 is swing, 0 is not moving
    legSwingStanceNotMove = zeros(size(legXVel,1), ...
        length(legInd));

    % stance is leg and fictrac moving in same direction (multiply to 1),
    %  swing is leg and fictrac moving in opposite directions (multiply 
    %  to -1)
    % assign not moving bouts to be 0
    % loop through all legs
    for i = 1:length(legInd)
        legFtDir = legXVelDir(:,legInd(i)) .* ftFwd;
        legFtDir(zeroVelInd) = 0;

        legStance(legFtDir==1,legInd(i))=true;
        legSwing(legFtDir==-1,legInd(i))=true;

        legSwingStanceNotMove(:,legInd(i)) = legFtDir;

    end

end