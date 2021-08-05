% callSwingStanceFrames_v2.m
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
%   legYVel - velocity of legs in Y direction, as matrix, output of
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
% CREATED: 12/7/20 - HHY
%
% UPDATED:
%   12/7/20 - HHY
%
function [legStance, legSwing, legSwingStanceNotMove] = ...
    callSwingStanceFrames_v2(fictrac, legVidFrameTimes, legXVel, ...
    legYVel, zeroVelInd, legInd)

    % determine stance vs. swing: only when fly is moving, compare to 
    %  direction of fictrac forward velocity
    % interpolate fictrac forward velocity to legVidFrametimes
    interpFtFwd = interp1(fictrac.t, fictrac.fwdVel, legVidFrameTimes);
    % interpolate fictrac yaw angular velocity to legVidFrametimes
    interpFtYaw = interp1(fictrac.t, fictrac.yawAngVel, legVidFrameTimes);

    % binarize fictrac forward velocity
    ftFwdInd = interpFtFwd >= 0; % logical
    % convert to +1 for moving forward, -1 for moving backwards
    ftFwd = zeros(size(legXVel,2),1);
    ftFwd(ftFwdInd) = 1; % fwd = 1
    ftFwd(~ftFwdInd) = -1; % not fwd = -1
    
    % binarize fictrac yaw angular velocity
    ftYawInd = interpFtYaw >= 0; % logical
    % convert to +1 for fly moving right, -1 for moving left
    ftYaw = zeros(size(legXVel,2),1);
    ftYaw(ftYawInd) = 1; % right = 1
    ftYaw(~ftYawInd) = -1; % left = -1

    % binarize leg x velocity
    binLegXVel = legXVel >=0;
    % convert to +1 for moving front to back, -1 for moving back to front
    legXVelDir = zeros(size(legXVel));
    legXVelDir(binLegXVel) = 1;
    legXVelDir(~binLegXVel) = -1;
    
    % binarize leg y velocity
    binLegYVel = legYVel <=0;
    % convert to +1 for moving right to left; -1 for moving left to right
    legYVelDir = zeros(size(legYVel));
    legYVelDir(binLegYVel) = 1; 
    legYVelDir(~binLegYVel) = -1;

    % preallocate stance and swing matricies
    legStance = false(size(legXVel,1), length(legInd));
    legSwing = false(size(legXVel,1), length(legInd));
    % non-logical that indicates swing, stance, or not moving in 1 vector:
    %  1 is stance, -1 is swing, 0 is not moving
    legSwingStanceNotMove = zeros(size(legXVel,1), ...
        length(legInd));

    % stance is leg and ball moving in same direction (with above sign 
    %  conventions, multiply to 1),
    %  swing is leg and ball moving in opposite directions (multiply 
    %  to -1)
    % assign not moving bouts to be 0
    % loop through all legs
    for i = 1:length(legInd)
        % compare leg x with ball forward
        legXFtDir = legXVelDir(:,legInd(i)) .* ftFwd;
        % compare leg y with ball yaw
        legYFtDir = legYVelDir(:,legInd(i)) .* ftYaw;
        
        % not moving = 0
        legXFtDir(zeroVelInd) = 0;
        legYFtDir(zeroVelInd) = 0;
        
        % vector incorporating both x and y
        legFtDir = legXFtDir; % initialize as x
        
        % check if swing/stance calls agree for leg x and y directions
        % indicies for where 2 are different
        diffLegXYFtDir = find((legXFtDir - legYFtDir) ~= 0);
        % loop through all indicies where 2 are different, 
        for j = 1:length(diffLegXYFtDir)
            % how much leg is moving in x and y (sign independent)
            xWeight = abs(legXVel(diffLegXYFtDir(j),legInd(i)));
            yWeight = abs(legYVel(diffLegXYFtDir(j),legInd(i)));
            
            % if leg is moving more in x, use x stance/swing call
            if (xWeight >= yWeight)
                legFtDir(diffLegXYFtDir(j)) = legXFtDir(diffLegXYFtDir(j));
            else % if leg is moving more in y, use y stance/swing call
                legFtDir(diffLegXYFtDir(j)) = legYFtDir(diffLegXYFtDir(j));
            end
        end

        legStance(legFtDir==1,legInd(i))=true;
        legSwing(legFtDir==-1,legInd(i))=true;

        legSwingStanceNotMove(:,legInd(i)) = legFtDir;

    end

end