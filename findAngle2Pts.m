% findAngle2Pts.m
%
% Function to find the angle between two points in x,y plane. Angle of the
%  vector that points from (x1, y1) to (x2, y2). In degrees
% Helper function for leg analysis, uses leg coordinate frame, where x axis
%  is direction parallel to the long axis of the body and y axis is
%  direction parallel to the short axis of the body. Negative x is the
%  head. Positive y is the fly's right.
% Angles defined such that angles to the right are positive and angles to
%  the left are negative. 0 degrees is pointing straight up to the head. 90
%  degrees is straight to the right side.
% Operates on single pair of points or vectors defining paired (x,y) values
% Returns NaN when 2 points are identical
%
% INPUTS:
%   x1 - vector of x1 coordinates of points (or single x1 value)
%   y1 - vector of y1 coordinates of points (or single y1 value)
%   x2 - vector of x2 coordinates of points (or single x2 value)
%   y2 - vector of y2 coordinates of points (or single y2 value)
%
% OUTPUTS:
%   ptAngs - vector of same length as inputs giving angles between all
%       points
%
% CREATED: 11/9/20 - HHY
%
% UPDATED:
%   11/9/20 - HHY
%
function ptAngs = findAngle2Pts(x1, y1, x2, y2)

    % zero the points
    x = x2 - x1;
    y = y2 - y1;
    
    % get absolute values, for calculating tangent between 0 and 90
    %  degrees; sign accounted for below
    xAbs = abs(x);
    yAbs = abs(y);
    
    % use atan to find angles
    ptAngs = atand(yAbs./xAbs);
    
    % map angles to correct sign given signs of inputs
    % negative x, positive y - don't do anything
    % positive x, positive y - 180 minus angle
    posXPosY = find((x > 0) & (y > 0));
    ptAngs(posXPosY) = 180 - ptAngs(posXPosY);
    % negative x, negative y - flip sign
    negXNegY = find((x <= 0) & (y < 0));
    ptAngs(negXNegY) = ptAngs(negXNegY) .* -1;
    % positve x, negative y - angle minus 180
    posXNegY = find((x > 0) & (y < 0));
    ptAngs(posXNegY) = ptAngs(posXNegY) - 180;
end