% rotatePts2D.m
%
% Function that takes in (x,y) coordiate of a point or set of points and
%  rotates them the specified angle (in degrees). Works in 2D. I.e. xy
%  plane
% Rotates point using rotation matrix
%
% INPUTS:
%   xPt - x coordinate(s) of point(s) as scalar or vector of values
%   yPt - y coordinate(s) of point(s) as scalar or vector of values paired
%       with xPt
%   rotAng - rotation angle, in degrees (positive rotation is
%       counterclockwise, as is Cartesian convention)
%
% OUTPUTS:
%   xPtRot - x coordiate(s) of rotated point(s)
%   yPtRot - y coordiate(s) of rotated point(s)
%
% CREATED: 10/19/20 - HHY
%
% UPDATED: 
%   10/19/20 - HHY
%
function [xPtRot, yPtRot] = rotatePts2D(xPt, yPt, rotAng)

    % check that input vectors are same length
    if ~(length(xPt) == length(yPt))
        xPtRot = NaN;
        yPtRot = NaN;
        disp('xPt and yPt must be same length');
        return; % end function here
    end
    
    % define rotation matrix
    rotMatrix = [cosd(rotAng) -sind(rotAng); sind(rotAng) cosd(rotAng)];
    
    % preallocate
    xPtRot = zeros(size(xPt));
    yPtRot = zeros(size(yPt));
    
    % rotate point(s)
    for i=1:length(xPt)
        point = [xPt(i) yPt(i)]'; % (x,y) point as column vector
        
        % multiply by rotation matrix to rotate point
        rotPoint = rotMatrix * point;
        
        % save rotated point into output
        xPtRot(i) = rotPoint(1);
        yPtRot(i) = rotPoint(2);
        
    end

end