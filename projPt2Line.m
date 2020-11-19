% projPt2Line.m
%
% Function that projects a point onto a line. I.e. returns point on line
%  that connects specified point to line with orthogonal line.
% User specifies slope and y-intercept of line, x and y coordinate of
%  point. If these are vectors, treats each set of indicies as one line and
%  point.
%
% INPUT:
%   slope - slope of line (can be vector of slopes of many lines)
%   yInt - y-intercept of line (can be vector of y-intercepts of many
%       lines)
%   ptX - x coordinate of point to project to line (can be vector)
%   ptY - y coordinate of point to project ot line (can be vector)
%
% OUTPUT:
%   projX - x coordinate of projected point on line
%   projY - y coordinate of projected point on line
%
% CREATED: 10/18/20 - HHY
%
% UPDATED:
%   10/18/20 - HHY
%
function [projX, projY] = projPt2Line(slope, yInt, ptX, ptY)
    % slope of projection line is orthogonal
    perpSlope = -1./slope;
    
    % y-intercept of orthogonal line (line goes through specified point)
    yIntNew = -perpSlope .* ptX + ptY;
    
    % solve for point on line that intercects orthogonal line
    projX = (yIntNew - yInt) ./ (slope - perpSlope); 
    % for some reason, solving for y using equation of perpendicular line
    %  doesn't always return correct answer, likely when slope too close to
    %  0
%     projY = perpSlope .* projX + yIntNew; 
    projY = slope .* projX + yInt;
end