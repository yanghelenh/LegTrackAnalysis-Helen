% distBtw2Pts.m
%
% Function to return the Euclidian distance between two points, using the
%  vector norm.
% If input x and y points are themselves vectors, returns the distance
%  between each set of points, paired at corresponding indicies
% Inputs must be of same length
%
% INPUTS:
%   x1 - x value of first point (can be vector of x values)
%   y1 - y value of first point (can be vector of y values)
%   x2 - x value of second point (can be vector of x values)
%   y2 - y value of second point (can be vector of y values)
%
% OUTPUTS:
%   dist - Euclidian distance between (x1,y1) and (x2,y2)
%
% CREATED: 10/18/20 - HHY
%
% UPDATED:
%   10/18/20 - HHY
%

function dist = distBtw2Pts(x1, y1, x2, y2)

    % check that all inputs are vectors or scalars, 1xn or nx1 vector where
    %  n is nonnegative integer
    if (~isvector(x1) || ~isvector(y1) || ~isvector(x2) || ~isvector(y2))
        dist = NaN;
        disp('Inputs must be vectors or scalars, not matrices');
        return; % end function here
    end

    % check that all input vectors are the same length
    lenVec = [length(x1) length(y1) length(x2) length(y2)];
    numUniqueLengths = length(unique(lenVec));
    
    % if number of unique vector lengths is not 1, then not all vectors are
    %  the same length
    if (numUniqueLengths ~= 1)
        dist = NaN; 
        disp('Input vectors to distBtw2Pts must be same length');
        return; % end function here
    end
      
    % check that all vectors are column vectors, make them column
    %  vectors if they're not
    if ~iscolumn(x1)
        x1 = x1';
    end
    if ~iscolumn(y1)
        y1 = y1';
    end
    if ~iscolumn(x2)
        x2 = x2';
    end
    if ~iscolumn(y2)
        y2 = y2';
    end
    
    % pair up x and y values of points, such that each is row in nx2 matrix
    pt1 = [x1 y1];
    pt2 = [x2 y2];
    
    % preallocate distance output vector
    dist = zeros(size(x1));
    
    % loop through all pairs of points, returning distance
    for i = 1:length(x1)
        dist(i) = norm(pt2(i,:)-pt1(i,:));
    end
        
end