% medFiltRmvOutliers.m
%
% Function to remove outliers from a 1D vector by median filtering
%  (dimension specified by user) and then replacing all outliers (defined
%  as values +/- specifed percentage away from median filtered value) with
%  median filtered value
% Helper function for preprocessLegTrack(). Run on leg position data.
%
% INPUTS:
%   vec - vector of data to filter, 1D array
%   medFiltDeg - degree of median filter, positive, odd integer
%   percntDev - percent deviation allowed from median filtered value;
%       values exceeding replaced; between 0 and 100
%   maxMinWin - number of samples around a point from which to extract the
%       max and the min, so that the % deviation can be turned into a threshold
%       value
%
% OUTPUTS:
%   filtVec - vector with outliers removed, same size as vec
%
% CREATED: 5/12/22 - HHY
%
% UPDATED:
%   5/12/22 - HHY
%
function filtVec = medFiltRmvOutliers(vec, medFiltDeg, percntDev, ...
    maxMinWin)

    % median filter vector
    medFiltVec = medfilt1(vec, medFiltDeg);

    % get half of maxMinWin, to nearest integer; for before and after
    winSize = ceil(maxMinWin/2);

    % preallocate
    accptDev = zeros(size(medFiltVec));

    % using values from median filtered vector, calculate threshold for
    %  deviation for outliers
    for i = 1:length(medFiltVec)
        % start and end indices for window around point for determining max
        %  and min values
        startInd = i - winSize;
        endInd = i + winSize;

        % if start index before vector start; start at vector start instead
        if (startInd < 1)
            startInd = 1;
        end

        % if end index after vector end; end at vector end instead
        if (endInd > length(medFiltVec))
            endInd = length(medFiltVec);
        end

        % window in which to determine max/min
        win = medFiltVec(startInd:endInd);

        % get acceptable deviation; to make threshold for outliers
        accptDev(i) = (max(win) - min(win)) * (percntDev/100);
    end

    % find deviation of actual data from median filtered data; positive
    %  values only
    filtActDiff = abs(vec - medFiltVec);

    % logical for values that deviate more than acceptable
    devLog = filtActDiff > accptDev;

    % replace values that deviate more than acceptable with median filtered
    %  value
    filtVec = vec;
    filtVec(devLog) = medFiltVec(devLog);
end