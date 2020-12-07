% shiftRotateNormalizeLegPos.m
%
% Function that takes in raw leg X and Y coordinates (output from
%  loadTrkFile) and returns coordinates that are shifted and rotated to
%  align all the frames to the center point of the fly (midpoint between 2
%  midlegs) and normalized so the units are body lengths (instead of
%  pixels).
%
% INPUTS:
%   legX - X coordinates of leg positions (frames x points), output from 
%       loadTrkFile
%   legY - Y coordinates of leg positions (frames x points), output from 
%       loadTrkFile
%   refPts - struct of indicies corresponding to reference points
%       lineFitInd - indicies for points for fitting line of body
%       midPtInd - index of midpoint of fly
%       headPtInd - index of head of fly
%       abdPtInd - index of abdomen tip of fly
%
% OUTPUTS:
%   srnLegX - X coordinates of leg positions (frames x points), after
%       shifting and rotating to align and normalizing to body lengths 
%       units
%   srnLegY - Y coordinates of leg positions (frames x points), after
%       shifting and rotating to align and normalizing to body lengths 
%       units
%
% CREATED: 11/16/20 - HHY
%
% UPDATED:
%   11/16/20 - HHY
%
function [srnLegX, srnLegY] = shiftRotateNormalizeLegPos(legX, legY, ...
    refPts)

    % fit line to 4 points defining head and 3 thorax-coxa joints; [7-10]

    % preallocate vector for saving coefficients of line fit
    bodyLineCoeffs = zeros(size(legX,1), 2);

    % get coefficients of line for each video frame; first is slope, second
    %  is y-intercept
    for i = 1:size(legX, 1)
        xPts = legX(i, refPts.lineFitInd);
        yPts = legY(i, refPts.lineFitInd);

        bodyLineCoeffs(i,:) = polyfit(xPts, yPts, 1);
    end

    % get projection of midleg thorax-coxa point onto body line (define 
    %  this as midpoint of fly)
    [midXProj, midYProj] = projPt2Line(bodyLineCoeffs(:,1), ...
        bodyLineCoeffs(:,2), legX(:,refPts.midPtInd), ...
        legY(:,refPts.midPtInd));

    % project head point and abdomen point onto body line, get distance
    %  between them to get body length
    [headXProj, headYProj] = projPt2Line(bodyLineCoeffs(:,1), ...
        bodyLineCoeffs(:,2), legX(:,refPts.headPtInd), ...
        legY(:,refPts.headPtInd));
    [abdXProj, abdYProj] = projPt2Line(bodyLineCoeffs(:,1), ...
        bodyLineCoeffs(:,2), legX(:,refPts.abdPtInd), ...
        legY(:,refPts.abdPtInd));

    bodyLen = distBtw2Pts(headXProj, headYProj, abdXProj, abdYProj);

    % shift all points so that they're relative to the midpoint, which is
    %  now defined as (0,0)
    % preallocate
    legXShift = zeros(size(legX));
    legYShift = zeros(size(legY));
    % shift each point
    for i = 1:size(legX,2)
        legXShift(:,i) = legX(:,i) - midXProj;
        legYShift(:,i) = legY(:,i) - midYProj;
    end

    % get body line coefficients, after shift
    bodyLineShiftCoeffs = zeros(size(legXShift,1),2);

    for i = 1:size(legXShift, 1)
        xPts = legXShift(i, refPts.lineFitInd);
        yPts = legYShift(i, refPts.lineFitInd);

        bodyLineShiftCoeffs(i,:) = polyfit(xPts, yPts, 1);
    end

    % for whole time series, define the body length and the body angle
    % for each frame, rotate around the midpoint by the body angle so x and
    %  y axes align with forward and lateral axes of fly
    % normalize units to body length (from pixels)

    % body length as median of body lengths throughout trial
    flyBodyLength = median(bodyLen);

    % tangent of body angle = median slope of body line, 
    %  post shift (y-int now ~0)
    flyBodyAng = atand(median(bodyLineShiftCoeffs(:,1)));

    % rotate all points about midpt, now (0,0), by inverse of fly body 
    %  angle
    % preallocate
    legXShiftRot = zeros(size(legXShift));
    legYShiftRot = zeros(size(legYShift));
    % loop through each marked point and rotate
    for i = 1:size(legXShift,2)
        [legXShiftRot(:,i), legYShiftRot(:,i)] = ...
            rotatePts2D(legXShift(:,i), legYShift(:,i), -1*flyBodyAng);
    end

    % normalize to body length
    srnLegX = legXShiftRot / flyBodyLength;
    srnLegY = legYShiftRot / flyBodyLength;

    % multiply srnLegY by -1 so that fly's right side has positive
    %  coordinates and left side has negative coordinates
    srnLegY = -1 * srnLegY;
end