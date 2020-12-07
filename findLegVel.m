% findLegVel.m
%
% Function that takes in shifted, rotated, and normalized leg position and
%  returns leg velocity.
% Performs Gaussian process smoothing on leg position prior to calculating
%  velocity
%
% INPUTS:
%   srnLeg - shifted, rotated, normalized leg position matrix (output of
%    shiftRotateNormalizeLegPos())
%   smoParams - parameters for Gaussian process smoothing
%       padLen
%       sigma
%
% OUTPUTS:
%   legVel - leg velocity, as matrix of same size as srnLeg
%
% CREATED: 11/18/20 - HHY
% 
% UPDATED:
%   11/18/20 - HHY
%
function legVel = findLegVel(srnLeg, smoParams)
    % preallocate
    srnLegSmo = zeros(size(srnLeg));

    % loop through all tracked points
    for i = 1:size(srnLegSmo,2)
        smoPosPy = py.proc_utils.safe_interp_conv_smooth_ball_data(...
            srnLeg(:,i)', smoParams.padLen, smoParams.sigma);
        % convert from python to matlab data format
        smoPos = cell2mat(cell(smoPosPy.tolist()));
        srnLegSmo(:,i) = smoPos';
    end

    % get leg velocities by taking gradient
    legVel = zeros(size(srnLegSmo));
    for i = 1:size(legVel,2)
        legVel(:,i) = gradient(srnLegSmo(:,i));
    end
end