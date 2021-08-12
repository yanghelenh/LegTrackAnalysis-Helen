% smoLegPos.m
%
% Function that takes in shifted, rotated, and normalized leg position and
%  smooths it using Gaussian process smoothing.
%
% INPUTS:
%   srnLeg - shifted, rotated, normalized leg position matrix (output of
%    shiftRotateNormalizeLegPos())
%   smoParams - parameters for Gaussian process smoothing
%       padLen
%       sigma
%
% OUTPUTS:
%   srnLegSmo - smoothed leg position, as matrix of same size as srnLeg
%
% CREATED: 8/6/21 - HHY
% 
% UPDATED:
%   8/6/21 - HHY
%
function srnLegSmo = smoLegPos(srnLeg, smoParams)
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
end