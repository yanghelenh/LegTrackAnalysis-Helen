% getLinProj.m
%
% Helper function to get linear projection of multiple variables onto
%  defined axis. Like PCA, but specify coefficients 
% Mean subtracts each variable but does not otherwise normalize values
%  across different variables 
%
% INPUTS:
%   varVals - matrix of size numSamps x numVars for sample values of each
%     variable at each pt
%   coeffs - coefficients that sum to 1 for weighting of each variable.
%     Size must be numVars
%
% OUTPUTS:
%   scores - vector of length numSamps, representing weighted projection of
%     variable values onto axis definted by coeffs
%
% CREATED: 8/8/23 - HHY
%
% UPDATED:
%   8/8/23 - HHY
%
function scores = getLinProj(varVals, coeffs)

    % number of samples, variables
    numSamps = size(varVals,1);
    numVars = size(varVals,2);
    
    % get mean for each variable
    varMeans = mean(varVals,1);

    % mean subtract variable values
    meanSubVarVals = varVals - repmat(varMeans,numSamps,1);

    % for each sample point, get its value as weighted projection for
    %  coefficient
    % make sure coeffs is column vector
    if ~iscolumn(coeffs)
        coeffs = coeffs';
    end

    scores = meanSubVarVals * coeffs;
end