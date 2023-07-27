% plotMultiCorr.m
%
% Function to plot output from saveStepParamMultiCorr(), 1 file
% Plots randomized as mean +/- SD and actual as single line
% 
% INPUTS:
%   datDir - path to folder containing correlation files
%   yScale - scale for correlations
%
% OUTPUTS:
%   none, but generates plot
%
% CREATED: 7/26/23 - HHY
%
% UPDATED:
%   7/26/23 - HHY
%
function plotMultiCorr(datDir, yScale)

    % prompt user to select cond_bout() file
    [condBoutFName, condBoutPath] = uigetfile('*.mat', ...
        'Select output file', datDir, 'MultiSelect', 'off');

    % load data from cond_bout file
    fullFilePath = [condBoutPath filesep condBoutFName];

    load(fullFilePath, 'allVarNames', 'corrVals', 'corrValsRand');

    % get mean and SD of random corr
    randCorrMean = mean(corrValsRand,2);
    randCorrStd = std(corrValsRand, [], 2);

    % plot
    figure;

    % indices for variable names
    xInd = 1:length(allVarNames);

    % plot rand dist
    errorbar(xInd,randCorrMean,randCorrStd, '.');

    hold on;

    % plot actual
    plot(xInd, corrVals, 'Marker', '_', 'LineStyle','none');

    % label x-axis
    xticks(xInd);
    xticklabels(allVarNames);

    % yScale
    ylim(yScale);

    % xScale
    xlim([xInd(1) - 0.5, xInd(end) + 0.5]);

    % title
    title('R-squared')
end
