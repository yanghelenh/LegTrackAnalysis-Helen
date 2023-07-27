% plotScatterPCAStepParams.m
%
% Function to plot scatterplot of output of pcaCompStepParams()
% Plots each turn bout set 1 PC1 vs set 2 PC1.
% Select output file of pcaCompStepParams() through GUI
%
% INPUTS:
%   datDir - path to folder containing pcaCompStepParams() output files
%
% OUTPUTS:
%   none, but generates plot
%
% CREATED: 7/27/23 - HHY
%
% UPDATED:
%   7/23/23 - HHY
%
function plotScatterPCAStepParams(datDir, numBins)

    % prompt user to select cond_bout() file
    [condBoutFName, condBoutPath] = uigetfile('*.mat', ...
        'Select output file', datDir, 'MultiSelect', 'off');

    % load data from cond_bout file
    fullFilePath = [condBoutPath filesep condBoutFName];

    load(fullFilePath, 'set1', 'set2');

    figure;

    % scatterplot all turn bouts in PC space
%     scatter(set1.score(:,1), set2.score(:,1), 'filled', ...
%         'MarkerFaceAlpha', 0.8, 'MarkerEdgeAlpha', 0.3);

    binscatter(set1.score(:,1), set2.score(:,1), numBins);

    hold on;

    % plot y = x line
    xLimits = xlim;
    yLimits = ylim;

    minPt = min([xLimits(1) yLimits(1)]);
    maxPt = max([xLimits(2) yLimits(2)]);

    plot([minPt maxPt], [minPt maxPt], 'k', 'LineWidth', 2);

    % plot x = 0 line
    line([0 0], [minPt maxPt], 'Color', 'k');
    % plot y = 0 line
    line([minPt maxPt], [0 0], 'Color', 'k');

    xAxisLabel = sprintf('%s, PC1: %.2f', set1.name, set1.explained(1));
    xlabel(xAxisLabel);

    yAxisLabel = sprintf('%s, PC1: %.2f', set2.name, set2.explained(1));
    ylabel(yAxisLabel);

    xlim(xLimits);
    ylim(yLimits);
end