% plotScatterPCAStepParams.m
%
% Function to plot scatterplot of output of pcaCompStepParams()
% Plots each turn bout set 1 PC1 vs set 2 PC1.
% Select output file of pcaCompStepParams() through GUI
%
% INPUTS:
%   datDir - path to folder containing pcaCompStepParams() output files
%   plotType - string for type of plot to make
%       'binScatter' - binned scatterplot, counts
%       'scatter' - regular scatterplot
%       'peakFwd' - scatter colored by bout peak forward velocity
%       'peakYaw' - scatter colored by  bout peak yaw velocity
%       'peakLat' - scatter colored by peak lateral velocity, free walk
%       'peakSlide' - scatter colored by peak slide velocity, ball walk
%   numBins - for 'binScatter' type plot, number of bins per dimension
%   xyLine - boolean for whether to plot y = x line
%
% OUTPUTS:
%   none, but generates plot
%
% CREATED: 7/27/23 - HHY
%
% UPDATED:
%   7/23/23 - HHY
%
function plotScatterPCAStepParams(datDir, plotType, numBins, xyLine)

    % prompt user to select cond_bout() file
    [condBoutFName, condBoutPath] = uigetfile('*.mat', ...
        'Select output file', datDir, 'MultiSelect', 'off');

    % load data from cond_bout file
    fullFilePath = [condBoutPath filesep condBoutFName];

    load(fullFilePath, 'set1', 'set2', 'boutPeakVel');

    figure;

    % scatterplot all turn bouts in PC space
    if (strcmpi(plotType, 'binScatter'))
        binscatter(set1.score(:,1), set2.score(:,1), numBins);
    elseif (strcmpi(plotType, 'peakFwd'))
        % convert to mm/s if needed
        if (median(boutPeakVel.fwd) < 1)
            fwdVel = boutPeakVel.fwd * 1000;
        else
            fwdVel = boutPeakVel.fwd;
        end
        scatter(set1.score(:,1), set2.score(:,1), [], ...
            fwdVel, 'filled', ...
            'MarkerFaceAlpha', 0.8, 'MarkerEdgeAlpha', 0.3);
        colormap(gca, 'cool');
%       caxis([5 10]);
        colorbar;
    elseif (strcmpi(plotType, 'peakYaw'))
        scatter(set1.score(:,1), set2.score(:,1), [], ...
            abs(boutPeakVel.yaw), 'filled', ...
            'MarkerFaceAlpha', 0.8, 'MarkerEdgeAlpha', 0.3);
        colormap(gca, 'cool');
        caxis([0 300]);
        colorbar;
    elseif (strcmpi(plotType, 'peakLat'))
        scatter(set1.score(:,1), set2.score(:,1), [], ...
            abs(boutPeakVel.lat * 1000), 'filled', ...
            'MarkerFaceAlpha', 0.8, 'MarkerEdgeAlpha', 0.3);
        colormap(gca, 'cool');
        colorbar;
    elseif (strcmpi(plotType, 'peakSlide'))
        scatter(set1.score(:,1), set2.score(:,1), [], ...
            abs(boutPeakVel.slide), 'filled', ...
            'MarkerFaceAlpha', 0.8, 'MarkerEdgeAlpha', 0.3);
        colormap(gca, 'cool');
        colorbar;
    else
        scatter(set1.score(:,1), set2.score(:,1), 'filled', ...
            'MarkerFaceAlpha', 0.8, 'MarkerEdgeAlpha', 0.3);
    end

    hold on;

    % plot reference lines
    xLimits = xlim;
    yLimits = ylim;

    minPt = min([xLimits(1) yLimits(1)]);
    maxPt = max([xLimits(2) yLimits(2)]);

    if (xyLine)
        plot([minPt maxPt], [minPt maxPt], 'k', 'LineWidth', 2);
    end

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