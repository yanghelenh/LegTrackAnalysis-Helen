% plotPairwiseCorr_heatmap.m
%
% Given one output file from saveStepParamPariwiseCorr_bouts(), plots the
%  actual correlation as well as the bootstrapped, randomized one as
%  heatmaps
%
% INPUTS:
%   datDir - path to folder containing correlation files
%   cScale - scale for color axis
%
% OUTPUTS:
%   none, but generates plot
%
% CREATED: 7/26/23 - HHY
%
% UPDATED:
%   7/26/23 - HHY
%
function plotPairwiseCorr_heatmap(datDir, cScale)

    % prompt user to select cond_bout() file
    [condBoutFName, condBoutPath] = uigetfile('*.mat', ...
        'Select output file', datDir, 'MultiSelect', 'off');

    % load data from cond_bout file
    fullFilePath = [condBoutPath filesep condBoutFName];

    load(fullFilePath, 'allVarNames', 'corrMatrix', 'corrMatrixRand');

    % get mean randomized corr matrix
    meanRandMatrix = mean(corrMatrixRand,3);


    % plot heat maps
    cFig = figure;
    heatmap(allVarNames, allVarNames, corrMatrix, 'Title','Correlations',...
        'Colormap',turbo, 'CellLabelColor','none');
    caxis(cScale);

    mFig = figure;
    heatmap(allVarNames, allVarNames, meanRandMatrix, 'Title', ...
        'Mean Randomized Correlations', 'Colormap', turbo, ...
        'CellLabelColor','none');
    caxis(cScale);
end