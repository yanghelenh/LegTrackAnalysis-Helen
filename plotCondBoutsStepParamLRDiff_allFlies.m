% plotCondBoutsStepParamLRDiff_allFlies.m
%
% Function to plot 1D step parameters for bouts aligned to yaw velocity
%  peaks. From output of saveLegStepParamCond_bouts()
% One figure with 6 subplots, 1 for each leg
% Can select more than 1 output file, will plot on same graph
% Has option to plot reference steps (from saveLegStepParamCond_indpt()).
%  Select one file through GUI. Plots mean and std/SEM (as specified by
%  input)
% Select output files through GUI
% Note: make sure all output files selected have same maxNumSteps. Plotting
%  will get confused otherwise
%
% INPUTS:
%   whichParam - name of parameter to plot
%   swingOrStance - boolean for whether to plot param during swing (false)
%       or stance (true)
%   datDir - directory with output files
%   yScale - y scale for plots, as [min max]
%   semError - boolean for whether to plot SEM; if false, plots std dev
%
% OUTPUTS:
%   none, but generates plot
%
% CREATED: 5/15/23 - HHY
%
% UPDATED:
%   5/15/23 - HHY
%   6/5/23 - HHY add reference steps
%   6/6/23 - HHY - for reference file, make sure stepDirections mean and
%       std use circular statistics
%   6/7/23 - HHY - remove outliers from reference file
%
function plotCondBoutsStepParamLRDiff_allFlies(whichParam, swingOrStance, ...
    datDir, yScale, semError)

    % all step parameters that are circular variables - need to use
    %  circular stats - 6/6/23 - HHY
    circStepParams = {'stepDirections'};

    % prompt user to select output files from saveLegStepParamCond_bouts()
    disp('Select output files from saveLegStepParamCond_bouts()');
    [outputFNames, outputPath] = uigetfile('*.mat', ...
        'Select cond_bouts files', datDir, 'MultiSelect', 'on');

    % if only 1 file selected, not cell array; make sure loop still
    %  works 
    % num flies is number of files
    if (iscell(outputFNames))
        numFiles = length(outputFNames);
    else
        numFiles = 1;
    end

    % preallocate
    allMeans = cell(numFiles, 1);
    allErrs = cell(numFiles, 1);
    numBouts = cell(numFiles, 1);
    legendStr = cell(numFiles, 1);

    % loop through all files, save output matrices for each file
    % save also, cond, numBouts
    for i = 1:numFiles
        % handle whether it's a cell array or not
        if (iscell(outputFNames))
            outName = outputFNames{i};
        else
            outName = outputFNames;
        end
        
        outputFullPath = [outputPath outName];

        % load data - with appropriate error and phase
        if (semError && swingOrStance) % SEM and Stance
            load(outputFullPath, 'stanceParamMeans', 'stanceParamSEM', ...
                'stanceParamN', 'maxNumSteps', 'cond');

            % if circular parameter and stance, wrap to 360 for plotting
            if(any(strcmpi(whichParam, circStepParams)))
                allMeans{i} = wrapTo360(stanceParamMeans.(whichParam));
            else % not circular
                allMeans{i} = stanceParamMeans.(whichParam);
            end

            allErrs{i} = stanceParamSEM.(whichParam);
            numBouts{i} = stanceParamN.(whichParam);
        elseif (semError && ~swingOrStance) % SEM and Swing
            load(outputFullPath, 'swingParamMeans', 'swingParamSEM', ...
                'swingParamN', 'maxNumSteps', 'cond');

            allMeans{i} = swingParamMeans.(whichParam);
            allErrs{i} = swingParamSEM.(whichParam);
            numBouts{i} = swingParamN.(whichParam);
        elseif (~semError && swingOrStance) % std dev and stance
            load(outputFullPath, 'stanceParamMeans', 'stanceParamStd', ...
                'stanceParamN', 'maxNumSteps', 'cond');

            % if circular parameter and stance, wrap to 360 for plotting
            if(any(strcmpi(whichParam, circStepParams)))
                allMeans{i} = wrapTo360(stanceParamMeans.(whichParam));
            else % not circular
                allMeans{i} = stanceParamMeans.(whichParam);
            end

            allErrs{i} = stanceParamStd.(whichParam);
            numBouts{i} = stanceParamN.(whichParam);
        elseif (~semError && ~swingOrStance) % std dev and swing
            load(outputFullPath, 'swingParamMeans', 'swingParamStd', ...
                'swingParamN', 'maxNumSteps', 'cond');

            allMeans{i} = swingParamMeans.(whichParam);
            allErrs{i} = swingParamStd.(whichParam);
            numBouts{i} = swingParamN.(whichParam);
        end

        % get legend str, named based on condition
        thisLegendStr = [];
        for j = 1:length(cond.whichParam)
            thisLegendStr = [thisLegendStr cond.whichParam{j}...
                cond.cond{j}];
            if (j~=length(cond.whichParam))
                thisLegendStr = [thisLegendStr ', '];
            end
        end
        legendStr{i} = thisLegendStr;
    end



    % get step time points (for x-axis): 0 for step at peak
    stepTPts = [fliplr((1:maxNumSteps) * -1), 0, 1:maxNumSteps];

    % initialize figure
    f = figure;
    c = colormap('lines');

    for i = 1:3
        subplot(1,3,i);
        hold on;

        % loop through all conditions
        for j = 1:numFiles
            % plot for this condition
            errorbar(stepTPts', allMeans{j}(:,i), allErrs{j}(:,i), ...
                'Marker', 'x','LineWidth',1, 'CapSize', 0, 'Color', c(j,:));

            hold on;

            % add n to legendStr
            thisLegendStr = [legendStr{j} '; n = ' ...
                num2str(numBouts{j}(:,i)')];

            % legend depends on if there's reference
            legendStrLegs{j,i} = thisLegendStr;
        end

        % axis scale and label
        ylim(yScale);
        xScale = xlim;
        xScale(1) = xScale(1) - (0.1 * (stepTPts(end)-stepTPts(1)));
        xScale(2) = xScale(2) + (0.1 * (stepTPts(end)-stepTPts(1)));
        xlim(xScale);

        line(xScale, [0 0], 'Color','k', 'LineWidth', 1);

        xticks(stepTPts);

        xlabel('Steps');
        ylabel(whichParam);

        % legend
        legend(legendStrLegs{:,i});
    end

    sgtitle(whichParam);
end