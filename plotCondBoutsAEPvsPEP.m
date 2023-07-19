% plotCondBoutsAEPvsPEP.m
%
% Function that takes output


function plotCondBoutsAEPvsPEP(datDir, xyScale, ref)

    % legs to subplot indices
    % puts left legs on left, and front legs on top
    subInd = [2 4 6 1 3 5]; 

    % prompt user to select output files from saveLegStepParamCond_bouts()
    disp('Select output files from saveLegStepParamCond_bouts()');
    [outputFNames, outputPath] = uigetfile('*.mat', 'Select AEP PEP files', ...
        datDir, 'MultiSelect', 'on');

    % if only 1 file selected, not cell array; make sure loop still
    %  works 
    % num flies is number of files
    if (iscell(outputFNames))
        numFiles = length(outputFNames);
    else
        numFiles = 1;
    end

    % preallocate
    allAEPXs = cell(numFiles, 1);
    allAEPYs = cell(numFiles, 1);
    allPEPXs = cell(numFiles, 1);
    allPEPYs = cell(numFiles, 1);
    legendStr = cell(numFiles, 1);

    % if reference
    if (ref)
        % prompt user to select reference file from 
        %  saveLegStepParamCond_indpt()
        disp('Select reference file');
        [refOutputName, refOutputPath] = uigetfile('*.mat', ...
            'Select reference AEP PEP file', datDir, 'MultiSelect', 'off');

        % load reference
        refOutputFullPath = [refOutputPath filesep refOutputName];
        load(refOutputFullPath, 'selStanceParams', 'legIDs');

        % number of legs
        numLegs = length(legIDs.ind);

        % preallocate
        refAEPxMeans = zeros(numLegs, 1);
        refAEPyMeans = zeros(numLegs, 1);
        refPEPxMeans = zeros(numLegs, 1);
        refPEPyMeans = zeros(numLegs, 1);
        
        refAEPxErrs = zeros(numLegs, 1);
        refAEPyErrs = zeros(numLegs, 1);
        refPEPxErrs = zeros(numLegs, 1);
        refPEPyErrs = zeros(numLegs, 1);

        refN = zeros(1, numLegs);

        % loop through each leg
        for j = 1:numLegs
            thisLeg = legIDs.ind(j);
            thisLegLog = selStanceParams.stepWhichLeg == thisLeg;
    
            % compute means, remove outliers
            refAEPxMeans(j) = mean(rmoutliers(selStanceParams.stepAEPX(thisLegLog)));
            refAEPyMeans(j) = mean(rmoutliers(selStanceParams.stepAEPY(thisLegLog)));
            refPEPxMeans(j) = mean(rmoutliers(selStanceParams.stepPEPX(thisLegLog)));
            refPEPyMeans(j) = mean(rmoutliers(selStanceParams.stepPEPY(thisLegLog)));
    
            % n for this leg
            thisN = length(rmoutliers(selStanceParams.stepAEPX(thisLegLog)));
            refN(j) = thisN;

            % compute errors, remove outliers
            refAEPxErrs(j) = std(rmoutliers(selStanceParams.stepAEPX(thisLegLog))) ...
                / sqrt(thisN);
            refAEPyErrs(j) = std(rmoutliers(selStanceParams.stepAEPY(thisLegLog)))...
                / sqrt(thisN);
            refPEPxErrs(j) = std(rmoutliers(selStanceParams.stepPEPX(thisLegLog)))...
                / sqrt(thisN);
            refPEPyErrs(j) = std(rmoutliers(selStanceParams.stepPEPY(thisLegLog)))...
                / sqrt(thisN);
    
        end
    end

    % loop through all files, save AEP and PEP matrices for each file
    % save also, cond, numBouts
    for i = 1:numFiles
        % handle whether it's a cell array or not
        if (iscell(outputFNames))
            outName = outputFNames{i};
        else
            outName = outputFNames;
        end
        
        outputFullPath = [outputPath outName];


        % load data
        load(outputFullPath, 'selStanceParams', 'maxNumSteps', 'cond');

        % other variables
        allAEPXs{i} = selStanceParams.stepAEPX;
        allAEPYs{i} = selStanceParams.stepAEPY;
        allPEPXs{i} = selStanceParams.stepPEPX;
        allPEPYs{i} = selStanceParams.stepPEPY;

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
    numStepTPts = length(stepTPts);

    % loop through all steps
    for i = 1:numStepTPts

        % initialize figures
        xFig = figure;
        yFig = figure;
        c = colormap('lines');
    
        % preallocate legend str
        if(ref)
            legendStrX = cell(numFiles + 1,1);
            legendStrY = cell(numFiles + 1,1);
        else
            legendStrX = cell(numFiles,1);
            legendStrY = cell(numFiles,1);
        end

        % AEPX vs PEPX plot
        figure(xFig);

        for k = 1:6
            subplot(3,2,subInd(k));
            hold on;



        % loop through all files
        for j = 1:numFiles
            % plot
            scatter(squeeze(allAEPXs{j}(i,k,:)), squeeze(allPEPXs{j}(i,k,:)),...
                'MarkerFaceColor',c(j,:), 'MarkerFaceAlpha',0.3, ...
                'MarkerEdgeAlpha', 0);
    
            hold on;

            % add n to legendStr
            thisLegendStr = [legendStr{j}];

            if(ref)
                legendStrX{j} = thisLegendStr;
            else
                legendStrX{j} = thisLegendStr;
            end
        end

                % if reference, plot those AEPs in black
        if(ref)
            errorbar(refAEPxMeans(k), refPEPxMeans(k), ...
                refPEPxErrs(k), refPEPxErrs(k), ...
                refAEPxErrs(k), refAEPxErrs(k), ...
                'Marker','x', 'LineStyle','none','Color','black', ...
                'LineWidth',1.5);
        
            % legend
            thisLegendStr = ['Reference; n = ' num2str(refN)];
    
            legendStrX{end} = thisLegendStr;
    
            hold on;
        end

        end

        % get title string
        ttlStr = ['AEPX vs. PEPX, Step Num = ' num2str(stepTPts(i))];
        sgtitle(ttlStr);

        % x and y lims
%         xlim(xyScale);
%         ylim(xyScale);
%     
%         axis('equal');
    
        % x and y labels
%         xlabel('Body Lengths <-L - R->');
%         ylabel('Body Lengths <-P - A->')
    
        % reverse y axis (x values) so head (neg vals) is at top
%         set(gca, 'YDir','reverse');
    
        legend(legendStrX); 


        % PEP plot
        figure(yFig);

        for k = 1:6
            subplot(3,2,subInd(k));
            hold on;



        % loop through all files
        for j = 1:numFiles
            % plot

            scatter(squeeze(allAEPYs{j}(i,k,:)), squeeze(allPEPYs{j}(i,k,:)),...
                'MarkerFaceColor',c(j,:), 'MarkerFaceAlpha',0.3);
  
            hold on;

            % add n to legendStr
            thisLegendStr = [legendStr{j}];

            if(ref)
                legendStrY{j} = thisLegendStr;
            else
                legendStrY{j} = thisLegendStr;
            end
        end

                % if reference, plot those AEPs in black
        if(ref)
            errorbar(refAEPyMeans(k), refPEPyMeans(k), ...
                refPEPyErrs(k), refPEPyErrs(k), ...
                refAEPyErrs(k), refAEPyErrs(k), ...
                'Marker','x', 'LineStyle','none','Color','black', ...
                'LineWidth',1.5);
    
            % legend
            thisLegendStr = ['Reference; n = ' num2str(refN)];
    
            legendStrY{end} = thisLegendStr;
    
            hold on;
        end

        end

        % get title string
        ttlStr = ['AEPY vs PEPY, Step Num = ' num2str(stepTPts(i))];
        sgtitle(ttlStr);

%         % x and y lims
%         xlim(xyScale);
%         ylim(xyScale);
%     
%         axis('equal');
%     
%         % x and y labels
%         xlabel('Body Lengths <-L - R->');
%         ylabel('Body Lengths <-P - A->')
%     
%         % reverse y axis (x values) so head (neg vals) is at top
%         set(gca, 'YDir','reverse');
    
        legend(legendStrY); 
    end
end