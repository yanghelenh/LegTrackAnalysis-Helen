% interactGetLegReversals.m
%
% Function for getting max and min in leg positions, when leg reverses
%  direction
% Brings up GUI, 1 leg at a time, plotting leg X position, x's for computed
%  min and max. Starts with sliders to adjust parameters. When user clicks
%  button, changes to mode to select points to delete or add.
%
% INPUTS:
%   legTrack - struct of leg tracking data, output of preprocessLegTrack
%   moveNotMove - struct of moving/not-moving bout indicies, starts/ends
%   legRevParams - struct of parameters for getting these leg reversals
%       minProm - parameter MinPeakProminence of findpeaks
%       minDist - parameter MinPeakDistance of findpeaks
%   legIDs - struct of parameters, IDs legs
%       ind - indicies of legs, into raw position matricies, carried
%           throughout
%       names - names of each leg, matching with ind
%
% OUTPUTS:
%   maxIndsAll - indicies (frames) of leg position maxima, all legs in 1
%       vector
%   minIndsAll - indicies of leg position minima
%   maxWhichLeg - number indicating which leg max indicies belong to, same
%       size & matched with maxIndsAll
%   minWhichLeg - which legs for min indicies
%   legRevParams - updated values
%   userSelVal - struct array (1 element for each leg) with user selected
%       values used for determining max and min indicies
%       legRevParams - updated parameter values
%       maxAdded - indices of maxima added
%       maxRmved - indices of maxima removed
%       minAdded - indices of minima added
%       minRmved - indices of minima removed
%
% CREATED: 9/29/21 - HHY
%
% UPDATED:
%   10/1/21 - HHY
%   7/7/22 - HHY - moveNotMove now contains both leg and FicTrac data,
%       update calls to reflect new field names
%   6/21/23 - HHY - update to use getLegReversals_1Leg() instead of
%       findLegReversals()
%   6/23/23 - HHY - bug fix; srnf instead of srn
%
function [maxIndsAll, minIndsAll, maxWhichLeg, minWhichLeg, ...
    legRevParams, userSelVal] = ...
    interactGetLegReversals(legTrack, moveNotMove, legRevParams, legIDs)

    % some parameters
    tRange = 10; % sec, amount of data to display
    xMax = legTrack.t(end); % max value
    % names of manual selection buttons
    buttonNames = {'Delete Max', 'Add Max', 'Delete Min', 'Add Min'};
    
    % slider parameters
    sldXPos = 1200;
    sldYPosStart = 850;
    sldYPosEnd = 50;
    sldHeight = 20;
    sldWidth = 300;
    numSld = length(fieldnames(legRevParams));
    sldYSpace = round((sldYPosStart - sldYPosEnd) / numSld);
    sldTMajorStep = (tRange/xMax) * 0.9;
    sldTMinorStep = 0.1 * sldTMajorStep;
    
    % button parameters
    btnXPosStart = 100;
    btnXSpace = 150;
    btnYPos = 120;
    btnHeight = 20;
    btnWidth = 80;
    
    % 6/21/23 - these are obsolete
%     % ranges for all legRevParams
%     lrpRanges.movAvgWinLen = [1 100];
%     lrpRanges.maxminWinLen = [1 100];
%     lrpRanges.adjThresh = [1 15];
%     lrpRanges.maxPosVelThresh = [0 1e-3];
%     lrpRanges.maxNegVelThresh = [-1e-2 0];
%     lrpRanges.minPosVelThresh = [0 1e-3];
%     lrpRanges.minNegVelThresh = [-1e-2 0];
%     lrpRanges.numNegVelFrames = [1 50];
%     lrpRanges.numPosVelFrames = [1 50];
    
    % ranges for legRevParams
    lrpRanges.minProm = [0 0.5];
    lrpRanges.minDist = [1 50];
    
    % preallocate userSelVal array of structs
    % single struct
    userSelVal1Ele.legRevParams = legRevParams;
    userSelVal1Ele.maxAdded = [];
    userSelVal1Ele.maxRmved = [];
    userSelVal1Ele.minAdded = [];
    userSelVal1Ele.minRmved = [];
    % repmat single struct to create array of structs, 1 for each leg
    userSelVal = repmat(userSelVal1Ele,1,length(legIDs.ind));
    
    % names for all parameters in legRevParams
    lrpNames = fieldnames(legRevParams);
    
    % initialize maxIndsAll, minIndsAll, maxWhichLeg, minWhichLeg
    maxIndsAll = [];
    minIndsAll = [];
    maxWhichLeg = [];
    minWhichLeg = [];
    
    
    
    % loop through all legs
    for i = 1:length(legIDs.ind)
        
        global selDone
        
        % initialize logical for next leg/done button press
        thisLegDone = 0;
        
        % initialize vectors for tracking max/min added/removed
        maxAdded = [];
        maxRmved = [];
        minAdded = [];
        minRmved = [];
        
        % index and name for this leg
        thisLegInd = legIDs.ind(i);
        thisLegName = legIDs.names{i};
        
        % save initial parameters
        legRevParamsInit = legRevParams;
        
        % leg reversal params struct that changes for this leg
        thisLegRevParams = legRevParams;
        
        
        % get max and min for this leg, with initial parameters
        % 6/21/23 - replace findLegReversals with getLegReversals_1Leg()
%         [maxInds, minInds] = findLegReversals(...
%             legTrack.srnLegX(:,thisLegInd), ...
%             legTrack.legXVel(:,thisLegInd), moveNotMove.legNotMoveInd, ...
%             legRevParamsInit);

        [maxInds, minInds] = getLegReversals_1Leg(...
            legTrack.srnfLegX(:,thisLegInd), ...
            moveNotMove.legNotMoveInd, legRevParamsInit);
        
        % get leg position values for these indicies
        initMaxVals = legTrack.srnfLegX(maxInds, thisLegInd);
        initMinVals = legTrack.srnfLegX(minInds, thisLegInd);
        
        % get shading for not moving
        notMovingX = [moveNotMove.legNotMoveBout(:,1)'; ...
            moveNotMove.legNotMoveBout(:,1)'; moveNotMove.legNotMoveBout(:,2)'; ...
            moveNotMove.legNotMoveBout(:,2)'];
        notMovingXT = legTrack.t(notMovingX);
        
        % get ylimits for this leg
        yMinLim = min(legTrack.srnfLegX(:,thisLegInd)) * 1.1;
        yMaxLim = max(legTrack.srnfLegX(:,thisLegInd)) * 1.1;
        
        y0 = ones(1,size(moveNotMove.legNotMoveBout,1)) * yMinLim;
        y1 = ones(1,size(moveNotMove.legNotMoveBout,1)) * yMaxLim;
        notMovingY = [y0; y1; y1; y0];
        
        
        % initialize figure
        f = figure('Position', [20 20 1600 920]);
        
        % plot this leg
        legPosAx = subplot('Position', [0.05 0.3 0.6 0.5]);
        % leg position
        plot(legTrack.t, legTrack.srnfLegX(:,thisLegInd));
        hold on;
        
        % plot max - x's 
        plot(legTrack.t(maxInds), initMaxVals, 'x', ...
            'LineStyle', 'none');
        % plot min - o's
        plot(legTrack.t(minInds), initMinVals, 'o', ...
            'LineStyle', 'none');
        
        % plot shading for not moving bouts
        patch(notMovingXT, notMovingY, 'black', 'FaceAlpha', 0.3);
        
        xlim([0 tRange]);
        ylim([yMinLim yMaxLim]);
        title(sprintf('Leg %s Position, Max and Min', thisLegName));
        xlabel('Time (s)');
        
        
        % text for names of parameter sliders
        % invisible axes for text
        txtAx = axes('Position', [0 0 1 1], 'Visible', 'off');
        set(gcf, 'CurrentAxes', txtAx);
        
        % text for slider labels
        for j = 1:length(lrpNames)
            thisTxtPos = sldYPosStart - sldYSpace * (j-1) + 10;
            
            text(sldXPos - 120, thisTxtPos, lrpNames{j}, ...
                'Units', 'pixels', 'FontSize', 12);
        end
        
        % text for values of parameter sliders
        allTxtH = {}; % handles to text obj
        for j = 1:length(lrpNames)
            thisTxtPos = sldYPosStart - sldYSpace * (j-1) + 10;
            
            thisDispVal = num2str(legRevParams.(lrpNames{j}));
            
            allTxtH{j} = text(sldXPos + 320, thisTxtPos, thisDispVal, ...
                'Units', 'pixels', 'FontSize', 12);
        end
        
        % slider adjusting view of leg positions
        tSlider = uicontrol(f, 'Style', 'slider', 'Position', ...
            [100 200 600 20]);
        tSlider.Value = 0;
        tSlider.SliderStep = [sldTMinorStep sldTMajorStep];
        tSlider.Callback = @updateTLim;
         

            
        % sliders for adjusting each of the parameter values
        % initialize cell array for slider objects
        allSld = {};
        % loop through all parameters, creating sliders
        for j = 1:length(lrpNames)
            allSld{j} = uicontrol(f, 'Style', 'slider');

            thisSldYPos = sldYPosStart - sldYSpace * (j - 1);

            allSld{j}.Position = [sldXPos thisSldYPos sldWidth sldHeight];

            allSld{j}.Value = thisLegRevParams.(lrpNames{j});
            thisParamRange = lrpRanges.(lrpNames{j});
            allSld{j}.Max = thisParamRange(2);
            allSld{j}.Min = thisParamRange(1);
            allSld{j}.Callback = {@updateGraph, j};  
        end
        
        % buttons to switch to manual selection: delete max, add max,
        %  delete min, add min
        % initialize cell array for button objects
        allButton = {};
        % loop through all buttons
        for j = 1:length(buttonNames)
            allButton{j} = uicontrol(f, 'Style', 'pushbutton');
            
            thisButtonXPos = btnXPosStart + btnXSpace * (j-1);
            
            allButton{j}.Position = ...
                [thisButtonXPos btnYPos btnWidth btnHeight];
            
            allButton{j}.String = buttonNames{j};
            allButton{j}.Callback = {@manualSelect, j};  
        end
        
        % text: to display which manual selection button pressed
        msTxt = uicontrol(f, 'Style', 'text');
        msTxt.Position = [thisButtonXPos + 200, btnYPos, 200 20];
        msTxt.FontSize = 12;
        msTxt.String = 'Automated Detection, use sliders';
        
        % button for next leg or to end (on last leg)
        % get string for text on this button
        if (i == length(legIDs.ind)) %for last leg
            legDoneStr = 'Done';
        else % for intermediate leg
            legDoneStr = 'Next Leg';
        end
        
        legDoneButton = uicontrol(f, 'Style', 'pushbutton', 'String', ...
            legDoneStr, 'Position', [100 30 50 20]);
        legDoneButton.Callback = @donePushed;
    

        
        % loop until user hits next leg/done button
        while ~thisLegDone
            pause(0.1);
        end
        
        % update all return values for this leg

        % return values for this, concatenate into vector for all legs
        maxIndsAll = [maxIndsAll; maxInds];
        minIndsAll = [minIndsAll; minInds];
        
        % which leg marker, generate vectors
        thisMaxWhichLeg = ones(size(maxInds)) * legIDs.ind(i);
        thisMinWhichLeg = ones(size(minInds)) * legIDs.ind(i);
        % concatenate into vector for all legs
        maxWhichLeg = [maxWhichLeg; thisMaxWhichLeg];
        minWhichLeg = [minWhichLeg; thisMinWhichLeg];
        
        % update userSelVal
        userSelVal(i).legRevParams = thisLegRevParams;
        userSelVal(i).maxAdded = maxAdded;
        userSelVal(i).maxRmved = maxRmved;
        userSelVal(i).minAdded = minAdded;
        userSelVal(i).minRmved = minRmved;
        
        % close this figure
        close(f);
        
    end
    
    
    
    % functions for figure, down here b/c can't be in loop
    
    % function to update display region for plot, every time that
    %  slider is moved
    function updateTLim(src, event)
        xlim(legPosAx, ...
            [tSlider.Value * (xMax-tRange),...
            tSlider.Value * (xMax-tRange) + tRange]);
    end

    % function for updating the figure, thisLegRevParams, max and min
    %  indices every time a slider is moved
    function updateGraph(src, event, nameInd)

        % for those parameters that must be integers, round slider value
        switch lrpNames{nameInd}
            case {'movAvgWinLen', 'maxminWinLen', 'adjThresh', ...
                    'numNegVelFrames', 'numPosVelFrames'}
                thisVal = round(allSld{nameInd}.Value);
            otherwise
                thisVal = allSld{nameInd}.Value;
        end

        % change the apppropriate parameter value in thisLegRevParams
        thisLegRevParams.(lrpNames{nameInd}) = thisVal;

        % get new max and min calls
        % 6/21/23 - findLegReversals() replaced by getLegReversals_1Leg()
%         [maxInds, minInds] = findLegReversals(...
%             legTrack.srnLegX(:,thisLegInd), ...
%             legTrack.legXVel(:,thisLegInd), moveNotMove.legNotMoveInd, ... 
%             thisLegRevParams);

        [maxInds, minInds] = getLegReversals_1Leg(...
            legTrack.srnfLegX(:,thisLegInd), ...
            moveNotMove.legNotMoveInd, thisLegRevParams);

        % get leg position values for these indicies
        maxVals = legTrack.srnfLegX(maxInds, thisLegInd);
        minVals = legTrack.srnfLegX(minInds, thisLegInd);

        % plot leg position, max and min
        cla(legPosAx);

        % leg position
        plot(legPosAx,legTrack.t, legTrack.srnfLegX(:,thisLegInd));
        hold on;

        % plot max - x's 
        plot(legPosAx,legTrack.t(maxInds), maxVals, 'x', ...
            'LineStyle', 'none');
        % plot min - o's
        plot(legPosAx,legTrack.t(minInds), minVals, 'o', ...
            'LineStyle', 'none');

        % plot shading for not moving bouts
        patch(legPosAx, notMovingXT, notMovingY, 'black', 'FaceAlpha', 0.3);

        % update display value around slider
        allTxtH{nameInd}.String = num2str(thisVal);
    end

    % function for when one of the manual selection buttons is pressed
    function manualSelect(src, event, nameInd)
        selDone = 1;
        thisButtonName = buttonNames{nameInd}; % as string
        
        % change text to display which button was pressed
        msTxt.String = thisButtonName;
        
        % clear parameter sliders, text if they're still present
        if (isvalid(allSld{1}))
            for k = 1:length(allSld)
                delete(allSld{k});
            end
            
            % delete all text
            cla(txtAx);
        end
        
        % initialize logical for stopping manual selection
        selDone = 0;
        
        % initialize vector tracking indicies
        theseInd = [];
        
        % ginput to get coordinates, loop until user presses keyboard key
        while ~selDone
            [thisX, thisY, thisButton] = ginput(1);
            
            % if user presses keyboard key (button ascii >=32), end loop
            if (thisButton >= 32)
                selDone = 1;
                
                % update text 
                msTxt.String = 'No manual selection active';
                
                % update vectors of max/min added/removed
                switch thisButton
                    case 'Delete Max'
                        maxRmved = [maxRmved; theseInd];
                    case 'Add Max'
                        maxAdded = [maxAdded; theseInd];
                    case 'Delete Min'
                        minRmved = [minRmved; theseInd];
                    case 'Add Min'
                        minAdded = [minAdded; theseInd];
                end
                
                continue;
            end
            
            % check if thisX or thisY out of bounds (click for interaction
            %  with other part of figure)
            
            % current X min and max
            curXMin = tSlider.Value * (xMax-tRange);
            curXMax = tSlider.Value * (xMax-tRange) + tRange;
            
            % only add/remove max/min when selected point within current
            %  axis bounds
            if ((thisX >= curXMin) && (thisX <= curXMax) && ...
                    (thisY >= yMinLim) && (thisY <= yMaxLim))
                % call to add/remove max/min, depending on which button
                %  selected
                switch thisButtonName
                    case 'Delete Max'
                        [thisInd, maxInds] = rmvIndAtSelX(maxInds, ...
                            legTrack.t, thisX);
                    case 'Add Max'
                        [thisInd, maxInds] = addIndAtSelX(maxInds, ...
                            legTrack.t, legTrack.srnfLegX(:,thisLegInd), ...
                            thisX, 'max');
                    case 'Delete Min'
                        [thisInd, minInds] = rmvIndAtSelX(minInds, ...
                            legTrack.t, thisX);
                    case 'Add Min'
                        [thisInd, minInds] = addIndAtSelX(minInds, ...
                            legTrack.t, legTrack.srnfLegX(:,thisLegInd), ...
                            thisX, 'min'); 
                end
                
                % update tracking of indices added/removed
                theseInd = [theseInd; thisInd];
                
                % update max and min vals
                maxVals = legTrack.srnfLegX(maxInds, thisLegInd);
                minVals = legTrack.srnfLegX(minInds, thisLegInd);
            
                % update plot to add/remove this max/min        
                cla(legPosAx);

                % leg position
                plot(legPosAx,legTrack.t, legTrack.srnfLegX(:,thisLegInd));
                hold on;

                % plot max - x's 
                plot(legPosAx,legTrack.t(maxInds), maxVals, 'x', ...
                    'LineStyle', 'none');
                % plot min - o's
                plot(legPosAx,legTrack.t(minInds), minVals, 'o', ...
                    'LineStyle', 'none');

                % plot shading for not moving bouts
                patch(legPosAx, notMovingXT, notMovingY, 'black', ...
                    'FaceAlpha', 0.3);
            end
        end 
    end

    % function for when next leg/done button is pressed
    function donePushed(src, event)
        % toggle logical to done, stops loop
        thisLegDone = 1;
    end


end
