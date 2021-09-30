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
%       movAvgWinLen - length of window, in frames, for moving average
%       maxminWinLen - length of window, in frames, for finding max/min
%       adjThresh - threshold for adjacent indicies
%       maxPosVelThresh - threshold of what 1st derivative (vel) should 
%           exceed in the positive direction for a position maximum
%       maxNegVelThresh - in negative direction
%       minPosVelThresh - in pos dir, for leg pos min
%       minNegVelThresh - in neg dir, for leg pos min
%       numNegVelFrames - num of frames before and after max/min to check 
%           for negative velocity thresh
%       numPosVelFrames - num frames to check for positive velocity thresh
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
%   9/29/21 - HHY
%
function [maxIndsAll, minIndsAll, maxWhichLeg, minWhichLeg, userSelVal] = ...
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
    
    % button parameters
    btnXPosStart = 100;
    btnXSpace = 150;
    btnYPos = 100;
    btnHeight = 20;
    btnWidth = 80;
    
    % ranges for all legRevParams
    lrpRanges.movAvgWinLen = [1 100];
    lrpRanges.maxminWinLen = [1 100];
    lrpRanges.adjThresh = [1 15];
    lrpRanges.maxPosVelThresh = [0 1e-3];
    lrpRanges.maxNegVelThresh = [-1e-2 0];
    lrpRanges.minPosVelThresh = [0 1e-3];
    lrpRanges.minNegVelThresh = [-1e-2 0];
    lrpRanges.numNegVelFrames = [1 50];
    lrpRanges.numPosVelFrames = [1 50];
    
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
    
    
    
    % loop through all legs
    for i = 1%:length(legIDs.ind)
        
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
        [maxInds, minInds] = findLegReversals(...
            legTrack.srnLegX(:,thisLegInd), ...
            legTrack.legXVel(:,thisLegInd), moveNotMove.notMoveInd, ...
            legRevParamsInit);
        
        % get leg position values for these indicies
        initMaxVals = legTrack.srnLegX(maxInds, thisLegInd);
        initMinVals = legTrack.srnLegX(minInds, thisLegInd);
        
        % get shading for not moving
        notMovingX = [moveNotMove.notMoveBout(:,1)'; ...
            moveNotMove.notMoveBout(:,1)'; moveNotMove.notMoveBout(:,2)'; ...
            moveNotMove.notMoveBout(:,2)'];
        notMovingXT = legTrack.t(notMovingX);
        
        % get ylimits for this leg
        yMinLim = min(legTrack.srnLegX(:,thisLegInd)) * 1.1;
        yMaxLim = max(legTrack.srnLegX(:,thisLegInd)) * 1.1;
        
        y0 = ones(1,size(moveNotMove.notMoveBout,1)) * yMinLim;
        y1 = ones(1,size(moveNotMove.notMoveBout,1)) * yMaxLim;
        notMovingY = [y0; y1; y1; y0];
        
        
        % initialize figure
        f = figure('Position', [20 20 1600 920]);
        
        % plot this leg
        legPosAx = subplot('Position', [0.05 0.3 0.6 0.5]);
        % leg position
        plot(legTrack.t, legTrack.srnLegX(:,thisLegInd));
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
            allButton{j} = uicontrol(f, 'Style', 'togglebutton');
            
            thisButtonXPos = btnXPosStart + btnXSpace * (j-1);
            
            allButton{j}.Position = ...
                [thisButtonXPos btnYPos btnWidth btnHeight];
            
            allButton{j}.String = buttonNames{j};
            allButton{j}.Callback = {@manualSelect, j};  
        end
        
        
        
        
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
        [maxInds, minInds] = findLegReversals(...
            legTrack.srnLegX(:,thisLegInd), ...
            legTrack.legXVel(:,thisLegInd), moveNotMove.notMoveInd, ... 
            thisLegRevParams);

        % get leg position values for these indicies
        maxVals = legTrack.srnLegX(maxInds, thisLegInd);
        minVals = legTrack.srnLegX(minInds, thisLegInd);

        % plot leg position, max and min
        cla(legPosAx);

        % leg position
        plot(legPosAx,legTrack.t, legTrack.srnLegX(:,thisLegInd));
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
        thisButton = buttonNames{nameInd}; % as string
        
        % initialize logical for stopping manual selection
        selDone = 0;
        
        % initialize vector tracking indicies
        theseInds = [];
        
        % ginput to get coordinates, loop until user presses keyboard key
        while ~selDone
            [thisX, thisY, thisButton] = ginput(1);
            
            % if user presses keyboard key (button ascii >=32), end loop
            if (thisButton >= 32)
                selDone = 1;
                
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
                
                % unselect button
                allButton{nameInd}.Value = allButton{nameInd}.Min;
                
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
                switch thisButton
                    case 'Delete Max'
                        [thisInd, maxInds] = addRmvIndAtSelX(...
                            maxInds, legTrack.t, thisX, 'remove', 'max');
                    case 'Add Max'
                        [thisInd, maxInds] = addRmvIndAtSelX(...
                            maxInds, legTrack.t, thisX, 'add', 'max');
                    case 'Delete Min'
                        [thisInd, minInds] = addRmvIndAtSelX(...
                            minInds, legTrack.t, thisX, 'remove', 'min');
                    case 'Add Min'
                        [thisInd, minInds] = addRmvIndAtSelX(...
                            minInds, legTrack.t, thisX, 'add', 'min');
                end
                
                % update tracking of indices added/removed
                theseInd = [theseInd; thisInd];
            end
            
        end
        
        
        
    end

    % function for when next leg/done button is pressed
    function donePushed(src, event)
        % toggle logical to done, stops loop
        thisLegDone = 1;
    end


end
