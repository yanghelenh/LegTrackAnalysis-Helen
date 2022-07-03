% interactGetNotMovingIndWFt.m
%
% Function for selecting not-moving times. 
% Brings up GUI with plot of mid-leg positions and FicTrac total speed (as
%  plot over time and as histogram)
%  Shaded regions of not moving & red dots on not-moving (for leg position)
%  and sliders to adjust all not moving parameters that determine 
%  not-moving portions.
% Adaptation of interactGetNotMovingInd() and determineMoveThresh()
%
% INPUTS:
%   legTrack - struct of leg tracking data, output of preprocessLegTrack
%   fictracProc - struct of processed fictrac data, output of
%       filtFictrac_all
%   notMoveParams - struct of initial values of parameters
%       medFiltNumSamps - number of samples for median filtering leg
%           velocity
%       zeroXVelThreshMed - zero velocity threshold (on speed, actually),
%           for median filtered leg X velocity
%       zeroYVelThreshMed - zero velocity threshold (on speed, actually),
%           for median filtered leg Y velocity
%       zeroXVelThresh - zero velocity threshold on leg X velocity, w/o
%           filtering
%       zeroYVelThresh - zero velocity threshold on leg Y velocity, w/o
%           filtering
%       movePosXVelThresh - threshold for movement, positive X velocity
%       moveNegXVelThresh - threshold for movement, negative X velocity
%       movePosYVelThresh - threshold for movement, positive Y velocity
%       moveNegYVelThresh - threshold for movement, negative Y velocity
%       minBoutLen - minimum length of not-moving bout, in samples
%       stepNegXVelThresh - threshold for step, negative X velocity
%       maxTimeFromStep - maximum number of samples from step, for midpoint
%           of not-moving bout
%       adjBoutSepInit - if not-moving bouts are less than this many 
%           samples apart initially, merge them
%       adjBoutSepEnd - if not-moving bouts are less than this many 
%           samples apart at the end of the analysis, merge them
%       ftTotSpdThresh - threshold on smoothed total speed, normalized
%       ftMinBoutLen - minimum bout length of movement or stopping, in
%           seconds
%       cmbMethod - method for combining leg and FicTrac calls, string of
%           'union', 'intersect', 'legOnly', 'fictracOnly'
%       sigma - standard deviation of Gaussian kernel used to smooth
%           velocity, in seconds
%   r2LegInd - index of right mid-leg
%   l2LegInd - index of left mid-leg
%
% OUTPUTS:
%   notMoveInd - indicies for when fly not moving
%   notMoveBout - n x 2 matrix of not move bout start (col 1) and end
%       (col 2) indicies
%   moveInd - indicies for when fly moving (inverse of notMoveInd)
%   moveBout - n x 2 matrix of move bout start (col 1) and end (col 2)
%       indicies
%   notMoveParams - struct of final values of parameters
%
% CREATED: 6/30/22 - HHY
%
% UPDATED:
%   7/1/22 - HHY
%
function [legNotMoveInd, legNotMoveBout, legMoveInd, legMoveBout, ...
    ftNotMoveInd, ftNotMoveBout, ftMoveInd, ftMoveBout, notMoveParams] = ...
    interactGetNotMovingIndWFt(legTrack, fictracProc, notMoveParams, ...
    r2LegInd, l2LegInd)

    % list of notMoveParams that aren't slider options
    nonSliderParams = {'cmbMethod','sigma'};

    % list of combined methods options
    cmbMetOptions = {'intersect', 'union', 'legOnly', 'fictracOnly'};

    % some parameters
    tRange = 30; % sec, amount of data to display
    xMax = legTrack.t(end); % max value
    
    % slider parameters
    sldXPos = 1200;
    sldYPosStart = 850;
    sldYPosEnd = 200;
    sldHeight = 20;
    sldWidth = 300;
    numSld = length(fieldnames(notMoveParams)) - length(nonSliderParams);
    sldYSpace = round((sldYPosStart - sldYPosEnd) / numSld);
    
    % ranges for all notMoveParams
    nmpRanges.medFiltNumSamps = [1 100];
    nmpRanges.zeroXVelThreshMed = [0 5e-3];
    nmpRanges.zeroXVelThresh = [0 0.2];
    nmpRanges.movePosXVelThresh = [0 0.1];
    nmpRanges.moveNegXVelThresh = [-0.2 0];
    nmpRanges.zeroYVelThreshMed = [0 0.05];
    nmpRanges.zeroYVelThresh = [0 0.5];
    nmpRanges.movePosYVelThresh = [0 0.1];
    nmpRanges.moveNegYVelThresh = [-0.2 0];
    nmpRanges.minBoutLen = [1 100];
    nmpRanges.stepNegXVelThresh = [-0.5 0];
    nmpRanges.maxTimeFromStep = [0 1000];
    nmpRanges.adjBoutSepInit = [1 100];
    nmpRanges.adjBoutSepEnd = [1 100];
    
    % initialize logical for done button press
    userDone = 0;
    
    % remember inital parameters
    initNotMoveParams = notMoveParams;

    % smooth FicTrac total speed again
    % interframe interval for fictrac
    ifi = median(diff(fictracProc.t));
    ftSampRate = 1/ifi; % sample rate for fictrac

    % smoothing parameters
    sigmaSamp = round(notMoveParams.sigma * ftSampRate);
    padLen = 3 * sigmaSamp; % pad length, should be longer than sigma

    % for computing histogram, remove when
    %  FicTrac dropped as well as sigmaSamp from edges
    validInd = 1:length(fictracProc.totSpd);
    validInd(fictracProc.dropInd) = [];
    validInd = validInd(sigmaSamp:(end-sigmaSamp));

    % normalize total speed
    maxSpd = max(fictracProc.totSpd(validInd));
    
    totSpdNorm = fictracProc.totSpd ./ maxSpd;

    % Gaussian process smooth normalized total speed
    smoTotSpdNorm = gaussSmooth(totSpdNorm,padLen, sigmaSamp);

    % ranges for FicTrac not move parameters
    nmpRanges.ftTotSpdThresh = [0 1];
    nmpRanges.ftMinBoutLen = [0 10];

    % convert min bout length in sec to samples
    ftMinBoutLenSamp = notMoveParams.ftMinBoutLen * ftSampRate;


    % get not moving calls with initial leg parameters
    [legNotMoveInd, legNotMoveStartInd, legNotMoveEndInd] = findFlyNotMoving(...
        legTrack.legXVel, legTrack.legYVel, notMoveParams, r2LegInd, ...
        l2LegInd);

    % get not moving calls with initial FicTrac parameters
    [ftNotMoveInd, ftNotMoveStartInd, ftNotMoveEndInd] = ...
        findFlyNotMovingFt(smoTotSpdNorm, ...
        notMoveParams.ftTotSpdThresh, ftMinBoutLenSamp);

    % get not moving, combined b/w leg and FicTrac calls
    [legNotMoveInd, legNotMoveStartInd, legNotMoveEndInd, ...
        ftNotMoveInd, ftNotMoveStartInd, ftNotMoveEndInd] = ...
        findLegFtCmbNotMove(legNotMoveStartInd, legNotMoveEndInd, ...
        legTrack.t, ftNotMoveStartInd, ftNotMoveEndInd, fictracProc.t, ...
        notMoveParams.cmbMethod);
    
    % get shading for not moving, for leg
    legNotMovingX = [legNotMoveStartInd; legNotMoveStartInd; ...
        legNotMoveEndInd; legNotMoveEndInd];
    legNotMovingXT = legTrack.t(legNotMovingX);
    legY0 = ones(size(legNotMoveStartInd)) * -0.6;
    legY1 = ones(size(legNotMoveStartInd)) * 0.6;
    legNotMovingY = [legY0; legY1; legY1; legY0];

    % get shading for not moving, for FicTrac
    ftNotMovingX = [ftNotMoveStartInd; ftNotMoveStartInd; ...
        ftNotMoveEndInd; ftNotMoveEndInd];
    ftNotMovingXT = fictracProc.t(ftNotMovingX);
    ftY0 = zeros(size(ftNotMoveStartInd));
    ftY1 = ones(size(ftNotMoveStartInd));
    ftNotMovingY = [ftY0; ftY1; ftY1; ftY0];
    
    
    % initialize figure
    f = figure('Position', [20 20 1600 920]);

    
    % plot right leg
    r2Ax = subplot('Position', [0.05 0.75 0.6 0.2]);
    plot(legTrack.t, legTrack.srnLegX(:,r2LegInd));
    hold on;
    
    % plot dots for all not move ind
    plot(legTrack.t(legNotMoveInd), legTrack.srnLegX(legNotMoveInd,r2LegInd), ...
        '.','LineStyle','none');
    % plot shading for not moving bouts
    patch(legNotMovingXT, legNotMovingY, 'black', 'FaceAlpha', 0.3');
    
    xlim([0 tRange]);
    title('R2 X position');
    xlabel('Time (s)');
    
    
    % plot left leg
    l2Ax = subplot('Position', [0.05 0.45 0.6 0.2]);
    plot(legTrack.t, legTrack.srnLegX(:,l2LegInd));
    hold on;
    
    % plot dots for all not move ind
    plot(legTrack.t(legNotMoveInd), legTrack.srnLegX(legNotMoveInd,l2LegInd), ...
        '.','LineStyle','none');
    % plot shading for not moving bouts
    patch(legNotMovingXT, legNotMovingY, 'black', 'FaceAlpha', 0.3);
    
    xlim([0 tRange]);
    title('L2 X position');
    xlabel('Time (s)');

    % plot FicTrac smoothed total speed
    ftAx = subplot('Position', [0.05 0.15 0.6 0.2]);
    plot(fictracProc.t, smoTotSpdNorm);
    hold on;
    % plot shading for not moving bouts
    patch(ftNotMovingXT, ftNotMovingY, 'black','FaceAlpha',0.3);

    xlim([0 tRange]);
    title('FicTrac Smoothed, Normalized Total Speed');
    xlabel('Time(s)');
    
    % link axes for left and right legs
    linkaxes([r2Ax l2Ax], 'xy');
    
    
    % get all notMoveParam names
    nmpNames = fieldnames(notMoveParams);
    % get all slider names
    sldNames = nmpNames;
    % loop through all names to remove
    for i = 1:length(nonSliderParams)
        isNonSlider = strcmp(nonSliderParams{i},sldNames);
        sldNames(isNonSlider) = [];
    end

    % text for names of parameter sliders
    % invisible axes for text
    txtAx = axes('Position',[0 0 1 1], 'Visible', 'off');
    set(gcf, 'CurrentAxes', txtAx);
    for i = 1:length(sldNames)
        thisSldYPos = sldYPosStart - sldYSpace * (i - 1);
        
        text(sldXPos - 120, thisSldYPos + 10, sldNames{i}, ...
             'Units', 'pixels', 'FontSize', 12);
    end
    % text for values of parameter sliders
    allTxtH = {}; % handles to text obj
    for i = 1:length(sldNames)
        thisSldYPos = sldYPosStart - sldYSpace * (i - 1);
        
        thisDispVal = num2str(notMoveParams.(sldNames{i}));
        
        allTxtH{i} = text(sldXPos + 320, thisSldYPos + 10, thisDispVal, ...
             'Units', 'pixels', 'FontSize', 12);
    end
    
    % slider adjusting view of leg positions
    tSlider = uicontrol(f, 'Style', 'slider', 'Position', [100 50 600 20]);
    tSlider.Value = 0;
    tSlider.Callback = @updateTLim;
    
    % function to update display region for plot, every time that slider is
    %  moved
    function updateTLim(src, event)
        xlim(r2Ax, ...
            [tSlider.Value * (xMax-tRange),...
            tSlider.Value * (xMax-tRange) + tRange]);
        xlim(l2Ax, ...
            [tSlider.Value * (xMax-tRange),...
            tSlider.Value * (xMax-tRange) + tRange]);
        xlim(ftAx, ...
            [tSlider.Value * (xMax-tRange),...
            tSlider.Value * (xMax-tRange) + tRange]);
    end

    % radio buttons for selecting method of combining FicTrac and leg not
    %  move methods
    cmbRadioButton = uibuttongroup(f,'Visible','off','Position', ...
        [0.75 0.1 0.15 0.1], 'SelectionChangedFcn',@cmbButtonSelection);
    cmbButtons = {};
    cmbButtons{1} = uicontrol(cmbRadioButton, 'Style', 'radiobutton',...
        'String',cmbMetOptions{1},'Position',[10 60 100 30],...
        'HandleVisibility','off');
    cmbButtons{2} = uicontrol(cmbRadioButton, 'Style', 'radiobutton',...
        'String',cmbMetOptions{2},'Position',[130 60 100 30],...
        'HandleVisibility','off');
    cmbButtons{3} = uicontrol(cmbRadioButton, 'Style', 'radiobutton',...
        'String',cmbMetOptions{3},'Position',[10 10 100 30],...
        'HandleVisibility','off');
    cmbButtons{4} = uicontrol(cmbRadioButton, 'Style', 'radiobutton',...
        'String',cmbMetOptions{4},'Position',[130 10 100 30],...
        'HandleVisibility','off');
    % get current method, as index
    currCmbMethodInd = find(strcmp(notMoveParams.cmbMethod,cmbMetOptions));
    % set radio button value to current one
    cmbRadioButton.SelectedObject = cmbButtons{currCmbMethodInd};
    cmbRadioButton.Visible = 'on';
    

    % function for updating the figure, notMoveParams, notMove indices when
    %  new radio button selected
    function cmbButtonSelection(src, event)
        % get new selection value
        selMethod = event.NewValue.String;

        % change parameter value
        notMoveParams.cmbMethod = selMethod;

        % get not moving calls with initial leg parameters
        [legNotMoveInd, legNotMoveStartInd, legNotMoveEndInd] = findFlyNotMoving(...
            legTrack.legXVel, legTrack.legYVel, notMoveParams, r2LegInd, ...
            l2LegInd);
    
        % get not moving calls with initial FicTrac parameters
        [ftNotMoveInd, ftNotMoveStartInd, ftNotMoveEndInd] = ...
            findFlyNotMovingFt(smoTotSpdNorm, ...
            notMoveParams.ftTotSpdThresh, ftMinBoutLenSamp);

        % get not moving, combined b/w leg and FicTrac calls
        [legNotMoveInd, legNotMoveStartInd, legNotMoveEndInd, ...
            ftNotMoveInd, ftNotMoveStartInd, ftNotMoveEndInd] = ...
            findLegFtCmbNotMove(legNotMoveStartInd, legNotMoveEndInd, ...
            legTrack.t, ftNotMoveStartInd, ftNotMoveEndInd, fictracProc.t, ...
            notMoveParams.cmbMethod);
        
        % get shading for not moving, for leg
        legNotMovingX = [legNotMoveStartInd; legNotMoveStartInd; ...
            legNotMoveEndInd; legNotMoveEndInd];
        legNotMovingXT = legTrack.t(legNotMovingX);
        legY0 = ones(size(legNotMoveStartInd)) * -0.6;
        legY1 = ones(size(legNotMoveStartInd)) * 0.6;
        legNotMovingY = [legY0; legY1; legY1; legY0];
    
        % get shading for not moving, for FicTrac
        ftNotMovingX = [ftNotMoveStartInd; ftNotMoveStartInd; ...
            ftNotMoveEndInd; ftNotMoveEndInd];
        ftNotMovingXT = fictracProc.t(ftNotMovingX);
        ftY0 = zeros(size(ftNotMoveStartInd));
        ftY1 = ones(size(ftNotMoveStartInd));
        ftNotMovingY = [ftY0; ftY1; ftY1; ftY0];
        
        % plot right leg
        cla(r2Ax);
        plot(r2Ax,legTrack.t, legTrack.srnLegX(:,r2LegInd));

        % plot dots for all not move ind
        plot(r2Ax,legTrack.t(legNotMoveInd), legTrack.srnLegX(legNotMoveInd,r2LegInd), ...
            '.','LineStyle','none');
        % plot shading for not moving bouts
        patch(r2Ax, legNotMovingXT, legNotMovingY, 'black', 'FaceAlpha', 0.3');


        % plot left leg
        cla(l2Ax);
        plot(l2Ax, legTrack.t, legTrack.srnLegX(:,l2LegInd));

        % plot dots for all not move ind
        plot(l2Ax,legTrack.t(legNotMoveInd), legTrack.srnLegX(legNotMoveInd,l2LegInd), ...
            '.','LineStyle','none');
        % plot shading for not moving bouts
        patch(l2Ax,legNotMovingXT, legNotMovingY, 'black', 'FaceAlpha', 0.3');

        % plot FicTrac
        cla(ftAx);
        plot(ftAx, fictracProc.t, smoTotSpdNorm);
        hold on;
        % plot shading for not moving bouts
        patch(ftAx, ftNotMovingXT, ftNotMovingY, 'black','FaceAlpha',0.3);

    end
    
    % sliders for adjusting each of the parameter values
    % initialize cell array for slider objects
    allSld = {};
    % loop through all parameters, creating sliders
    for i = 1:length(sldNames)
        allSld{i} = uicontrol(f, 'Style', 'slider');
        
        thisSldYPos = sldYPosStart - sldYSpace * (i - 1);
        
        allSld{i}.Position = [sldXPos thisSldYPos sldWidth sldHeight];
        
        allSld{i}.Value = notMoveParams.(sldNames{i});
        thisParamRange = nmpRanges.(sldNames{i});
        allSld{i}.Max = thisParamRange(2);
        allSld{i}.Min = thisParamRange(1);
        allSld{i}.Callback = {@updateGraph, i};  
    end
    
    % function for updating the figure, notMoveParams, notMove indices
    %  every time a slider is moved
    function updateGraph(src, event, nameInd)
        
        % for those parameters that must be integers, round slider value
        switch sldNames{nameInd}
            case {'medFiltNumSamps', 'minBoutLen', 'maxTimeFromStep', ...
                    'adjBoutSepInit', 'adjBoutSepEnd'}
                thisVal = round(allSld{nameInd}.Value);
            otherwise
                thisVal = allSld{nameInd}.Value;
        end

        % change the apppropriate parameter value in notMoveParams
        notMoveParams.(sldNames{nameInd}) = thisVal;
        
        % get not moving calls with initial leg parameters
        [legNotMoveInd, legNotMoveStartInd, legNotMoveEndInd] = findFlyNotMoving(...
            legTrack.legXVel, legTrack.legYVel, notMoveParams, r2LegInd, ...
            l2LegInd);
    
        % get not moving calls with initial FicTrac parameters
        [ftNotMoveInd, ftNotMoveStartInd, ftNotMoveEndInd] = ...
            findFlyNotMovingFt(smoTotSpdNorm, ...
            notMoveParams.ftTotSpdThresh, ftMinBoutLenSamp);
    
        % get not moving, combined b/w leg and FicTrac calls
        [legNotMoveInd, legNotMoveStartInd, legNotMoveEndInd, ...
            ftNotMoveInd, ftNotMoveStartInd, ftNotMoveEndInd] = ...
            findLegFtCmbNotMove(legNotMoveStartInd, legNotMoveEndInd, ...
            legTrack.t, ftNotMoveStartInd, ftNotMoveEndInd, fictracProc.t, ...
            notMoveParams.cmbMethod);
        
        % get shading for not moving, for leg
        legNotMovingX = [legNotMoveStartInd; legNotMoveStartInd; ...
            legNotMoveEndInd; legNotMoveEndInd];
        legNotMovingXT = legTrack.t(legNotMovingX);
        legY0 = ones(size(legNotMoveStartInd)) * -0.6;
        legY1 = ones(size(legNotMoveStartInd)) * 0.6;
        legNotMovingY = [legY0; legY1; legY1; legY0];
    
        % get shading for not moving, for FicTrac
        ftNotMovingX = [ftNotMoveStartInd; ftNotMoveStartInd; ...
            ftNotMoveEndInd; ftNotMoveEndInd];
        ftNotMovingXT = fictracProc.t(ftNotMovingX);
        ftY0 = zeros(size(ftNotMoveStartInd));
        ftY1 = ones(size(ftNotMoveStartInd));
        ftNotMovingY = [ftY0; ftY1; ftY1; ftY0];
        
        % plot right leg
        cla(r2Ax);
        plot(r2Ax,legTrack.t, legTrack.srnLegX(:,r2LegInd));

        % plot dots for all not move ind
        plot(r2Ax,legTrack.t(legNotMoveInd), legTrack.srnLegX(legNotMoveInd,r2LegInd), ...
            '.','LineStyle','none');
        % plot shading for not moving bouts
        patch(r2Ax, legNotMovingXT, legNotMovingY, 'black', 'FaceAlpha', 0.3');


        % plot left leg
        cla(l2Ax);
        plot(l2Ax, legTrack.t, legTrack.srnLegX(:,l2LegInd));

        % plot dots for all not move ind
        plot(l2Ax,legTrack.t(legNotMoveInd), legTrack.srnLegX(legNotMoveInd,l2LegInd), ...
            '.','LineStyle','none');
        % plot shading for not moving bouts
        patch(l2Ax,legNotMovingXT, legNotMovingY, 'black', 'FaceAlpha', 0.3');

        % plot FicTrac
        cla(ftAx);
        plot(ftAx, fictracProc.t, smoTotSpdNorm);
        hold on;
        % plot shading for not moving bouts
        patch(ftAx, ftNotMovingXT, ftNotMovingY, 'black','FaceAlpha',0.3);
        
        % update display value around slider
        allTxtH{nameInd}.String = num2str(thisVal);
    end

    % button for user to reset parameter values to original ones
    resetButton = uicontrol(f, 'Style', 'pushbutton', 'String', 'Reset', ...
        'Position', [1100 50 50 20]);
    resetButton.Callback = @resetPushed;
    
    % function for resetting notMoveParams when reset button pushed
    function resetPushed(src, event)
        notMoveParams = initNotMoveParams;
        
        % get not moving calls with initial leg parameters
        [legNotMoveInd, legNotMoveStartInd, legNotMoveEndInd] = findFlyNotMoving(...
            legTrack.legXVel, legTrack.legYVel, notMoveParams, r2LegInd, ...
            l2LegInd);
    
        % get not moving calls with initial FicTrac parameters
        [ftNotMoveInd, ftNotMoveStartInd, ftNotMoveEndInd] = ...
            findFlyNotMovingFt(smoTotSpdNorm, ...
            notMoveParams.ftTotSpdThresh, ftMinBoutLenSamp);
    
        % get not moving, combined b/w leg and FicTrac calls
        [legNotMoveInd, legNotMoveStartInd, legNotMoveEndInd, ...
            ftNotMoveInd, ftNotMoveStartInd, ftNotMoveEndInd] = ...
            findLegFtCmbNotMove(legNotMoveStartInd, legNotMoveEndInd, ...
            legTrack.t, ftNotMoveStartInd, ftNotMoveEndInd, fictracProc.t, ...
            notMoveParams.cmbMethod);
        
        % update all slider values
        for j = 1:length(sldNames)
            allSld{j}.Value = notMoveParams.(sldNames{j});
            allTxtH{j}.String = num2str(notMoveParams.(sldNames{j}));
        end

        % update radio buttons
        metButInd = find(strcmp(notMoveParams.cmbMethod,cmbMetOptions));
        % set radio button value to current one
        cmbRadioButton.SelectedObject = cmbButtons{metButInd};
        
        % update plot
        % update patches
        % get shading for not moving, for leg
        legNotMovingX = [legNotMoveStartInd; legNotMoveStartInd; ...
            legNotMoveEndInd; legNotMoveEndInd];
        legNotMovingXT = legTrack.t(legNotMovingX);
        legY0 = ones(size(legNotMoveStartInd)) * -0.6;
        legY1 = ones(size(legNotMoveStartInd)) * 0.6;
        legNotMovingY = [legY0; legY1; legY1; legY0];
    
        % get shading for not moving, for FicTrac
        ftNotMovingX = [ftNotMoveStartInd; ftNotMoveStartInd; ...
            ftNotMoveEndInd; ftNotMoveEndInd];
        ftNotMovingXT = fictracProc.t(ftNotMovingX);
        ftY0 = zeros(size(ftNotMoveStartInd));
        ftY1 = ones(size(ftNotMoveStartInd));
        ftNotMovingY = [ftY0; ftY1; ftY1; ftY0];
        
        % plot right leg
        cla(r2Ax);
        plot(r2Ax,legTrack.t, legTrack.srnLegX(:,r2LegInd));

        % plot dots for all not move ind
        plot(r2Ax,legTrack.t(legNotMoveInd), legTrack.srnLegX(legNotMoveInd,r2LegInd), ...
            '.','LineStyle','none');
        % plot shading for not moving bouts
        patch(r2Ax, legNotMovingXT, legNotMovingY, 'black', 'FaceAlpha', 0.3');


        % plot left leg
        cla(l2Ax);
        plot(l2Ax, legTrack.t, legTrack.srnLegX(:,l2LegInd));

        % plot dots for all not move ind
        plot(l2Ax,legTrack.t(legNotMoveInd), legTrack.srnLegX(legNotMoveInd,l2LegInd), ...
            '.','LineStyle','none');
        % plot shading for not moving bouts
        patch(l2Ax,legNotMovingXT, legNotMovingY, 'black', 'FaceAlpha', 0.3');

        % plot FicTrac
        cla(ftAx);
        plot(ftAx, fictracProc.t, smoTotSpdNorm);
        hold on;
        % plot shading for not moving bouts
        patch(ftAx, ftNotMovingXT, ftNotMovingY, 'black','FaceAlpha',0.3);
    end


    % button for when user is done
    doneButton = uicontrol(f, 'Style', 'pushbutton', 'String', 'Done', ...
        'Position', [100 30 50 20]);
    doneButton.Callback = @donePushed;
    
    % function for when user is done
    function donePushed(src, event)
        % toggle logical to done, stops loop
        userDone = 1; 
    end
    
    % loop until user hits done button
    while ~userDone
        pause(0.1);
    end




    % get all return values, leg
    legNotMoveBout = [legNotMoveStartInd, legNotMoveEndInd];
    
    allLegInd = (1:length(legTrack.t))'; % all indicies in trial
    
    % move indicies are inverse of not move indicies
    legMoveInd = allLegInd(~ismember(allLegInd,legNotMoveInd));
    % convert indicies to bout starts and ends
    [legMoveBoutStartInd, legMoveBoutEndInd, ~] = findBouts(legMoveInd);
    legMoveBout = [legMoveBoutStartInd, legMoveBoutEndInd];

    % get all return values, FicTrac
    ftNotMoveBout = [ftNotMoveStartInd, ftNotMoveEndInd];

    allFtInd = (1:length(fictracProc.t))';

    ftMoveInd = allFtInd(~ismember(allFtInd,ftNotMoveInd));
    [ftMoveBoutStartInd, ftMoveBoutEndInd, ~] = findBouts(ftMoveInd);
    ftMoveBout = [ftMoveBoutStartInd, ftMoveBoutEndInd];
    
    close(f);
end