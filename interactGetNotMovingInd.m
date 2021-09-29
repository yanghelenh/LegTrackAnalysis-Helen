% interactGetNotMovingInd.m
%
% Function for selecting not-moving times. 
% Brings up GUI with plot of mid-leg positions, shaded regions of not 
%  moving & red dots on not-moving, and sliders to adjust all parameters 
%  that determine not-moving portions.
%
% INPUTS:
%   legTrack - struct of leg tracking data, output of preprocessLegTrack
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
% CREATED: 9/24/21 - HHY
%
% UPDATED:
%   9/28/21 - HHY
%
function [notMoveInd, notMoveBout, moveInd, moveBout, notMoveParams] = ...
    interactGetNotMovingInd(legTrack, notMoveParams, r2LegInd, l2LegInd)

    % some parameters
    tRange = 30; % sec, amount of data to display
    xMax = legTrack.t(end); % max value
    
    % slider parameters
    sldXPos = 1200;
    sldYPosStart = 850;
    sldYPosEnd = 50;
    sldHeight = 20;
    sldWidth = 300;
    numSld = length(fieldnames(notMoveParams));
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

    % get not moving calls with initial parameters
    [zeroVelInd, notMoveStartInd, notMoveEndInd] = findFlyNotMoving(...
        legTrack.legXVel, legTrack.legYVel, notMoveParams, r2LegInd, ...
        l2LegInd);
    
    % get shading for not moving
    notMovingX = [notMoveStartInd'; notMoveStartInd'; notMoveEndInd'; ...
        notMoveEndInd'];
    notMovingXT = legTrack.t(notMovingX);
    y0 = ones(size(notMoveStartInd)) * -0.6;
    y1 = ones(size(notMoveStartInd)) * 0.6;
    notMovingY = [y0'; y1'; y1'; y0'];
    
    
    % initialize figure
    f = figure('Position', [20 20 1600 920]);

    
    % plot right leg
    r2Ax = subplot('Position', [0.05 0.65 0.6 0.3]);
    plot(legTrack.t, legTrack.srnLegX(:,r2LegInd));
    hold on;
    
    % plot dots for all not move ind
    plot(legTrack.t(zeroVelInd), legTrack.srnLegX(zeroVelInd,r2LegInd), ...
        '.','LineStyle','none');
    % plot shading for not moving bouts
    patch(notMovingXT, notMovingY, 'black', 'FaceAlpha', 0.3');
    
    xlim([0 tRange]);
    title('R2 X position');
    xlabel('Time (s)');
    
    
    % plot left leg
    l2Ax = subplot('Position', [0.05 0.25 0.6 0.3]);
    plot(legTrack.t, legTrack.srnLegX(:,l2LegInd));
    hold on;
    
    % plot dots for all not move ind
    plot(legTrack.t(zeroVelInd), legTrack.srnLegX(zeroVelInd,l2LegInd), ...
        '.','LineStyle','none');
    % plot shading for not moving bouts
    patch(notMovingXT, notMovingY, 'black', 'FaceAlpha', 0.3');
    
    xlim([0 tRange]);
    title('L2 X position');
    xlabel('Time (s)');
    
    
    % link axes for left and right legs
    linkaxes([r2Ax l2Ax], 'xy');
    
    
    % get all notMoveParam names
    nmpNames = fieldnames(notMoveParams);
    
    
    % text for names of parameter sliders
    % invisible axes for text
    txtAx = axes('Position',[0 0 1 1], 'Visible', 'off');
    set(gcf, 'CurrentAxes', txtAx);
    for i = 1:length(nmpNames)
        thisSldYPos = sldYPosStart - sldYSpace * (i - 1);
        
        text(sldXPos - 120, thisSldYPos + 10, nmpNames{i}, ...
             'Units', 'pixels', 'FontSize', 12);
    end
    % text for values of parameter sliders
    allTxtH = {}; % handles to text obj
    for i = 1:length(nmpNames)
        thisSldYPos = sldYPosStart - sldYSpace * (i - 1);
        
        thisDispVal = num2str(notMoveParams.(nmpNames{i}));
        
        allTxtH{i} = text(sldXPos + 320, thisSldYPos + 10, thisDispVal, ...
             'Units', 'pixels', 'FontSize', 12);
    end
    
    % slider adjusting view of leg positions
    tSlider = uicontrol(f, 'Style', 'slider', 'Position', [100 150 600 20]);
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
    end
    
    
    % sliders for adjusting each of the parameter values
    % initialize cell array for slider objects
    allSld = {};
    % loop through all parameters, creating sliders
    for i = 1:length(nmpNames)
        allSld{i} = uicontrol(f, 'Style', 'slider');
        
        thisSldYPos = sldYPosStart - sldYSpace * (i - 1);
        
        allSld{i}.Position = [sldXPos thisSldYPos sldWidth sldHeight];
        
        allSld{i}.Value = notMoveParams.(nmpNames{i});
        thisParamRange = nmpRanges.(nmpNames{i});
        allSld{i}.Max = thisParamRange(2);
        allSld{i}.Min = thisParamRange(1);
        allSld{i}.Callback = {@updateGraph, i};  
    end
    
    % function for updating the figure, notMoveParams, notMove indicies
    %  every time a slider is moved
    function updateGraph(src, event, nameInd)
        
        % for those parameters that must be integers, round slider value
        switch nmpNames{nameInd}
            case {'medFiltNumSamps', 'minBoutLen', 'maxTimeFromStep', ...
                    'adjBoutSepInit', 'adjBoutSepEnd'}
                thisVal = round(allSld{nameInd}.Value);
            otherwise
                thisVal = allSld{nameInd}.Value;
        end

        % change the apppropriate parameter value in notMoveParams
        notMoveParams.(nmpNames{nameInd}) = thisVal;
        
        % find new not moving regions
        [zeroVelInd, notMoveStartInd, notMoveEndInd] = findFlyNotMoving(...
            legTrack.legXVel, legTrack.legYVel, notMoveParams, r2LegInd, ...
            l2LegInd);
        
        % update patches
        notMovingX = [notMoveStartInd'; notMoveStartInd'; notMoveEndInd'; ...
            notMoveEndInd'];
        notMovingXT = legTrack.t(notMovingX);
        y0 = ones(size(notMoveStartInd)) * -0.6;
        y1 = ones(size(notMoveStartInd)) * 0.6;
        notMovingY = [y0'; y1'; y1'; y0'];
        
        % plot right leg
        cla(r2Ax);
        plot(r2Ax,legTrack.t, legTrack.srnLegX(:,r2LegInd));

        % plot dots for all not move ind
        plot(r2Ax,legTrack.t(zeroVelInd), legTrack.srnLegX(zeroVelInd,r2LegInd), ...
            '.','LineStyle','none');
        % plot shading for not moving bouts
        patch(r2Ax, notMovingXT, notMovingY, 'black', 'FaceAlpha', 0.3');


        % plot left leg
        cla(l2Ax);
        plot(l2Ax, legTrack.t, legTrack.srnLegX(:,l2LegInd));

        % plot dots for all not move ind
        plot(l2Ax,legTrack.t(zeroVelInd), legTrack.srnLegX(zeroVelInd,l2LegInd), ...
            '.','LineStyle','none');
        % plot shading for not moving bouts
        patch(l2Ax,notMovingXT, notMovingY, 'black', 'FaceAlpha', 0.3');
        
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
        
        % find new not moving regions
        [zeroVelInd, notMoveStartInd, notMoveEndInd] = findFlyNotMoving(...
            legTrack.legXVel, legTrack.legYVel, notMoveParams, r2LegInd, ...
            l2LegInd);
        
        % update all slider values
        for j = 1:length(nmpNames)
            allSld{j}.Value = notMoveParams.(nmpNames{j});
            allTxtH{j}.String = num2str(notMoveParams.(nmpNames{j}));
        end
        
        % update plot
        % update patches
        notMovingX = [notMoveStartInd'; notMoveStartInd'; notMoveEndInd'; ...
            notMoveEndInd'];
        notMovingXT = legTrack.t(notMovingX);
        y0 = ones(size(notMoveStartInd)) * -0.6;
        y1 = ones(size(notMoveStartInd)) * 0.6;
        notMovingY = [y0'; y1'; y1'; y0'];
        
        % plot right leg
        cla(r2Ax);
        plot(r2Ax,legTrack.t, legTrack.srnLegX(:,r2LegInd));

        % plot dots for all not move ind
        plot(r2Ax,legTrack.t(zeroVelInd), legTrack.srnLegX(zeroVelInd,r2LegInd), ...
            '.','LineStyle','none');
        % plot shading for not moving bouts
        patch(r2Ax, notMovingXT, notMovingY, 'black', 'FaceAlpha', 0.3');


        % plot left leg
        cla(l2Ax);
        plot(l2Ax, legTrack.t, legTrack.srnLegX(:,l2LegInd));

        % plot dots for all not move ind
        plot(l2Ax,legTrack.t(zeroVelInd), legTrack.srnLegX(zeroVelInd,l2LegInd), ...
            '.','LineStyle','none');
        % plot shading for not moving bouts
        patch(l2Ax,notMovingXT, notMovingY, 'black', 'FaceAlpha', 0.3');
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

    % get all parameters
    notMoveInd = zeroVelInd;
    notMoveBout = [notMoveStartInd, notMoveEndInd];
    
    allInd = (1:length(legTrack.t))'; % all indicies in trial
    
    % move indicies are inverse of not move indicies
    moveInd = allInd(~ismember(allInd,zeroVelInd));
    % convert indicies to bout starts and ends
    [moveBoutStartInd, moveBoutEndInd, ~] = findBouts(moveInd);
    moveBout = [moveBoutStartInd, moveBoutEndInd];
    
    close(f);
end