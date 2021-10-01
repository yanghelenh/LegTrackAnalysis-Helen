% processLegTrack.m
%
% Call this function to process leg tracking data output from APT
%
% Takes .trk file output from APT and returns processed leg tracking data.
%  Leg positions, velocities, step parameters
%
% NOTE: work in progress as additional analyses added
%
% INPUTS:
%   none, but GUI prompts throughout
%
% OUTPUTS:
%   none, but saves processed data into pData file for experiment it comes
%       from (see 2P or ephys analysis code)
% 
% CREATED: 9/24/21 - HHY
%
% UPDATED:
%   9/24/21 - HHY
%
function processLegTrack()

    % SOME CONSTANT PARAMETERS
    
    % indicies for specific body parts
    % fit line to head and 3 pts defining midpt b/w legs
    refPts.lineFitInd = 7:10; 
    refPts.midPtInd = 9;
    refPts.headPtInd = 7;
    refPts.abdPtInd = 11;
    
    % indicies for specific legs
    r2LegInd = 2;
    l2LegInd = 5;
    
    % names and indicies for legs
    legIDs.ind = 1:6; % indicies into raw position matricies for legs
    legIDs.names = {'R1', 'R2', 'R3', 'L1', 'L2', 'L3'};
    
    % parameters for smoothing for determining leg velocities
    smoParams.padLen = 50; % pad length in samples
    smoParams.sigma = 10; % in samples
    
    % INITIAL PARAMETER VALUES - these change throughout this code
    
    % parameters - move/not move determination
    notMoveParams.medFiltNumSamps = 10;
    notMoveParams.zeroXVelThreshMed = 0.0004; %0.004;
    notMoveParams.zeroXVelThresh = 0.02;
    notMoveParams.movePosXVelThresh = 0.01;
    notMoveParams.moveNegXVelThresh = -0.02;
    notMoveParams.zeroYVelThreshMed = 0.004; %0.008;
    notMoveParams.zeroYVelThresh = 0.04;
    notMoveParams.movePosYVelThresh = 0.01;
    notMoveParams.moveNegYVelThresh = -0.02;
    notMoveParams.minBoutLen = 10; % in samples
    notMoveParams.stepNegXVelThresh = - 0.035;
    notMoveParams.maxTimeFromStep = 100; % in samples
    % merge any not-moving bouts less than this many samples apart
    notMoveParams.adjBoutSepInit = 50;  
    notMoveParams.adjBoutSepEnd = 50;
    
    % parameters - max/min determination
    % length of window, in frames, for moving average
    legRevParams.movAvgWinLen = 30; 
    % length of window, in frames, for finding max/min
    legRevParams.maxminWinLen = 30; 
    legRevParams.adjThresh = 4; % threshold for adjacent indicies
    % threshold of what 1st derivative (vel) should exceed in the positive
    %  direction for a position maximum
    legRevParams.maxPosVelThresh = 0.0001; 
    legRevParams.maxNegVelThresh = -0.001; % in negative direction
    legRevParams.minPosVelThresh = 0.0001; % in pos dir, for leg pos min
    legRevParams.minNegVelThresh = -0.001; % in neg dir, for leg pos min
    % num of frames before and after max/min to check for velocity thresh
    legRevParams.numNegVelFrames = 8;
    legRevParams.numPosVelFrames = 12;
    
    

    % prompt user for .trk file; defaults to folder containing trk files
    disp('Select .trk file to process');
    [trkName, trkPath] = uigetfile('*.trk', 'Select .trk file', trkPath());
    % full path to .trk file
    trkFullPath = [trkPath filesep trkName];
    
    % find matching pData file
    % prompt user for pData folder, through GUI
    disp('Select pData folder');
    pDataPath = uigetdir(pwd, 'Select pData folder');
    
    % get pData file corresponding to this .trk file
    pDataFilename = fromTrkGetpDataName(trkName, pDataPath);
    
    % load pData
    pDataFilePath = [pDataPath filesep pDataFilename]; % full path
    % load
    load(pDataFilePath, 'fictrac', 'fictracProc', 'fictracParams', 'leg');
    
    % get variables saved in pData file
    pDatVars = whos('-file', pDataFilePath);
    
    % convert pDatVars into cell array of just names
    for i = 1:length(pDatVars)
        pDatVarsNames{i} = pDatVars(i).name;
    end
    
    % check if leg tracking analysis was already run on this file (leg
    %  tracking variables exist in pData)
    if (any(strcmpi(pDatVarsNames, 'legTrack')))
        disp('Leg Tracking analysis run on this pData file before');
        % ask if user wants to continue
        contStr = input('Continue? y/n ', 's');
        % if not, end function here
        if (~strcmpi(contStr, 'y'))
            disp('Ending processLegTrack');
            return;
        end
    end
        
    % preprocess .trk file - leg positions, velocities
    legTrack = preprocessLegTrack(trkFullPath, leg.frameTimes, ...
        refPts, smoParams);
    
    % check if move/not move determined already and if yes, whether user
    %  wants to redo it
    if (any(strcmpi(pDatVarsNames, 'moveNotMove')))
        disp('Moving/not moving bouts already determined');
        % ask if user wants to rerun
        contStr = input('Rerun moving/not moving selection? y/n ', 's');
        % yes, rerun
        if (strcmpi(contStr, 'y'))
    
            % get not moving and moving indicies and bout starts/ends
            [notMoveInd, notMoveBout, moveInd, moveBout, notMoveParams] = ...
                interactGetNotMovingInd(legTrack, notMoveParams, r2LegInd, ...
                l2LegInd);
            % save into moveNotMove struct
            moveNotMove.notMoveInd = notMoveInd;
            moveNotMove.notMoveBout = notMoveBout;
            moveNotMove.moveInd = moveInd;
            moveNotMove.moveBout = moveBout;
            moveNotMove.notMoveParams = notMoveParams;
            
            % update pData
            save(pDataFilePath, 'legTrack', 'moveNotMove', '-append', ...
                '-v7.3');
            
        % if not rerunning, load in previous
        else
            disp('Loading in previous moveNotMove');
            load(pDataFilePath, 'moveNotMove');
        end
    % move/not move not yet determined
    else
        % get not moving and moving indicies and bout starts/ends
        [notMoveInd, notMoveBout, moveInd, moveBout, notMoveParams] = ...
            interactGetNotMovingInd(legTrack, notMoveParams, r2LegInd, ...
            l2LegInd);

        % save into moveNotMove struct
        moveNotMove.notMoveInd = notMoveInd;
        moveNotMove.notMoveBout = notMoveBout;
        moveNotMove.moveInd = moveInd;
        moveNotMove.moveBout = moveBout;
        moveNotMove.notMoveParams = notMoveParams;
            
        % update pData with leg info at this point
        save(pDataFilePath, 'legTrack', 'moveNotMove', '-append', ...
            '-v7.3');
    end
    
    
    % get max and min points of leg position, for step determination
    
    % check if steps already computed
    if (any(strcmpi(pDatVarsNames, 'legSteps')))
        disp('Steps already computed');
        % ask if user wants to redo
        contStr = input('Rerun max/min selection? y/n ', 's');
        
        % yes, redo
        if (strcmpi(contStr, 'y'))
            
            % compute max/min
            [maxIndsAll, minIndsAll, maxWhichLeg, minWhichLeg, ...
                userSelVal] = interactGetLegReversals(legTrack, ...
                moveNotMove, legRevParams, legIDs);
            
            % save into legSteps struct
            legSteps.maxIndsAll = maxIndsAll;
            legSteps.minIndsAll = minIndsAll;
            legSteps.maxWhichLeg = maxWhichLeg;
            legSteps.minWhichLeg = minWhichLeg;
            legSteps.userSelVal = userSelVal;
            
            
            % update pData file
            save(pDataFilePath, 'legSteps', '-append', '-v7.3');
            
        % no, don't redo
        else
            disp('Loading in previous legSteps');
            load(pDataFilePath, 'legSteps');
        end
    % steps not yet computed
    else
        % compute max/min
        [maxIndsAll, minIndsAll, maxWhichLeg, minWhichLeg, ...
            userSelVal] = interactGetLegReversals(legTrack, ...
            moveNotMove, legRevParams, legIDs);

        % save into legSteps struct
        legSteps.maxIndsAll = maxIndsAll;
        legSteps.minIndsAll = minIndsAll;
        legSteps.maxWhichLeg = maxWhichLeg;
        legSteps.minWhichLeg = minWhichLeg;
        legSteps.userSelVal = userSelVal;
        
        % update pData file
        save(pDataFilePath, 'legSteps', '-append', '-v7.3');
    end
    
    
    % compute steps
    [stepInds, stepWhichLeg] = minMaxPosToSteps(legSteps.maxIndsAll, ...
        legSteps.minIndsAll, legSteps.maxWhichLeg, legSteps.minWhichLeg,...
        moveNotMove.notMoveBout, moveNotMove.moveBout);
    
    % save into legSteps struct
    legSteps.stepInds = stepInds;
    legSteps.stepWhichLeg = stepWhichLeg;
    
    
    % compute step parameters
    
    
    
    
    % update pData file
    save(pDataFilePath, 'legSteps', '-append', '-v7.3'); 
end