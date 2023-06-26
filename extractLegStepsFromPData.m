% extractLegStepsFromPData.m
%
% Function that computes leg step parameters from raw leg tracking data
% User selects single pData file through GUI. This pData file must contain
%  the legTrack struct (output of loadTrk2PData()), the fictrac_proc struct
%  (output of filtFictrac_all().
% Calls two interactive figures, one to determine moving/not-moving bouts
%  and the other to correct the leg position max/min calls (needed for
%  determining steps)
% Saves back into same pData file new structs legSteps, stanceStepParams,
%  and swingStepParams
%
% Current function. Effectively replaces processLegTrack()
%
% INPUTS:
%   none, but expects pData file with legTrack and fictrac_proc
%
% OUTPUTS:
%
% CREATED: 7/6/22 - HHY
%
% UPDATED:
%   7/6/22 - HHY
%   7/19/22 - HHY - add times to moveNotMove struct
%   7/27/22 - HHY - bug fixes, option for non-interactive max/min leg
%       position selection
%   2/13/23 - HHY - update compute legSteps to include ephys data
%   6/21/23 - HHY - update to use MATLAB findpeaks() instead of
%       findLegReversals() to find max and min position
%   6/22/23 - HHY - add phase info
%   6/23/23 - HHY - when redoing move/not-move, start from previous values
%
function extractLegStepsFromPData()

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
    % FicTrac notMoveParams
    notMoveParams.ftTotSpdThresh = 0.1;
    notMoveParams.ftMinBoutLen = 0.2;
    notMoveParams.cmbMethod = 'union';
    notMoveParams.sigma = 0.2;
    
    % as of 6/21/23 - outdated parameters
%     % parameters - max/min determination
%     % length of window, in frames, for moving average
%     legRevParams.movAvgWinLen = 30; 
%     % length of window, in frames, for finding max/min
%     legRevParams.maxminWinLen = 30; 
%     legRevParams.adjThresh = 4; % threshold for adjacent indicies
%     % threshold of what 1st derivative (vel) should exceed in the positive
%     %  direction for a position maximum
%     legRevParams.maxPosVelThresh = 0.0001; 
%     legRevParams.maxNegVelThresh = -0.001; % in negative direction
%     legRevParams.minPosVelThresh = 0.0001; % in pos dir, for leg pos min
%     legRevParams.minNegVelThresh = -0.001; % in neg dir, for leg pos min
%     % num of frames before and after max/min to check for velocity thresh
%     legRevParams.numNegVelFrames = 8;
%     legRevParams.numPosVelFrames = 12;

    % default parameters for extracting peaks - max/min position
    %  determination
    legRevParams.minProm = 0.05; % MinPeakProminence of findpeaks
    legRevParams.minDist = 6; % MinPeakDistance of findpeaks


    % prompt user for pData file; defaults to folder containing pData files
    disp('Select pData file to process');
    [pDataName, pDataPath] = uigetfile('*.mat', 'Select pData file', ...
        pDataDir());
    % full path to pData file
    pDataFullPath = [pDataPath filesep pDataName];

    fprintf('Selected %s\n', pDataName);

    % get variables saved in pData file
    pDatVars = whos('-file', pDataFullPath);

    pDatVarsNames = cell(size(pDatVars));
    
    % convert pDatVars into cell array of just names
    for i = 1:length(pDatVars)
        pDatVarsNames{i} = pDatVars(i).name;
    end

    % check that this pData file has all the variables needed for later
    %  analysis
    if (~any(strcmpi(pDatVarsNames, 'legTrack')))
        fprintf('%s does not contain the legTrack struct. Ending \n', ...
            pDataName);
        return;
    elseif (~any(strcmpi(pDatVarsNames, 'fictracProc')))
        fprintf('%s does not contain the fictracProc struct. Ending \n', ...
            pDataName);
        return;
    end

    % whether the pData file has ephysSpikes (only trials w/ephys data)
    if (any(strcmpi(pDatVarsNames, 'ephysSpikes'))) % yes
        % load pData file
        load(pDataFullPath, 'legTrack', 'fictracProc', 'ephysSpikes');
    else % no ephysSpikes, make empty vector
        load(pDataFullPath, 'legTrack', 'fictracProc');
        ephysSpikes = [];
    end

    % Moving/not moving %

    % check if move/not move determined already and if yes, whether user
    %  wants to redo it
    if (any(strcmpi(pDatVarsNames, 'moveNotMove')))
        fprintf('Moving/not moving bouts already determined for %s \n', ...
            pDataName);
        % ask if user wants to redo
        contStr = input('Rerun moving/not moving selection? y/n ', 's');
        % yes, rerun
        if (strcmpi(contStr, 'y'))

            % load previous move/not move params
            load(pDataFullPath, 'moveNotMove');
            notMoveParams = moveNotMove.notMoveParams;
    
            % get not moving and moving indicies and bout starts/ends for
            %  leg and Fictrac
            [legNotMoveInd, legNotMoveBout, legMoveInd, legMoveBout, ...
                ftNotMoveInd, ftNotMoveBout, ftMoveInd, ftMoveBout, ...
                notMoveParams] = interactGetNotMovingIndWFt(legTrack, ...
                fictracProc, notMoveParams, r2LegInd, l2LegInd);

            % save into moveNotMove struct
            moveNotMove.legNotMoveInd = legNotMoveInd;
            moveNotMove.legNotMoveBout = legNotMoveBout;
            moveNotMove.legMoveInd = legMoveInd;
            moveNotMove.legMoveBout = legMoveBout;
            moveNotMove.ftNotMoveInd = ftNotMoveInd;
            moveNotMove.ftNotMoveBout = ftNotMoveBout;
            moveNotMove.ftMoveInd = ftMoveInd;
            moveNotMove.ftMoveBout = ftMoveBout;
            moveNotMove.notMoveParams = notMoveParams;
            
            % update pData
            save(pDataFullPath, 'moveNotMove', '-append');
            
        % if not rerunning, load in previous
        else
            disp('Loading in previous moveNotMove');
            load(pDataFullPath, 'moveNotMove');
        end
    else
        % get not moving and moving indicies and bout starts/ends for
        %  leg and Fictrac
        [legNotMoveInd, legNotMoveBout, legMoveInd, legMoveBout, ...
            ftNotMoveInd, ftNotMoveBout, ftMoveInd, ftMoveBout, ...
            notMoveParams] = interactGetNotMovingIndWFt(legTrack, ...
            fictracProc, notMoveParams, r2LegInd, l2LegInd);

        % save into moveNotMove struct
        moveNotMove.legNotMoveInd = legNotMoveInd;
        moveNotMove.legNotMoveBout = legNotMoveBout;
        moveNotMove.legMoveInd = legMoveInd;
        moveNotMove.legMoveBout = legMoveBout;
        moveNotMove.legT = legTrack.t;
        moveNotMove.ftNotMoveInd = ftNotMoveInd;
        moveNotMove.ftNotMoveBout = ftNotMoveBout;
        moveNotMove.ftMoveInd = ftMoveInd;
        moveNotMove.ftMoveBout = ftMoveBout;
        moveNotMove.ftT = fictracProc.t;
        moveNotMove.notMoveParams = notMoveParams;
        
        % update pData
        save(pDataFullPath, 'moveNotMove', '-append');
    end


    % Leg Position Reversals %

    % check if steps already computed
    if (any(strcmpi(pDatVarsNames, 'legSteps')))
        fprintf('%s already contains legSteps struct \n', pDataName);
        % ask if user wants to redo
        contStr = input('Rerun max/min selection? y/n ', 's');
        
        % yes, redo
        if (strcmpi(contStr, 'y'))  
            % ask if user wants to use interactive leg reversal selection
            intStr = input('Interactive max/min selection? y/n ', 's');
            if (strcmpi(intStr,'y'))
                % compute max/min interactively
                % updated 6/21/23
                [maxIndsAll, minIndsAll, maxWhichLeg, minWhichLeg, ...
                    userSelVal] = interactGetLegReversals(legTrack, ...
                    moveNotMove, legRevParams, legIDs);
            else
                % compute max/min, non-interactively
                % replaced 6/21/23
%                 [maxIndsAll, minIndsAll, maxWhichLeg, minWhichLeg] = ...
%                     nonInteractGetLegReversals(legTrack, moveNotMove, ...
%                     legRevParams, legIDs);

                [maxIndsAll, minIndsAll, maxWhichLeg, minWhichLeg] = ...
                    getLegReversals(legTrack, moveNotMove, ...
                    legRevParams, legIDs);

                userSelVal = [];
            end
            
            % save into legSteps struct
            legSteps.maxIndsAll = maxIndsAll;
            legSteps.minIndsAll = minIndsAll;
            legSteps.maxWhichLeg = maxWhichLeg;
            legSteps.minWhichLeg = minWhichLeg;
            legSteps.legRevParams = legRevParams;
            legSteps.userSelVal = userSelVal;
            legSteps.legIDs = legIDs;
            
            
            % update pData file
            save(pDataFullPath, 'legSteps', '-append');
            
        % no, don't redo
        else
            disp('Loading in previous legSteps');
            load(pDataFullPath, 'legSteps');
        end
    % steps not yet computed
    else
        % ask if user wants to use interactive leg reversal selection
        intStr = input('Interactive max/min selection? y/n ', 's');
        if (strcmpi(intStr,'y'))
            % compute max/min interactively
            % updated 6/21/23 
            [maxIndsAll, minIndsAll, maxWhichLeg, minWhichLeg, ...
                userSelVal] = interactGetLegReversals(legTrack, ...
                moveNotMove, legRevParams, legIDs);
        else
            % compute max/min, non-interactively
            % 6/21/23 - replace with getLegReversals()
%             [maxIndsAll, minIndsAll, maxWhichLeg, minWhichLeg] = ...
%                 nonInteractGetLegReversals(legTrack, moveNotMove, ...
%                 legRevParams, legIDs);

            [maxIndsAll, minIndsAll, maxWhichLeg, minWhichLeg] = ...
                getLegReversals(legTrack, moveNotMove, legRevParams, ...
                legIDs);

            userSelVal = [];
        end
        
        % save into legSteps struct
        legSteps.maxIndsAll = maxIndsAll;
        legSteps.minIndsAll = minIndsAll;
        legSteps.maxWhichLeg = maxWhichLeg;
        legSteps.minWhichLeg = minWhichLeg;
        legSteps.legRevParams = legRevParams;
        legSteps.userSelVal = userSelVal;
        legSteps.legIDs = legIDs;
        
        
        % update pData file
        save(pDataFullPath, 'legSteps', '-append');
    end

    % compute steps
    [stepInds, stepWhichLeg] = minMaxPosToSteps(legSteps.maxIndsAll, ...
        legSteps.minIndsAll, legSteps.maxWhichLeg, legSteps.minWhichLeg,...
        moveNotMove.legNotMoveBout, moveNotMove.legMoveBout);
    
    % save into legSteps struct
    legSteps.stepInds = stepInds;
    legSteps.stepWhichLeg = stepWhichLeg;
    
    
    % compute step parameters
    legSteps = computeStepParameters(legSteps, legTrack, fictracProc, ...
        ephysSpikes);
    
    % get swing/stance, currently, use duration method
    legSteps = callSwingStanceSteps(legSteps, legTrack, ...
        'duration', fictracProc);
    
    % step parameters during swing/stance
    [stanceStepParams, swingStepParams] = getStepParamsSwingStance(...
        legSteps);

    % get phase, using method of fitting each half step
    legPhase = getLegPhaseFromSteps(legSteps.stepInds, ...
        legSteps.stepWhichLeg, legTrack.srnfLegX(:,legIDs.ind), ...
        moveNotMove.legNotMoveBout);

    % get phase differences b/w legs
    phaseDiffs = getLegPhaseDiffs(legPhase, legIDs, 'degrees');

    % update pData file
    save(pDataFullPath, 'legSteps', 'stanceStepParams', ...
        'swingStepParams', 'legPhase', 'phaseDiffs', '-append');

end