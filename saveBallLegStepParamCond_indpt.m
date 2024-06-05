% saveBallLegStepParamCond_indpt.m
%
% Function that saves legStep parameters for all steps that meet the user 
%  specified conditions on FicTrac
% Treats all steps as independent, as long as FictracSmo condition(s) are
%  met, unlike saveBallLegStepParamCond_bouts.m
% Steps are included if at least for half their duration, the FictracSmo 
%  condition(s) are met. Step is both swing and stance
% If L/R asymmetric conditions are included, pools left and right (i.e.
%  angVel > 50 pools left and right turns). If there are no L/R asymmetric
%  conditions and flipLegsLR is true, swaps left/right leg assignments and
%  inverts L/R asymmetric step parameters. Equivalent to taking the
%  mirror-image for that fly.
% User selects one or more pData files through GUI
%
% Adaptation of saveLegStepParamCond_indpt.m from free walking analysis
% 
% INPUTS:
%   cond - struct of conditions, if multiple conditions, treats it as AND
%       whichParam - cell array (even if 1 element) on which FicTrac field
%           to condition on, one for each condition
%       cond - cell array of strings to condition on, for eval(); same size
%           as whichParam
%   postStimExclDur - additional time, in sec, after stimulation (opto or 
%       iInj to exclude turning bouts from overlapping with
%   flipLegsLR - boolean for whether to flip legs left right. Only works
%       when cond doesn't include L/R asymmetric parameters
%   pDataFNames - cell array of pData file names or [] if select through
%       GUI
%   pDataPath - full path to pData directory
%   saveFilePath - directory in which to save output file
%   saveFileName - name of output file, without .mat part
%
% OUTPUTS:
%   none, but saves output file with name saveFileName in saveFilePath
%       selLegSteps - all steps, all pData files, that meet conditions
%       selStanceParams - all steps, all pData files, that meet condition,
%           stance only
%       selSwingParams - all steps, all pData files, that meet condition,
%           swing only
%       pDataFiles - struct of info on pData files
%           names - name of each pData file with at least 1 valid step, as
%               cell array
%           inds - indices (corresponding to indices in legSteps) that
%               belong to each pData file, as start and end indices
%       legIDs - legSteps.legIDs
%       cond - same as INPUT
%
% CREATED: 7/5/23 - HHY
%
% UPDATED:
%   7/5/23 - HHY
%   8/1/23 - HHY - didn't actually remove stimulation periods, fix
%   4/25/24 - HHY - add flipLegsLR option
%   5/2/24 - HHY - add pDataFNames input, also, save file name is name of
%       first pData fly name + saveFileName as suffix
%
function saveBallLegStepParamCond_indpt(cond, postStimExclDur, ...
    flipLegsLR, pDataFNames, pDataPath, saveFilePath, saveFileName)

    % names of all step parameters to save
    stepParamNames = {'stepLengths', 'stepXLengths', 'stepWhichLeg', ...
        'stepYLengths', 'stepDirections', 'stepDurations', 'stepSpeeds',...
        'stepVelX', 'stepVelY', 'stepAEPX', 'stepAEPY', 'stepPEPX', ...
        'stepPEPY', 'stepFtFwd', 'stepFtLat', 'stepFtYaw'};
    % all the step parameters where values need to be * -1 for left turns
    flipStepParams = {'stepVelY', 'stepAEPY', 'stepPEPY', 'stepFtLat',...
        'stepFtYaw', 'stepDirections'};
    % fictracSmo parameters that are L/R asymmetric contain one of the
    %  following
    lrAsymBTParams = {'yawAngVel', 'slideVel'};

    % prompt user to select pData files
    if isempty(pDataFNames)
        [pDataFNames, pDataDirPath] = uigetfile('*.mat', ...
            'Select pData files', pDataPath, 'MultiSelect', 'on');
    else
        pDataDirPath = pDataPath;
    end
    
    % if only 1 pData file selected, not cell array; make sure loop still
    %  works 
    if (iscell(pDataFNames))
        numPDataFiles = length(pDataFNames);
    else
        numPDataFiles = 1;
    end

    % preallocate 
    pDataFiles.names = pDataFNames;
    pDataFiles.inds = [];

    for i = 1:length(stepParamNames)
        selLegSteps.(stepParamNames{i}) = [];
        selStanceParams.(stepParamNames{i}) = [];
        selSwingParams.(stepParamNames{i}) = [];
    end

    rmvInd = []; % indices of pData files to remove


    % loop through all pData files
    for i = 1:numPDataFiles
    
        % handle whether it's a cell array or not
        if (iscell(pDataFNames))
            pDataName = pDataFNames{i};
        else
            pDataName = pDataFNames;
        end

        pDataFullPath = [pDataDirPath filesep pDataName];

        % save fly name as first pDataName's date, fly, cell (19 characters)
        if (i == 1)
            flyName = pDataName(1:19);
        end


        % get variables saved in pData file
        pDatVars = whos('-file', pDataFullPath);
    
        pDatVarsNames = cell(size(pDatVars));
        
        % convert pDatVars into cell array of just names
        for j = 1:length(pDatVars)
            pDatVarsNames{j} = pDatVars(j).name;
        end

        % check if this pData file has legSteps, fictracSmo, moveNotMove
        %  structs, if not, skip
        if (~any(strcmpi(pDatVarsNames, 'legSteps')) || ...
                ~any(strcmpi(pDatVarsNames, 'fictracSmo')) || ...
                ~any(strcmpi(pDatVarsNames, 'moveNotMove')))
            continue;
        end

        % load variables from pData
        % check if opto or iInj trial, load those as appropriate
        % also, save opto or iInj times as stim times; no stim, empty
        if (any(strcmpi(pDatVarsNames, 'opto'))) % opto stim
            load(pDataFullPath, 'legTrack', 'moveNotMove', 'fictracSmo', ...
                'legSteps', 'stanceStepParams', 'swingStepParams', 'opto');

            stimStartTimes = opto.stimStartTimes;
            stimEndTimes = opto.stimEndTimes;
        elseif (any(strcmpi(pDatVarsNames, 'iInj'))) % current inj
            load(pDataFullPath, 'legTrack', 'moveNotMove', 'fictracSmo', ...
                'legSteps', 'stanceStepParams', 'swingStepParams', 'iInj');

            stimStartTimes = iInj.startTimes;
            stimEndTimes = iInj.endTimes;
        else % no stimulation
            load(pDataFullPath, 'legTrack', 'moveNotMove', 'fictracSmo', ...
                'legSteps', 'stanceStepParams', 'swingStepParams');

            stimStartTimes = [];
            stimEndTimes = [];
        end

        % convert legTrack.refPts to legIDs
        legIDs.ind = legSteps.legIDs.ind;
        legIDs.name = legSteps.legIDs.names;
        

        % get matching b/w corresponding left and right legs
        rightLegInd = find(contains(legIDs.name, 'R'));
        leftLegInd = find(contains(legIDs.name, 'L'));
        matchedLegInd = zeros(length(rightLegInd),2);

        for j = 1:length(rightLegInd)
            thisLegNum = legIDs.name{rightLegInd(j)}(end);
            thisLeftInd = find(contains(legIDs.name(leftLegInd),thisLegNum));
            matchedLegInd(j,1) = rightLegInd(j);
            matchedLegInd(j,2) = leftLegInd(thisLeftInd);
        end

        % check if we're conditioning on any L/R asymmetric parameters
        % if yes, allocate left and right turns separately
        if(any(contains(cond.whichParam, lrAsymBTParams))) 
            sepLR = true; % logical for whether we're separating L/R
            % preallocate logical for valid indices of fictracSmo
            btRightValInd = true(size(fictracSmo.t));
            btLeftValInd = true(size(fictracSmo.t));

            % initialize - keep track of indices of steps to include
            inclStepIndRight = [];
            inclStepIndLeft = [];
        else % if no, no need to separate left and right
            sepLR = false;
            % preallocate logical for valid indices of fictracSmo
            btValInd = true(size(fictracSmo.t));

            % initialize - keep track of indices of steps to include
            inclStepInd = [];
        end

        % loop through all conditions
        for j = 1:length(cond.whichParam)
            % check if conditioning on speed, if yes, get absolute value of
            %  corresponding velocity param
            if (contains(cond.whichParam{j}, 'Spd'))
                condWhichParam = replace(cond.whichParam{j}, 'Spd', 'Vel');
                thisCondVar = abs(fictracSmo.(condWhichParam));
            % otherwise, just copy over paaram    
            else
                thisCondVar = fictracSmo.(cond.whichParam{j});
            end

            % if we're separating left and right (at least 1 param is
            %  asymmetric)
            if(sepLR)
                % get valid indices for right turns
                thisValIndRight = eval(['thisCondVar' cond.cond{j}]);
                % combine with all other conditions
                btRightValInd = btRightValInd & thisValIndRight;

                % if this is a parameter to be flipped, flip it for left
                if(contains(cond.whichParam{j},lrAsymBTParams))
                    thisCondVar = thisCondVar * -1;
                end

                % left turns
                thisValIndLeft = eval(['thisCondVar' cond.cond{j}]);
                btLeftValInd = btLeftValInd & thisValIndLeft;
            else % not separating left and right
                % bodytraj indices where this condition is true
                thisValInd = eval(['thisCondVar' cond.cond{j}]);
    
                % combine with all other conditions
                btValInd = btValInd & thisValInd;
            end
        end

        % filter for moving bouts only
        if(sepLR)
            btRightValInd(moveNotMove.ftNotMoveInd) = false;
            btLeftValInd(moveNotMove.ftNotMoveInd) = false;
        else
            btValInd(moveNotMove.ftNotMoveInd) = false;
        end

        % filter for not during stimulation 
        for j = 1:length(stimStartTimes)
            thisStartTime = stimStartTimes(j);
            thisEndTime = stimEndTimes(j) + postStimExclDur;

            exclInd = find((fictracSmo.t >= thisStartTime) & ...
                (fictracSmo.t <= thisEndTime));

            if(sepLR)
                btRightValInd(exclInd) = false;
                btLeftValInd(exclInd) = false;
            else
                btValInd(exclInd) = false;
            end
        end


        % fictrac IFI
        ftIFI = median(diff(fictracSmo.t));

        % loop through all steps of legSteps, check if they're to be
        %  included
        for j = 1:length(legSteps.stepInds)
            % get duration of this step
            thisStepDur = legSteps.stepT(j,3) - legSteps.stepT(j,1);

            % this step start and end as nearest FicTrac indices
            thisStartInd = find(fictracSmo.t >= legSteps.stepT(j,1), 1, ...
                'first');
            thisEndInd = find(fictracSmo.t <= legSteps.stepT(j,3), 1, ...
                'last');

            if(sepLR)
                % right
                numInclIndRight = sum(btRightValInd(...
                    thisStartInd:thisEndInd)) * ftIFI;
    
                % save this step index if more than half of the time during
                %  cond
                if (numInclIndRight >= thisStepDur/2)
                    inclStepIndRight = [inclStepIndRight; j];
                end

                % left
                numInclIndLeft = sum(btLeftValInd(...
                    thisStartInd:thisEndInd)) * ftIFI;
    
                % save this step index if more than half of the time during
                %  cond
                if (numInclIndLeft >= thisStepDur/2)
                    inclStepIndLeft = [inclStepIndLeft; j];
                end
            else
                numInclInd = sum(btValInd(...
                    thisStartInd:thisEndInd)) * ftIFI;
    
                % save this step index if more than half of the time during
                %  cond
                if (numInclInd >= thisStepDur/2)
                    inclStepInd = [inclStepInd; j];
                end
            end
        end

        if(sepLR)
            % if there are valid steps, save info for this trial
            if (~isempty(inclStepIndRight) || ~isempty(inclStepIndLeft)) 
                % index into selLegSteps, selStanceParams, selSwingParams for start
                %  of steps corresponding to this trial
                trialStartInd = length(selLegSteps.stepWhichLeg) + 1;
        
                % save legStep parameters for this trial
                for j = 1:length(stepParamNames)
                    % for stepWhichLeg, need to swap left and right legs
                    %  for left indices
                    if(strcmpi(stepParamNames{j},'stepWhichLeg'))
                        % right
                        thisLegStep = legSteps.(stepParamNames{j});
                        thisLegStepVals = thisLegStep(inclStepIndRight,:);
                        selLegSteps.(stepParamNames{j}) = ...
                            [selLegSteps.(stepParamNames{j}); thisLegStepVals];
            
                        thisStanceStep = stanceStepParams.(stepParamNames{j});
                        thisStanceVals = thisStanceStep(inclStepIndRight,:);
                        selStanceParams.(stepParamNames{j}) = ...
                            [selStanceParams.(stepParamNames{j}); thisStanceVals];
            
                        thisSwingStep = swingStepParams.(stepParamNames{j});
                        thisSwingVals = thisSwingStep(inclStepIndRight,:);
                        selSwingParams.(stepParamNames{j}) = ...
                            [selSwingParams.(stepParamNames{j}); thisSwingVals];

                        % left
                        thisLegStep = legSteps.(stepParamNames{j});
                        thisLegStepVals = thisLegStep(inclStepIndLeft,:);
                        % swap left and right legs
                        for k = 1:length(thisLegStepVals)
                            thisLegStepVals(k) = swapLeftRightLegs(...
                                thisLegStepVals(k), matchedLegInd);
                        end
                        % add to running paramVals
                        selLegSteps.(stepParamNames{j}) = ...
                            [selLegSteps.(stepParamNames{j}); thisLegStepVals];
            
                        thisStanceStep = stanceStepParams.(stepParamNames{j});
                        thisStanceVals = thisStanceStep(inclStepIndLeft,:);
                        % swap left and right legs
                        for k = 1:length(thisStanceVals)
                            thisStanceVals(k) = swapLeftRightLegs(...
                                thisStanceVals(k), matchedLegInd);
                        end
                        selStanceParams.(stepParamNames{j}) = ...
                            [selStanceParams.(stepParamNames{j}); thisStanceVals];
            
                        thisSwingStep = swingStepParams.(stepParamNames{j});
                        thisSwingVals = thisSwingStep(inclStepIndLeft,:);
                        % swap left and right legs
                        for k = 1:length(thisSwingVals)
                            thisSwingVals(k) = swapLeftRightLegs(...
                                thisSwingVals(k), matchedLegInd);
                        end
                        selSwingParams.(stepParamNames{j}) = ...
                            [selSwingParams.(stepParamNames{j}); thisSwingVals];

                    % if this is a parameter where the value needs to
                    %  be inverted
                    elseif(contains(stepParamNames{j}, flipStepParams))
                        % right - no inversion
                        thisLegStep = legSteps.(stepParamNames{j});
                        thisLegStepVals = thisLegStep(inclStepIndRight,:);
                        selLegSteps.(stepParamNames{j}) = ...
                            [selLegSteps.(stepParamNames{j}); thisLegStepVals];
            
                        thisStanceStep = stanceStepParams.(stepParamNames{j});
                        thisStanceVals = thisStanceStep(inclStepIndRight,:);
                        selStanceParams.(stepParamNames{j}) = ...
                            [selStanceParams.(stepParamNames{j}); thisStanceVals];
            
                        thisSwingStep = swingStepParams.(stepParamNames{j});
                        thisSwingVals = thisSwingStep(inclStepIndRight,:);
                        selSwingParams.(stepParamNames{j}) = ...
                            [selSwingParams.(stepParamNames{j}); thisSwingVals];

                        % left - invert
                        thisLegStep = legSteps.(stepParamNames{j});
                        thisLegStepVals = thisLegStep(inclStepIndLeft,:) * -1;
                        selLegSteps.(stepParamNames{j}) = ...
                            [selLegSteps.(stepParamNames{j}); thisLegStepVals];
            
                        thisStanceStep = stanceStepParams.(stepParamNames{j});
                        thisStanceVals = thisStanceStep(inclStepIndLeft,:) * -1;
                        selStanceParams.(stepParamNames{j}) = ...
                            [selStanceParams.(stepParamNames{j}); thisStanceVals];
            
                        thisSwingStep = swingStepParams.(stepParamNames{j});
                        thisSwingVals = thisSwingStep(inclStepIndLeft,:) * -1;
                        selSwingParams.(stepParamNames{j}) = ...
                            [selSwingParams.(stepParamNames{j}); thisSwingVals];

                    % no need to invert this parameter, not stepWhichLeg
                    else 
                        % combine left and right indices
                        thisInclStepInd = [inclStepIndRight; inclStepIndLeft];

                        thisLegStep = legSteps.(stepParamNames{j});
                        thisLegStepVals = thisLegStep(thisInclStepInd,:);
                        selLegSteps.(stepParamNames{j}) = ...
                            [selLegSteps.(stepParamNames{j}); thisLegStepVals];
            
                        thisStanceStep = stanceStepParams.(stepParamNames{j});
                        thisStanceVals = thisStanceStep(thisInclStepInd,:);
                        selStanceParams.(stepParamNames{j}) = ...
                            [selStanceParams.(stepParamNames{j}); thisStanceVals];
            
                        thisSwingStep = swingStepParams.(stepParamNames{j});
                        thisSwingVals = thisSwingStep(thisInclStepInd,:);
                        selSwingParams.(stepParamNames{j}) = ...
                            [selSwingParams.(stepParamNames{j}); thisSwingVals];
                    end
                end
        
                % end index for steps for this trial
                trialEndInd = length(selLegSteps.stepWhichLeg);
        
                thisTrialInds = [trialStartInd trialEndInd];
        
                pDataFiles.inds = [pDataFiles.inds; thisTrialInds];
            % otherwise, remove this pData file from list included    
            else
                rmvInd = [rmvInd; i];
            end
        else

            % if there are valid steps, save info for this trial
            if (~isempty(inclStepInd)) 
                % index into selLegSteps, selStanceParams, selSwingParams for start
                %  of steps corresponding to this trial
                trialStartInd = length(selLegSteps.stepWhichLeg) + 1;
        
                % save legStep parameters for this trial
                for j = 1:length(stepParamNames)
                    % if flipping left and right legs
                    if (flipLegsLR)
                        % for stepWhichLeg, swap indices of left and right
                        %  legs
                        if (strcmpi(stepParamNames{j}, 'stepWhichLeg'))
                            thisLegStep = legSteps.(stepParamNames{j});
                            thisLegStepVals = thisLegStep(inclStepInd,:);
                            % swap left and right legs
                            for k = 1:length(thisLegStepVals)
                                thisLegStepVals(k) = swapLeftRightLegs(...
                                    thisLegStepVals(k), matchedLegInd);
                            end
                            % add to running paramVals
                            selLegSteps.(stepParamNames{j}) = ...
                                [selLegSteps.(stepParamNames{j}); thisLegStepVals];
                
                            thisStanceStep = stanceStepParams.(stepParamNames{j});
                            thisStanceVals = thisStanceStep(inclStepInd,:);
                            % swap left and right legs
                            for k = 1:length(thisStanceVals)
                                thisStanceVals(k) = swapLeftRightLegs(...
                                    thisStanceVals(k), matchedLegInd);
                            end
                            selStanceParams.(stepParamNames{j}) = ...
                                [selStanceParams.(stepParamNames{j}); thisStanceVals];
                
                            thisSwingStep = swingStepParams.(stepParamNames{j});
                            thisSwingVals = thisSwingStep(inclStepInd,:);
                            % swap left and right legs
                            for k = 1:length(thisSwingVals)
                                thisSwingVals(k) = swapLeftRightLegs(...
                                    thisSwingVals(k), matchedLegInd);
                            end
                            selSwingParams.(stepParamNames{j}) = ...
                                [selSwingParams.(stepParamNames{j}); thisSwingVals];
                        % if this is a parameter where the value needs to
                        %  be inverted
                        elseif(contains(stepParamNames{j}, flipStepParams))
                            thisLegStep = legSteps.(stepParamNames{j});
                            thisLegStepVals = thisLegStep(inclStepInd,:) * -1;
                            selLegSteps.(stepParamNames{j}) = ...
                                [selLegSteps.(stepParamNames{j}); thisLegStepVals];
                
                            thisStanceStep = stanceStepParams.(stepParamNames{j});
                            thisStanceVals = thisStanceStep(inclStepInd,:) * -1;
                            selStanceParams.(stepParamNames{j}) = ...
                                [selStanceParams.(stepParamNames{j}); thisStanceVals];
                
                            thisSwingStep = swingStepParams.(stepParamNames{j});
                            thisSwingVals = thisSwingStep(inclStepInd,:) * -1;
                            selSwingParams.(stepParamNames{j}) = ...
                                [selSwingParams.(stepParamNames{j}); thisSwingVals];
                        % no inversion of this parameter    
                        else
                            thisLegStep = legSteps.(stepParamNames{j});
                            thisLegStepVals = thisLegStep(inclStepInd,:);
                            selLegSteps.(stepParamNames{j}) = ...
                                [selLegSteps.(stepParamNames{j}); thisLegStepVals];
                
                            thisStanceStep = stanceStepParams.(stepParamNames{j});
                            thisStanceVals = thisStanceStep(inclStepInd,:);
                            selStanceParams.(stepParamNames{j}) = ...
                                [selStanceParams.(stepParamNames{j}); thisStanceVals];
                
                            thisSwingStep = swingStepParams.(stepParamNames{j});
                            thisSwingVals = thisSwingStep(inclStepInd,:);
                            selSwingParams.(stepParamNames{j}) = ...
                                [selSwingParams.(stepParamNames{j}); thisSwingVals];
                        end

                    % not flipping legs left-right  
                    % no need to flip leg indices in stepWhichLeg or to
                    %  invert values
                    else
                        thisLegStep = legSteps.(stepParamNames{j});
                        thisLegStepVals = thisLegStep(inclStepInd,:);
                        selLegSteps.(stepParamNames{j}) = ...
                            [selLegSteps.(stepParamNames{j}); thisLegStepVals];
            
                        thisStanceStep = stanceStepParams.(stepParamNames{j});
                        thisStanceVals = thisStanceStep(inclStepInd,:);
                        selStanceParams.(stepParamNames{j}) = ...
                            [selStanceParams.(stepParamNames{j}); thisStanceVals];
            
                        thisSwingStep = swingStepParams.(stepParamNames{j});
                        thisSwingVals = thisSwingStep(inclStepInd,:);
                        selSwingParams.(stepParamNames{j}) = ...
                            [selSwingParams.(stepParamNames{j}); thisSwingVals];
                    end
                end
        
                % end index for steps for this trial
                trialEndInd = length(selLegSteps.stepWhichLeg);
        
                thisTrialInds = [trialStartInd trialEndInd];
        
                pDataFiles.inds = [pDataFiles.inds; thisTrialInds];
            % otherwise, remove this pData file from list included    
            else
                rmvInd = [rmvInd; i];
            end
        end
    end

    % remove unused pData files
    pDataFiles.names(rmvInd) = [];

    % copy over legIDs
    legIDs = legSteps.legIDs;

    % save output file
    fullSavePath = [saveFilePath filesep flyName '_' saveFileName '.mat'];

    save(fullSavePath, 'selLegSteps', 'selStanceParams',...
        'selSwingParams', 'pDataFiles', 'cond', 'legIDs', ...
        'postStimExclDur', 'flipLegsLR', '-v7.3');
end