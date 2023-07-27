% saveStepParamPairwiseCorr_bouts.m
%
% Function that takes output from saveBallLegStepParamCond_bouts() or
%  saveLegStepParamCond_bouts() and gets pairwise correlations (Pierson)
%  between all specified leg step parameters (specific param, leg, step).
% Also saves all points that go into correlations.
% Select input file through GUI
%
% INPUTS:
%   corrVars - struct of all vars to correlate, 1 entry in each field per
%       variable
%     params - name of leg step parameters to correlate
%     legs - name of legs (R1, R2, R3, L1, L2, L3)
%     whichStep - which step (0 for at yaw peak, neg for before, pos for
%       after)
%     whichPhase - which phase ('swing' or 'stance')
%   datDir - full path to input
%   bootN - number of times to bootstrap
%   saveFileName - name of output file
%   saveFileDir - full path to directory to save output file
%
% OUTPUTS:
%   none, but saves file with output correlations
%
% CREATED: 7/23/23 - HHY
%
% UPDATED:
%   7/23/23 - HHY
%
function saveStepParamPairwiseCorr_bouts(corrVars, datDir, bootN, saveFileName, ...
    saveFileDir)

    legIDs.name = {'R1', 'R2', 'R3', 'L1', 'L2', 'L3'};
    legIDs.ind = 1:6;

    circStepParams = {'stepDirections'};

    % prompt user to select cond_bout() file
    [condBoutFName, condBoutPath] = uigetfile('*.mat', ...
        'Select cond_bout file', datDir, 'MultiSelect', 'off');

    % load data from cond_bout file
    fullFilePath = [condBoutPath filesep condBoutFName];

    load(fullFilePath, 'selStanceParams', 'selSwingParams', 'maxNumSteps');

    
    % number of variables for pairwise correlation
    numVars = length(corrVars.params);

    % number of turn bouts
    numTurnBouts = size(selStanceParams.stepLengths, 3);

    % preallocate
    allVarNames = cell(numVars, 1);
    corrMatrix = zeros(numVars, numVars);
    corrMatrixRand = zeros(numVars, numVars, bootN);
    indivTurns = zeros(numTurnBouts, numVars); 
    rmvTurnInd = [];

    % get all variable names
    for i = 1:numVars
        % variable name is all 4 conditions
        thisVarName = [corrVars.params{i} '\_' corrVars.legs{i} '\_' ...
            'step' num2str(corrVars.whichStep(i)) '\_'...
            corrVars.whichPhase{i}];

        allVarNames{i} = thisVarName;
    end

   
    % loop through all turn bouts
    for i = 1:numTurnBouts
        % loop through all variables to correlate
        for j = 1:numVars

            % convert leg string to index
            thisLegInd = legIDs.ind(strcmpi(corrVars.legs{j}, legIDs.name));

            % convert step designation to index
            thisStepInd = maxNumSteps + 1 + corrVars.whichStep(j);
    
            % get the value for this variable, this turn
            % swing or stance
            if (strcmpi(corrVars.whichPhase{j}, 'stance'))
                thisVarParam = selStanceParams.(corrVars.params{j});
            elseif (strcmpi(corrVars.whichPhase{j}, 'swing'))
                thisVarParam = selSwingParams.(corrVars.params{j});
            end

            thisVarVal = thisVarParam(thisStepInd, thisLegInd, i);

            % check if this variable value is NaN
            if (isnan(thisVarVal))
                rmvTurnInd = [rmvTurnInd; i];
                continue;
            end

            % update matrix
            indivTurns(i,j) = thisVarVal;
        end
    end


    % remove all turn bouts where at least one variable value is NaN
    indivTurns(rmvTurnInd, :) = [];

    % remove outliers
    indivTurns = rmoutliers(indivTurns);


    % compute all pairwise correlations
    for i = 1:numVars
        % check if this variable is circular
        iCirc = any(strcmpi(corrVars.params{i}, circStepParams));
        for j = 1:numVars
            % check if this variable is circular
            jCirc = any(strcmpi(corrVars.params{j}, circStepParams));

            % compute correlations, handle circular appropriately
            if (iCirc && jCirc) % both circular
                [thisCorrCoef,~] = circ_corrcc(deg2rad(indivTurns(:,i)),...
                    deg2rad(indivTurns(:,j)));
            elseif (iCirc && ~jCirc) % one circular, one linear
                [thisCorrCoef,~] = circ_corrcl(deg2rad(indivTurns(:,i)),...
                    indivTurns(:,j));
            elseif (~iCirc && jCirc) % other circular, one linear
                [thisCorrCoef,~] = circ_corrcl(deg2rad(indivTurns(:,j)),...
                    indivTurns(:,i));
            else % both linear
                thisCorrCoef = corrcoef(indivTurns(:,i), indivTurns(:,j));
                thisCorrCoef = thisCorrCoef(1,2);
            end

            corrMatrix(i,j) = thisCorrCoef;
        end
    end

    % compute pairwise correlations, turnIDs random draws
    for k = 1:bootN
        for i = 1:numVars
            % check if this variable is circular
            iCirc = any(strcmpi(corrVars.params{i}, circStepParams));
            for j = 1:numVars
                % check if this variable is circular
                jCirc = any(strcmpi(corrVars.params{j}, circStepParams));

                % random draws
                var1 = indivTurns(:,i);
                var2 = indivTurns(:,j);

                randIndVar1 = randi(length(var1), length(var1), 1);
                randIndVar2 = randi(length(var2), length(var2), 1);

                randVar1 = var1(randIndVar1);
                randVar2 = var2(randIndVar2);
    
                % compute correlations, handle circular appropriately
                if (iCirc && jCirc) % both circular
                    [thisCorrCoef,~] = circ_corrcc(deg2rad(randVar1),...
                        deg2rad(randVar2));
                elseif (iCirc && ~jCirc) % one circular, one linear
                    [thisCorrCoef,~] = circ_corrcl(deg2rad(randVar1),...
                        randVar2);
                elseif (~iCirc && jCirc) % other circular, one linear
                    [thisCorrCoef,~] = circ_corrcl(deg2rad(randVar2),...
                        randVar1);
                else % both linear
                    thisCorrCoef = corrcoef(randVar1, randVar2);
                    thisCorrCoef = thisCorrCoef(1,2);
                end
    
                corrMatrixRand(i,j,k) = thisCorrCoef;
            end
        end
    end


    % save all output vars
    saveFileFullPath = [saveFileDir filesep saveFileName '.mat'];

    save(saveFileFullPath, 'allVarNames', 'corrMatrix', 'indivTurns', ...
        'corrMatrixRand', 'fullFilePath', '-v7.3');
end