% saveStepParamMultiCorr.m
%
% Function that takes output from saveBallLegStepParamCond_bouts() or
%  saveLegStepParamCond_bouts() and gets R2 from multiple linear regression
%  between specified leg step parameters. For all specified variables,
%  calculate regression where each one is independent variable being
%  predicted by the rest as dependent variables
% Also saves all points that go into correlations.
% Select input file through GUI
%
% NOTE: doesn't work on circular parameters
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
% CREATED: 7/25/23 - HHY
%
% UPDATED:
%   7/25/23 - HHY
%
function saveStepParamMultiCorr(corrVars, datDir, bootN, saveFileName, ...
    saveFileDir)

    legIDs.name = {'R1', 'R2', 'R3', 'L1', 'L2', 'L3'};
    legIDs.ind = 1:6;

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
    corrVals = zeros(numVars, 1);
    corrValsRand = zeros(numVars, bootN);
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


    % compute all multiple regressions
    for i = 1:numVars
        indptVar = indivTurns(:,i);
        deptInd = setdiff(1:numVars,i);
        deptVar = indivTurns(:,deptInd);

        % add column of ones to deptVar
        deptVar = [deptVar ones(size(indptVar))];

        % compute multiple regression
        [~,~,~,~,stats] = regress(indptVar, deptVar);

        corrVals(i) = stats(1);
    end

    % compute multiple regressions, turnIDs random draws
    for k = 1:bootN
        for i = 1:numVars
            % randomize independent var
            indptVar = indivTurns(:,i);
            randInd = randi(length(indptVar), length(indptVar), 1);
            randIndptVar = indptVar(randInd);

            % randomize dependent var
            deptInd = setdiff(1:numVars,i);
            deptVar = indivTurns(:,deptInd);

            randIndDept = randi(size(deptVar,1), size(deptVar,1), ...
                size(deptVar,2));

            for j = 1:size(deptVar,2)
                randDeptVar(:,j) = deptVar(randIndDept(:,j),j);
            end

            % add column of ones to deptVar
            randDeptVar = [randDeptVar ones(size(randIndptVar))];

            % compute multiple regression
            [~,~,~,~,stats] = regress(randIndptVar, randDeptVar);

            corrValsRand(i,k) = stats(1);
        end
    end


    % save all output vars
    saveFileFullPath = [saveFileDir filesep saveFileName '.mat'];

    save(saveFileFullPath, 'allVarNames', 'corrVals', 'indivTurns', ...
        'corrValsRand', 'fullFilePath', '-v7.3');
end