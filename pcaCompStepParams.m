% pcaCompStepParams.m
%
% Function that takes in two sets of step parameters, performs PCA on each,
%  and compares the two against each other.
% Operates on output from saveBallLegStepParamCond_bouts() or
%  saveLegStepParamCond_bouts(). One point per turn bout
%
% INPUTS:
%   set1 - struct of first set of step params
%     params - cell array of step param names, 1 per variable
%     legs - cell array of leg names (R1-3, L1-3), 1 per variable
%     whichStep - which step (0 for at yaw peak, neg for before, pos for
%       after)
%     whichPhase - which phase ('swing' or 'stance')
%     params - cell array of step param names, 1 per variable
%     legs - cell array of leg names (R1-3, L1-3), 1 per variable
%     whichStep - which step (0 for at yaw peak, neg for before, pos for
%       after)
%     whichPhase - which phase ('swing' or 'stance')
%   set2 - struct of second set of step params
%   datDir - full path to input
%   saveFileName - name of output file
%   saveFileDir - full path to directory to save output file
%
% OUTPUTS:
%   none, but saves to output file
%
% CREATED: 7/27/23 - HHY
%
% UPDATED:
%   7/27/23 - HHY
%
function pcaCompStepParams(set1, set2, datDir, saveFileName, saveFileDir)

    legIDs.name = {'R1', 'R2', 'R3', 'L1', 'L2', 'L3'};
    legIDs.ind = 1:6;

    % prompt user to select cond_bout() file
    [condBoutFName, condBoutPath] = uigetfile('*.mat', ...
        'Select cond_bout file', datDir, 'MultiSelect', 'off');

    % load data from cond_bout file
    fullFilePath = [condBoutPath filesep condBoutFName];

    load(fullFilePath, 'selStanceParams', 'selSwingParams', 'maxNumSteps');

    % number of step parameter variables
    set1NumVars = length(set1.params);
    set2NumVars = length(set2.params);
    totNumVars = set1NumVars + set2NumVars;

    % number of turn bouts
    numTurnBouts = size(selStanceParams.stepLengths, 3);

    % preallocate
    indivTurns = zeros(numTurnBouts, totNumVars); 
    rmvTurnInd = [];

   
    % loop through all turn bouts
    for i = 1:numTurnBouts
        % loop through all set1 variables
        for j = 1:set1NumVars

            % convert leg string to index
            thisLegInd = legIDs.ind(strcmpi(set1.legs{j}, legIDs.name));

            % convert step designation to index
            thisStepInd = maxNumSteps + 1 + set1.whichStep(j);
    
            % get the value for this variable, this turn
            % swing or stance
            if (strcmpi(set1.whichPhase{j}, 'stance'))
                thisVarParam = selStanceParams.(set1.params{j});
            elseif (strcmpi(set1.whichPhase{j}, 'swing'))
                thisVarParam = selSwingParams.(set1.params{j});
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

        % loop through all set2 variables
        for j = 1:set2NumVars
            % convert leg string to index
            thisLegInd = legIDs.ind(strcmpi(set2.legs{j}, legIDs.name));

            % convert step designation to index
            thisStepInd = maxNumSteps + 1 + set2.whichStep(j);
    
            % get the value for this variable, this turn
            % swing or stance
            if (strcmpi(set2.whichPhase{j}, 'stance'))
                thisVarParam = selStanceParams.(set2.params{j});
            elseif (strcmpi(set2.whichPhase{j}, 'swing'))
                thisVarParam = selSwingParams.(set2.params{j});
            end

            thisVarVal = thisVarParam(thisStepInd, thisLegInd, i);

            % check if this variable value is NaN
            if (isnan(thisVarVal))
                rmvTurnInd = [rmvTurnInd; i];
                continue;
            end

            % update matrix
            indivTurns(i,j + set1NumVars) = thisVarVal;
        end
    end


    % remove all turn bouts where at least one variable value is NaN
    rmvTurnInd = unique(rmvTurnInd); % no repeats
    indivTurns(rmvTurnInd, :) = [];

    % remove outliers
    indivTurns = rmoutliers(indivTurns);

    
    % separate indivTurns matrix into two matricies, set1 and set2
    set1.indivTurns = indivTurns(:, 1:set1NumVars);
    set2.indivTurns = indivTurns(:, ...
        (set1NumVars + 1):(set1NumVars + set2NumVars));

    % perform PCA
    [set1.coeff, set1.score, set1.latent, set1.tsquared, ...
        set1.explained, set1.mu] = pca(set1.indivTurns);

    [set2.coeff, set2.score, set2.latent, set2.tsquared, ...
        set2.explained, set2.mu] = pca(set2.indivTurns);

    % save output
    saveFileFullPath = [saveFileDir filesep saveFileName '.mat'];

    save(saveFileFullPath, 'set1', 'set2', '-v7.3');

    % print some outputs to screen
    fprintf('Set 1, PC1, variance explained = %.2f%%\n', set1.explained(1));
    fprintf('Set 2, PC1, variance explained = %.2f%%\n', set2.explained(1));
end