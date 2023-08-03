% computeEnvelopeStepParamsPData.m
%
% Function that computes a continuous estimate of AEP, PEP, step length,
%  and step direction on multiple selected pData files (selected through 
%  GUI). 
% pData file must have legTrack, legSteps, and moveNotMove structs
% This functionality has now been added to extractLegStepsFromPData() and
%  extractStepsFromPData() (as of 8/3/23), so this function is just to
%  compute these values on older data. Does not need to be run on newer
%  data
%
% INPUTS:
%   pDataPath - path to folder with pData files
%
% OUTPUTS:
%   none, but saves legStepsCont struct into pData file
%
% CREATED: 8/3/23 - HHY
%
% UPDATED: 
%   8/3/23 - HHY
%
function computeEnvelopeStepParamsPData(pDataPath)

    % prompt user to select pData files
    [pDataFNames, pDataDirPath] = uigetfile('*.mat', ...
        'Select pData files', pDataPath, 'MultiSelect', 'on');
    
    % if only 1 pData file selected, not cell array; make sure loop still
    %  works 
    if (iscell(pDataFNames))
        numPDataFiles = length(pDataFNames);
    else
        numPDataFiles = 1;
    end

    % loop through all pData files
    for i = 1:numPDataFiles
    
        % handle whether it's a cell array or not
        if (iscell(pDataFNames))
            pDataName = pDataFNames{i};
        else
            pDataName = pDataFNames;
        end

        pDataFullPath = [pDataDirPath filesep pDataName];

        % get variables saved in pData file
        pDatVars = whos('-file', pDataFullPath);
    
        pDatVarsNames = cell(size(pDatVars));
        
        % convert pDatVars into cell array of just names
        for j = 1:length(pDatVars)
            pDatVarsNames{j} = pDatVars(j).name;
        end

        % check if this pData file has legTrack, legSteps, moveNotMove
        %  structs, if not, skip
        if (~any(strcmpi(pDatVarsNames, 'legTrack')) || ...
                ~any(strcmpi(pDatVarsNames, 'legSteps')) || ...
                ~any(strcmpi(pDatVarsNames, 'moveNotMove')))
            continue;
        end

        % load data
        load(pDataFullPath, 'legTrack', 'legSteps', 'moveNotMove');

        % display which file is being analyzed
        fprintf('Processing %s\n', pDataName);

        % get interpolated step parameters
        legStepsCont = getEnvelopeStepParams(legTrack, legSteps, ...
            moveNotMove);

        % update pData file
        save(pDataFullPath, 'legStepsCont', '-append');
    end
end
