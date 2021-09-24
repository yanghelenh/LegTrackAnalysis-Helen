% fromTrkGetpDataName.m
%
% Function that takes in .trk file name and returns the pData file name for
%  the same experiment
% Assumes standard naming convention for experiments, either 2P imaging or
%  ephys: date_fly_fov/cell_trial
%
% INPUTS:
%   trkFilename - name of .trk file
%   pDataFilepath - path to pData folder
%
% OUTPUTS:
%   pDataFilename - name of pData file
%
% CREATED: 9/24/21 - HHY
%
% UPDATED:
%   9/24/21 - HHY
%
function pDataFilename = fromTrkGetpDataName(trkFilename, pDataFilepath)

    % file name parts separated by _, get indicies of _
    unInd = strfind(trkFilename, '_');
    
    % part of file name we care about is stuff before 4th _ 
    %  date_fly_fov/cell_trial starts file name for both .trk file and
    %  pData file
    trialName = extracBefore(trkFilename, unInd(4));
    
    % find pData file with this trial name
    % get contents of pData folder
    pDataFiles = dir(pDataFilepath);
    
    % cell array for dir output
    pDataDirCell = struct2cell(pDataFiles);
    
    % cell array for just names of pData files - names are first row,
    %  exclude 2 corresponding to . and ..
    pDataNamesCell = pDataDirCell(1,3:end);
    
    % find index corresponding to matching pData file
    pDataInd = find(contains(pDataNamesCell, trialName));
    
    % there should be 1 and only 1 matching pData file
    % too many pData files, throw error
    if (length(pDataInd) > 1)
        ME = MException('multiplePDataFiles', ...
            'More than one pData file names match .trk file name');
        throw(ME);
    % no matching pData files, throw error    
    elseif (isempty(pDataInd))
        ME = MException('noPDataFiles', ...
            'No pData files match this .trk file');
        throw(ME);
    % get pData file name from index
    else
        pDataFilename = pDataNamesCell{pDataInd};
    end
end