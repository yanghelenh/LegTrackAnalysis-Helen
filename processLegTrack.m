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

    % parameters for smoothing for determining leg velocities
    smoParams.padLen = 50; % pad length in samples
    smoParams.sigma = 10; % in samples
    

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
    
    
    % preprocess .trk file - leg positions, velocities
    legTrack = preprocessLegTrack(trkFullPath, leg.frameTimes, ...
        refPts, smoParams);
end