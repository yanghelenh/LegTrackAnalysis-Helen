% addMyPaths.m
% 
% function to add relevant ephys data analysis paths, including subfolders
%
% Users: change the following paths to match those on your local computer

function addMyPaths() 
    % Determine which computer this code is running on
    comptype = computer; % get the string describing the computer type
    PC_STRING = 'PCWIN64';  % string for PC on 2P rig
    MAC_STRING = 'MACI64'; %string for macbook
    
    % paths for windows
    if strcmp(comptype, PC_STRING)
        APTpath = 'C:\Users\WilsonLab\Documents\MATLAB\APT'; % APT
        pythonPath = 'C:\ProgramData\Anaconda3\python.exe';
        a2libPath = 'C:\Users\WilsonLab\Documents\HelenExperimentalCode\LegTrackAnalysis-Helen\a2lib';
        analysisPath = 'C:\Users\WilsonLab\Documents\HelenExperimentalCode\LegTrackAnalysis-Helen';
        
    % paths for mac
    elseif strcmp(comptype, MAC_STRING)
        APTpath = '/Users/hyang/Documents/MATLAB/APT'; % APT
        pythonPath = '/anaconda3/bin/python';
        a2libPath = '/Users/hyang/Documents/LegTrackAnalysis-Helen/a2lib';
        analysisPath = '/Users/hyang/Documents/LegTrackAnalysis-Helen';
        
    else
        disp('ERROR: Paths not added. Unrecognized computer type');
    end
    
    %% Add the paths
    addpath(genpath(APTpath));
    addpath(genpath(analysisPath));
    
    %% Python functions, add to python system path
    
    [~, ~, pyLoaded] = pyversion;
    
    if ~(pyLoaded)
        pyversion(pythonPath);
    end
    
    P = py.sys.path;
    if count(P,a2libPath) == 0
        insert(P,int32(0),a2libPath);
    end
    
    pyversion

end