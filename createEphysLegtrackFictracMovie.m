% createEphysLegtrackFictracMovie.m
%
% Quick and dirty script to generate movie showing ephys recording, FicTrac
%  variables, and leg movie. Starts from raw data.
% Made 9/20/21, used for video for U19 year 4 review talk



% full path to data folder


pDataFilePath = ['/Users/hyang/Dropbox (HMS)/EphysAnalysis-Helen/pData/'...
    '200821_fly01_cell01_trial01.mat'];

% full path to leg video
% vidPath = '/Users/hyang/Desktop';
% vidName = 'fly01_cell01_trial01_legVid';

vidPath = '/Users/hyang/Dropbox (HMS)/APTlegTracking/200821_fly01_cell01_trial01_legVid_crop.mp4';

%
% full path to .trk leg tracking data
trkFilePath = ['/Users/hyang/Dropbox (HMS)/APTlegTracking/'...
    'Tracking_v1_mdn_labeled2798_ephys/' ...
    '200821_fly01_cell01_trial01_legVid_crop_caLegTrack_v1_mdn.trk'];

tVid = [234 240]; % start and end times, in seconds

%%
pDataFilePath = ['/Users/hyang/Dropbox (HMS)/EphysAnalysis-Helen/pData/'...
    '200826_fly01_cell01_trial01.mat'];

% full path to leg video
% vidPath = '/Users/hyang/Desktop';
% vidName = 'fly01_cell01_trial01_legVid';

vidPath = '/Users/hyang/Dropbox (HMS)/APTlegTracking/200826_fly01_cell01_trial01_legVid_crop.mp4';

%
% full path to .trk leg tracking data
trkFilePath = ['/Users/hyang/Dropbox (HMS)/APTlegTracking/'...
    'Tracking_v1_mdn_labeled2798_ephys/' ...
    '200826_fly01_cell01_trial01_legVid_crop_caLegTrack_v1_mdn.trk'];

tVid = [50 56]; % start and end times, in seconds

%%



% indicies for specific body parts
% fit line to head and 3 pts defining midpt b/w legs
refPts.lineFitInd = 7:10; 
refPts.midPtInd = 9;
refPts.headPtInd = 7;
refPts.abdPtInd = 11;


% load pData
load(pDataFilePath, 'ephysData', 'ephysSpikes', 'fictrac', ...
    'fictracProc', 'fictracParams', 'leg');

% load .trk file; get leg positions
[legX, legY] = loadTrkFile(trkFilePath);


% some parameters about .trk file
numFrames = size(legX, 1); % number of frames in this trial
numPoints = size(legX, 2); % number of tracked points

% get leg positions after alignment to fly midpoint and normalization 
%  to body length units
[srnLegX, srnLegY] = shiftRotateNormalizeLegPos(legX, legY, refPts);



% FicTrac variables
ficTracTimes = fictrac.t;

sampRate = 1/(median(diff(ficTracTimes)));
avgWindow = 0.1;
smoYawAngVel = moveAvgFilt(fictrac.yawAngVel, sampRate, avgWindow);
smoFwdVel = moveAvgFilt(fictrac.fwdVel, sampRate, avgWindow);


% legVid variables
legVidFrameTimes = leg.frameTimes;

legStartInd = find(legVidFrameTimes >= tVid(1), 1, 'first');
legEndInd = find(legVidFrameTimes < tVid(2), 1, 'last');
legVidImgSize = [416 416];

% video reader
vidReader = VideoReader(vidPath);

legVidImg = read(vidReader,[legStartInd legEndInd]);

legVidImgGray = zeros(legVidImgSize(1),legVidImgSize(2),(legEndInd-legStartInd + 1));

% convert to grayscale
for i = 1:(legEndInd-legStartInd + 1)
    legVidImgGray(:,:,i) = rgb2gray(legVidImg(:,:,:,i));
end
    

legVidMinInt = 30;
legVidMaxInt = 160;

legVidFrameRate = 1/(median(diff(legVidFrameTimes)));

% get indicies to display on each frame
ficTracDispInd = zeros(1,legEndInd - legStartInd + 1);
for i = 1:(legEndInd - legStartInd + 1)
    ficTracDispInd(i) = find(ficTracTimes == ...
        legVidFrameTimes(i + legStartInd - 1), 1, 'first');    
end
ephysDispInd = ficTracDispInd;


ficTracStartInd = ficTracDispInd(1);
ficTracFigTimes = ficTracTimes(ficTracDispInd) - ficTracTimes(ficTracStartInd);
smoFwdVelFig = smoFwdVel(ficTracDispInd);
smoAngVelFig = smoYawAngVel(ficTracDispInd);
ephysVFig = ephysData.scaledVoltage;

fwdVelScale = [-5 10];
angVelScale = [-200 200];
ephysVScale = [-55 -25];
legXScale = [-1 1];
legYScale = [-1 1];

timeScale = [0 6];

% colormap for points
cmap = lines(size(legX,2));
% swap out last 4 as same color for body
bodyCol = cmap(1,:);
cmap(1:6,:) = cmap(2:7,:);
cmap(7:end,:) = repmat(bodyCol,5,1);

% make figure
figFrames(legEndInd - legStartInd + 1) = struct('cdata',[],'colormap',[]);
% figFrames(200) = struct('cdata',[],'colormap',[]);
f = figure;
set(f,'Position',[10,10,1000,700]);
for i = (legEndInd - legStartInd + 1)%1:(legEndInd - legStartInd + 1)
    ficTracInd = ficTracDispInd(i);
    legTrackInd = i+legStartInd - 1;
        
    % leg video
    legAx = subplot(5, 3, [1, 4, 7, 10, 13]);
    imagesc(legVidImgGray(:,:,i),[legVidMinInt legVidMaxInt]);
    hold on;
    scatter(legX(legTrackInd,:),legY(legTrackInd,:),50,cmap,'x', 'LineWidth', 1.5);
    axis equal;
    axis tight;
    colormap(legAx,'gray');
    set(gca,'XTick',[],'YTick',[]);
    legPos = get(legAx,'Position');
    newLegPos = [legPos(1)-0.06, legPos(2)-0.2, legPos(3)*1.5, legPos(4)*1.5];
    set(legAx,'Position',newLegPos);
    
    % ephys trace
    ephysTime = ephysData.t(ephysDispInd(1):ephysDispInd(i)) - ephysData.t(ephysDispInd(1));
    ephysAx = subplot(5, 3, [2,3]);
    plot(ephysTime, ephysVFig(ephysDispInd(1):ephysDispInd(i)), 'b');
    ylim(ephysVScale);
    xlim(timeScale);
    yticks([-55 -45 -35 -25]);
    ylabel('Membrane Potential (mV)');
    ephysPos = get(ephysAx,'Position');
    newEphysPos = [ephysPos(1)+0.06, ephysPos(2), ephysPos(3), ephysPos(4)*0.9];
    set(ephysAx,'Position',newEphysPos);
    title('Electrophysiology Recording');
    
    % leg time
    legTime = leg.frameTimes(legStartInd:legTrackInd) - leg.frameTimes(legStartInd);
    % leg X position
    legXAx = subplot(5,3,[5,6]);
    hold on;
    for j = 1:6
        plot(legTime, srnLegX(legStartInd:legTrackInd,j),'Color',cmap(j,:));
    end
    ylim(legXScale);
    xlim(timeScale);
    yticks([-1 -0.5 0 0.5 1]);
    ylabel('body lengths');
    legXPos = get(legXAx, 'Position');
    newLegXPos = [legXPos(1)+0.06, legXPos(2), legXPos(3), legXPos(4)*0.9];
    set(legXAx, 'Position', newLegXPos);
    set(legXAx, 'Ydir', 'reverse');
    title('Leg Position (Front-Back Axis)');
    
    
    % leg Y position
    legYAx = subplot(5,3,[8,9]);
    hold on;
    for j = 1:6
        plot(legTime, srnLegY(legStartInd:legTrackInd,j),'Color',cmap(j,:));    
    end
    ylim(legYScale);
    xlim(timeScale);
    yticks([-1 -0.5 0 0.5 1]);
    ylabel('body lengths');
    legYPos = get(legYAx, 'Position');
    newLegYPos = [legYPos(1)+0.06, legYPos(2), legYPos(3), legYPos(4)*0.9];
    set(legYAx, 'Position', newLegYPos);
    title('Leg Position (Left-Right Axis)');
      
    
    % FicTrac forward velocity
    fwdVelAx = subplot(5,3,[11,12]);
    plot(ficTracFigTimes(1:i),...
        smoFwdVelFig(1:i), 'k',...
        'LineWidth',1.5);
    ylim(fwdVelScale);
    xlim(timeScale);
    ylabel('mm/s');
    fwdVelPos = get(fwdVelAx,'Position');
    newFVPos = [fwdVelPos(1)+0.06, fwdVelPos(2), fwdVelPos(3), fwdVelPos(4)*0.9];
    set(fwdVelAx,'Position',newFVPos);   
    title('Forward Velocity');
    
    % FicTrac yaw velocity
    yawAx = subplot(5,3,[14,15]);
    plot(ficTracFigTimes(1:i),...
        smoAngVelFig(1:i), 'k',...
        'LineWidth',1.5);
    ylim(angVelScale);
    xlim(timeScale);
    xlabel('Time (s)');
    ylabel('deg/s');
    yawPos = get(yawAx,'Position');
    newYawPos = [yawPos(1)+0.06, yawPos(2), yawPos(3), yawPos(4)*0.9];
    set(yawAx,'Position',newYawPos);   
    title('Rotational Velocity');   
    
    drawnow;
    figFrames(i) = getframe(gcf); 
end

%% save video

vidPath = '/Users/hyang/Dropbox (Personal)/Wilson Lab/Conferences/U19 Review 2021';
vidName = [vidPath filesep 'figureVid2.mp4'];

v = VideoWriter(vidName, 'MPEG-4');
v.FrameRate = legVidFrameRate;
open(v);
writeVideo(v, figFrames);
close(v);