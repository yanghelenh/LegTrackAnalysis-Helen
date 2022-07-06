% callSwingStanceSteps.m
%
% Function that takes in steps (defined by start, mid, end indices),
%  step parameters (calculated in computeStepParameters()), and FicTrac
%  info to call each step's half step as swing or stance
% Option to call swing/stance based on half step duration: shorter duration
%  of two is swing, other is stance
% Option to call swing/stance based on FicTrac velocity: same direction is
%  stance, opposite direction is swing
% Select one of two options
%
% INPUTS:
%   legSteps - struct of leg step data, output from processLegTrack(), with
%       fields:
%       maxIndsAll
%       minIndsAll
%       maxWhichLeg
%       minWhichLeg
%       userSelVal
%       stepInds
%       stepWhichLeg
%       stepLengths
%       stepXLengths
%       stepYLengths
%       stepDirections
%       stepDurations
%       stepSpeeds
%       stepVelX
%       stepVelY
%       stepAEPX
%       stepAEPY
%       stepPEPX
%       stepPEPY
%       stepFtFwd
%       stepFtLat
%       stepFtYaw
%   legTrack - struct of leg tracking data, output of preprocessLegTrack
%   whichMethod - string for which method to determine swing/stance:
%       'duration' or 'fictrac'
%   fictrac - struct of fictrac data
%
% OUTPUTS:
%   legSteps - struct of leg step data, updated with additional fields:
%       stepSwingStance - n x 2 matrix for swing/stance calls for each half
%           step; 1 for stance, -1 for swing
%       swingStanceMethod - whichMethod from inputs
%
% CREATED: 10/1/21 - HHY
%
% UPDATED:
%   10/2/21 - HHY
%
function legSteps = callSwingStanceSteps(legSteps, legTrack, ...
    whichMethod, fictrac)

    % preallocate: n x 2
    stepSwingStance = zeros(size(legSteps.stepInds,1),2);
    
    switch whichMethod
        case 'duration'
            % loop through all steps, determine whether swing or stance
            for i = 1:size(legSteps.stepInds, 1)
                % compare durations of each half step; assign swing/stance
                %  appropriately: shorter is swing (-1)
                if (legSteps.stepDurations(i,1) >= ...
                        legSteps.stepDurations(i,2))
                    stepSwingStance(i,1) = 1;
                    stepSwingStance(i,2) = -1;
                else
                    stepSwingStance(i,1) = -1;
                    stepSwingStance(i,2) = 1;
                end
            end
        case 'fictrac'
            % interpolate FicTrac position to leg frame times (spline)
            interpFtFwdPos = interp1(fictrac.t, fictrac.fwdCumPos, ...
                legTrack.t, 'spline');

            % loop through all steps, determine whether swing or stance
            for i = 1:size(legSteps.stepInds, 1)

                % FicTrac fwd position for step start point
                stepStartFtFwdPos = interpFtFwdPos(legSteps.stepInds(i,1));
                % FicTrac fwd position for step mid point
                stepMidFtFwdPos = interpFtFwdPos(legSteps.stepInds(i,2));
                % FicTrac fwd position for step end point
                stepEndFtFwdPos = interpFtFwdPos(legSteps.stepInds(i,3));

                % get FicTrac velocity (delta fwd position / delta time), 
                %  binarized +1 for moving forward, -1 for moving backwards 
                thisStepFtFwdVel = (stepMidFtFwdPos - stepStartFtFwdPos)...
                    / legSteps.stepDurations(i,1);
                if (thisStepFtFwdVel >= 0)
                    thisStepBinFtFwdVel = 1; % fwd = 1
                else
                    thisStepBinFtFwdVel = -1; % not fwd = -1
                end

                % get binarized leg velocity
                thisStepLegXVel = legSteps.stepVelX(i,1);
                if (thisStepLegXVel >= 0)
                    thisStepBinLegXVel = 1; % fwd/front-to-back = 1
                else
                    thisStepBinLegXVel = -1; % not fwd/back-to-front = -1
                end

                % determine whether this half step is swing or stance by 
                %  multiplying together 2 binarized velocities
                % with how directions are defined, same direction 
                %  multiplies to 1 and is stance, opposite direction 
                %  multiplies to -1 and is swing
                stepSwingStance(i,1) = thisStepBinFtFwdVel * ...
                    thisStepBinLegXVel;


                % repeat for second half step

                % get FicTrac velocity (delta fwd position / delta time),
                %  binarized +1 for moving forward, -1 for moving backwards 
                thisStepFtFwdVel = (stepEndFtFwdPos - stepMidFtFwdPos) / ...
                    legSteps.stepDurations(i,2);
                if (thisStepFtFwdVel >= 0)
                    thisStepBinFtFwdVel = 1; % fwd = 1
                else
                    thisStepBinFtFwdVel = -1; % not fwd = -1
                end

                % get binarized leg velocity
                thisStepLegXVel = legSteps.stepVelX(i,2);
                if (thisStepLegXVel >= 0)
                    thisStepBinLegXVel = 1; % fwd/front-to-back = 1
                else
                    thisStepBinLegXVel = -1; % not fwd/back-to-front = -1
                end

                % determine whether this half step is swing or stance by 
                %  multiplying together 2 binarized velocities
                % with how directions are defined, same direction 
                %  multiplies to 1 and is stance, opposite direction 
                %  multiplies to -1 and is swing
                stepSwingStance(i,2) = thisStepBinFtFwdVel * ...
                    thisStepBinLegXVel;
            end
        otherwise
            disp('Invalid method for swing/stance call chosen');
            return;
    end


    % add swing/stance call to legSteps struct
    legSteps.stepSwingStance = stepSwingStance;
    legSteps.swingStanceMethod = whichMethod;
    
    
end