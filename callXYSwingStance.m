% callXYSwingStance.m
% 
% Function to call swing, stance, and not moving based on comparing the
%  movement of the leg in the xy plane to the movement of the ball in the
%  xy plane. If the leg and the ball are moving in the same direction,
%  stance; leg and ball moving in opposite directions, swing. Not-moving
%  specified by findFlyNotMovingXY()
% Direction of movement of leg and ball computed over user-specified
%  window. Fit line to points, get vector
%
% INPUTS:
%   legXPos - matrix of leg x positions over time; num frames x num tracked
%       pts
%   legYPos - matrix of leg Y positions over time; num frames x num tracked
%       pts
%   legT - time for each leg position
%   zeroVelInd - indicies of leg time points when fly isn't moving
%   legInd - indicies of tracked legs
%   fwdPos - FicTrac fwdCumPos, forward cumulative position over time
%   slidePos - FicTrac slideCumPos, lateral cumulative position over time
%   yawAngPos - FicTrac yawAngCumPos, yaw angular cumulative position over
%       time
%   fictracT - time for each FicTrac position
%   winLen - length in seconds of window to compute leg and ball direction
%       over
%
% OUTPUTS:
%   legStance - logical for whether leg is in stance, for all legs; 
%       frames x legs matrix
%   legSwing - logical for whether leg is in swing, for all legs; 
%       frames x legs matrix
%   legSwingStanceNotMove - matrix for swing (-1), stance (1), not moving
%       (0)
%
% CREATED: 12/7/20 - HHY
%
% UPDATED:
%   12/7/20 - HHY
%
function [legStance, legSwing, legSwingStanceNotMove] = ...
    callXYSwingStance(legXPos, legYPos, legT, zeroVelInd, legInd, ...
    fwdPos, slidePos, yawAngPos, fictracT, winLen)

    % determine number of samples to compute direction over, for legs and
    %  FicTrac
    legFrameRate = 1/median(diff(legT));
    ftFrameRate = 1/median(diff(fictracT));
    
    numLegSamps = round(winLen * legFrameRate);
    % if even number of samples, make odd number
    if ~(mod(numlegSamps,2))
        numLegSamps = numLegSamps + 1;
    end
    
    numFtSamps = round(winLen * ftFrameRate);
    if ~(mod(numFtSamps,2))
        numFtSamps = numFtSamps + 1;
    end
    

    % loop over all legs
    for i = 1:legnth(legInd)
        % loop over all leg tracking time points
        for j = 1:size(legXPos,1)
            % get vector for leg movement direction
            % window brackets index j +/- (numLegSamps-1)/2 
            legStartInd =  j - (numLegSamps - 1)/2;
            if (legStartInd < 1) % at beginning, clip to start at 1
                legStartInd = 1;
            end
            legEndInd = j + (numLegSamps - 1)/2;
            if (legEndInd > size(legXPos,1)) % at end, clip to end ind
                legEndInd = size(legXPos,1);
            end
            
            % x and y coordinates of current window
            xCord = legXPos(legStartInd:legEndInd, legInd(i));
            yCord = legYPos(legStartInd:legEndInd, legInd(i));
            
            % fit line
            lineCoeffs = polyfit(xCord, yCord, 1);
            
            % project first and last point to line
            [projx1, projy1] = projPt2Line(lineCoeffs(1), lineCoeffs(2), ...
                xCord(1), yCord(1));
            [projx2, projy2] = projPt2Line(lineCoeffs(1), lineCoeffs(2), ...
                xCord(end), yCord(end));
            % zero vector
            legXVec = projx2 - projx1;
            legYVec = projy2 - projy1;
            
            % get vector for FicTrac movement direction
            midTime = legT(j);
            repST = repmat(midTime, size(fictracT));
            [~, midInd] = min(abs(fictracT - repST));
            ftStartInd = midInd - (numFtSamps - 1)/2;
            if (ftStartInd < 1)
                ftStartInd = 1;
            end
            ftEndInd = midInd - (numFtSamps - 1)/2;
            if (ftEndInd > size(fictracT,1))
                ftEndInd = size(fictracT,1);
            end
            
            % fwd and slide coordinates of fictrac of current window
            fwdCord = fwdPos(ftStartInd:ftEndInd);
            slideCord = slidePos(ftStartInd:ftEndInd);
            
            % fit line
            ftLineCoeffs = polyfit(fwdCord, slideCord, 1);
            
            % project first and last point to line
            [projf1, projs1] = projPt2Line(ftLineCoeffs(1), ...
                ftLineCoeffs(2), fwdCord(1), slideCord(1));
            [projf2, projs2] = projPt2Line(ftLineCoeffs(1), ...
                ftLineCoeffs(2), fwdCord(end), slideCord(end));
            % zero vector
            fwdVec = projf2 - projf1;
            slideVec = projs2 - projs1;
            
            % rotate vector by amount of rotation in current window
            degRot = yawAngPos(ftEndInd) - yawAngPos(ftStartInd);
            [fwdVecRot, slideVecRot] = rotatePts2D(fwdVec, slideVec, ...
                degRot);
          
            
        end
    end
    
    % indicate when fly isn't moving

end