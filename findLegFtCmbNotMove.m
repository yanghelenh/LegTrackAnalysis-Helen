% findLegFtCmbNotMove.m
% 
% Function that takes in start and end indices for not moving bouts, as
%  calculated on leg data and FicTrac data, and returns start and end
%  indices for not moving bouts on both datasets, as calculated as
%  intersect, union, leg only, or FicTrac only
%
% Helper function for interactGetNotMovingIndWFt()
%
% INPUTS:
%   legStartInd - start indices of not moving bouts, for leg data
%   legEndInd - end indices for not moving bouts, for leg data
%   legT - time vector for leg data, in sec
%   ftStartInd - start indices of not moving bouts, for FicTrac data
%   ftEndInd - end indices for not moving bouts, for FicTrac data
%   ftT - time vector for FicTrac data, in sec
%   whichMethod - which method for computing combined not moving bouts. 
%       Options are: 'intersect', 'union', 'legOnly', 'FictracOnly'
%
% OUTPUTS:
%   legNotMoveIndNew - indices of all not move points, for leg data
%   legStartIndNew - updated start indices of not moving bouts, for leg
%   legEndIndNew - updated end indices for not moving bouts, for leg data
%   ftNotMoveIndNew - indices of all not move points, for FicTrac data
%   ftStartInd - updated start indices of not moving bouts, for FicTrac
%   ftEndInd - updated end indices of not moving bouts, for FicTrac data
%
% CREATED: 7/1/22 - HHY
%
% UPDATED:
%   7/1/22 - HHY
%
function [legNotMoveIndNew, legStartIndNew, legEndIndNew, ...
    ftNotMoveIndNew, ftStartIndNew, ftEndIndNew] = findLegFtCmbNotMove(...
    legStartInd, legEndInd, legT, ftStartInd, ftEndInd, ftT, whichMethod)

    % initialize not move logicals for leg and fictrac
    legNotMoveLog = false(size(legT));
    ftNotMoveLog = false(size(ftT));

    % not move bout start and end times
    legNotMoveStartTimes = legT(legStartInd);
    legNotMoveEndTimes = legT(legEndInd);
    ftNotMoveStartTimes = ftT(ftStartInd);
    ftNotMoveEndTimes = ftT(ftEndInd);

    switch whichMethod
        % loop through all indices for leg and fictrac, save indices for 
        %  when both not moving
        case 'intersect'
            % leg loop
            for i = 1:length(legNotMoveLog)
                thisT = legT(i);
                % get index of not move bout that this point could belong
                %  to, for leg
                thisLegBoutInd = find(thisT >= legNotMoveStartTimes,...
                    1,'first'); 
                % if time is less than end time of bout, then in not moving
                %  bout
                if (thisT <= legNotMoveEndTimes(thisLegBoutInd))
                    thisLegNotMove = true;
                else
                    thisLegNotMove = false;
                    legNotMoveLog(i) = false;
                    continue; % don't need to check farther for this point
                end

                % get index of not move bout this point could belong to,
                %  for FicTrac
                thisFtBoutInd = find(thisT >= ftNotMoveStartTimes, ...
                    1, 'first');
                % if time is less than end time of bout, then in not moving
                %  bout
                if (thisT <= ftNotMoveEndTimes(thisFtBoutInd))
                    thisFtNotMove = true;
                else
                    thisFtNotMove = false;
                    legNotMoveLog(i) = false;
                    continue; % don't need to check farther
                end

                % if both not moving, then this also not moving
                if (thisLegNotMove && thisFitNotMove)
                    legNotMoveLog(i) = true;
                else
                    legNotMoveLog(i) = false;
                end
            end

            % FicTrac loop
            for i = 1:length(ftNotMoveLog)
                thisT = ftT(i);
                % get index of not move bout that this point could belong
                %  to, for leg
                thisLegBoutInd = find(thisT >= legNotMoveStartTimes,...
                    1,'first'); 
                % if time is less than end time of bout, then in not moving
                %  bout
                if (thisT <= legNotMoveEndTimes(thisLegBoutInd))
                    thisLegNotMove = true;
                else
                    thisLegNotMove = false;
                    ftNotMoveLog(i) = false;
                    continue; % don't need to check farther for this point
                end

                % get index of not move bout this point could belong to,
                %  for FicTrac
                thisFtBoutInd = find(thisT >= ftNotMoveStartTimes, ...
                    1, 'first');
                % if time is less than end time of bout, then in not moving
                %  bout
                if (thisT <= ftNotMoveEndTimes(thisFtBoutInd))
                    thisFtNotMove = true;
                else
                    thisFtNotMove = false;
                    ftNotMoveLog(i) = false;
                    continue; % don't need to check farther
                end

                % if both not moving, then this also not moving
                if (thisLegNotMove && thisFitNotMove)
                    ftNotMoveLog(i) = true;
                else
                    ftNotMoveLog(i) = false;
                end
            end

        % loop through all not moving bouts for both leg and fictrac, save 
        %  indices for when either is not moving
        case 'union'
            % leg loop
            for i = 1:length(legNotMoveStartTimes)
                % for leg logical, flip leg bits
                legNotMoveLog(legStartInd(i):legEndInd(i)) = true;

                % for FicTrac logical, get indices corresponding to this
                %  bout's start and end times
                thisStartInd = find(ftT >= legNotMoveStartTimes(i),1,...
                    'first');
                thisEndInd = find(ftT <= legNotMoveEndTimes(i),1,'last');
                % flip appropriate bits
                ftNotMoveLog(thisStartInd:thisEndInd) = true;
            end

            % FicTrac loop
            for i = 1:length(ftNotMoveStartTimes)
                % for leg logical, get indices corresponding to this bout's
                %  start and end times
                thisStartInd = find(legT >= ftNotMoveStartTimes(i),1,'first');
                thisEndInd = find(legT <= ftNotMoveEndTimes(i),1,'last');
                % flip bits
                legNotMoveLog(thisStartInd:thisEndInd) = true;

                % for FicTrac logical, flip FicTrac bits
                ftNotMoveLog(ftStartInd(i):ftEndInd(i)) = true;
            end

        % convert leg not moving bouts to FicTrac time base
        case 'legOnly'
            % loop through all leg not moving bouts
            for i = 1:length(legNotMoveStartTimes)
                % for leg logical, flip bits
                legNotMoveLog(legStartInd(i):legEndInd(i)) = true;

                % for FicTrac logical, get indices corresponding to this
                %  bout's start and end times
                thisStartInd = find(ftT >= legNotMoveStartTimes(i),1,...
                    'first');
                thisEndInd = find(ftT <= legNotMoveEndTimes(i),1,'last');
                % flip appropriate bits
                ftNotMoveLog(thisStartInd:thisEndInd) = true;
            end

        % convert FicTrac not moving bouts to leg time base
        case 'fictracOnly'
            % loop through all FicTrac not moving bouts
            for i = 1:length(ftNotMoveStartTimes)
                % for leg logical, get indices corresponding to this bout's
                %  start and end times
                thisStartInd = find(legT >= ftNotMoveStartTimes(i),1,'first');
                thisEndInd = find(legT <= ftNotMoveEndTimes(i),1,'last');
                % flip bits
                legNotMoveLog(thisStartInd:thisEndInd) = true;

                % for FicTrac logical, flip FicTrac bits
                ftNotMoveLog(ftStartInd(i):ftEndInd(i)) = true;
            end

        % return empty vectors if invalid method
        otherwise
            disp('Invalid method for combining Leg and FicTrac bouts');
            legNotMoveIndNew = [];
            legStartIndNew = [];
            legEndIndNew = [];
            ftNotMoveIndNew = [];
            ftStartInd = [];
            ftEndInd = [];
            return; % end function here
    end

    % convert logicals to start and end indices of bouts, and all indices
    [legStartIndNew, legEndIndNew, legNotMoveIndNew] = ...
        convertNotMoveLogToBouts(legNotMoveLog);
    [ftStartIndNew, ftEndIndNew, ftNotMoveIndNew] = ...
        convertNotMoveLogToBouts(ftNotMoveLog);

    if (isrow(legNotMoveIndNew))
        legNotMoveIndNew = legNotMoveIndNew';
    end
    if (isrow(legStartIndNew))
        legStartIndNew = legStartIndNew';
    end
    if (isrow(legEndIndNew))
        legEndIndNew = legEndIndNew';
    end
    if (isrow(ftNotMoveIndNew))
        ftNotMoveIndNew = ftNotMoveIndNew';
    end
    if (isrow(ftStartIndNew))
        ftStartIndNew = ftStartIndNew';
    end
    if (isrow(ftEndIndNew))
        ftEndIndNew = ftEndIndNew';
    end
end