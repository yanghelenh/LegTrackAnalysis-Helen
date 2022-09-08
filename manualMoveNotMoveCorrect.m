% manualMoveNotMoveCorrect.m
% quick script to modify moveNotMove to make not moving anything after
% cutoffT (in sec)
% 8/5/22

cutoffT = 120;
cutoffLegInd = find(legTrack.t >= cutoffT,1,'first');
cutoffFtInd = find(fictracProc.t >= cutoffT, 1, 'first');

legNotMoveInd = unique([moveNotMove.legNotMoveInd; (cutoffLegInd:length(legTrack.t))']);

legNotMoveLog = false(size(legTrack.t'));
legNotMoveLog(legNotMoveInd) = true;

    [legStartIndNew, legEndIndNew, legNotMoveIndNew] = ...
        convertNotMoveLogToBouts(legNotMoveLog);

    legNotMoveStartInd = legStartIndNew';
    legNotMoveEndInd = legEndIndNew';

ftNotMoveInd = unique([moveNotMove.ftNotMoveInd; (cutoffFtInd:length(fictracProc.t))']);

ftNotMoveLog = false(size(fictracProc.t));
ftNotMoveLog(ftNotMoveInd) = true;


    [ftStartIndNew, ftEndIndNew, ftNotMoveIndNew] = ...
        convertNotMoveLogToBouts(ftNotMoveLog);

    ftNotMoveStartInd = ftStartIndNew';
    ftNotMoveEndInd = ftEndIndNew';


        % get all return values, leg
    legNotMoveBout = [legNotMoveStartInd, legNotMoveEndInd];
    
    allLegInd = (1:length(legTrack.t))'; % all indicies in trial
    
    % move indicies are inverse of not move indicies
    legMoveInd = allLegInd(~ismember(allLegInd,legNotMoveInd));
    % convert indicies to bout starts and ends
    [legMoveBoutStartInd, legMoveBoutEndInd, ~] = findBouts(legMoveInd);
    legMoveBout = [legMoveBoutStartInd, legMoveBoutEndInd];

    % get all return values, FicTrac
    ftNotMoveBout = [ftNotMoveStartInd, ftNotMoveEndInd];

    allFtInd = (1:length(fictracProc.t))';

    ftMoveInd = allFtInd(~ismember(allFtInd,ftNotMoveInd));
    [ftMoveBoutStartInd, ftMoveBoutEndInd, ~] = findBouts(ftMoveInd);
    ftMoveBout = [ftMoveBoutStartInd, ftMoveBoutEndInd];

            moveNotMove.legNotMoveInd = legNotMoveInd;
        moveNotMove.legNotMoveBout = legNotMoveBout;
        moveNotMove.legMoveInd = legMoveInd;
        moveNotMove.legMoveBout = legMoveBout;
        moveNotMove.legT = legTrack.t;
        moveNotMove.ftNotMoveInd = ftNotMoveInd;
        moveNotMove.ftNotMoveBout = ftNotMoveBout;
        moveNotMove.ftMoveInd = ftMoveInd;
        moveNotMove.ftMoveBout = ftMoveBout;
        moveNotMove.ftT = fictracProc.t;
        moveNotMove.notMoveParams = moveNotMove.notMoveParams;