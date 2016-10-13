function [rVector, vVector, fVector, tVector] = buildMotorVectors(obj, trialIdx, tData, fData)

% [~, zInd] = intersect(obj.dataTMaze.SESSION.allTrials(1).pospars, 'Z');
% [~, thInd] = intersect(obj.dataTMaze.SESSION.allTrials(1).pospars, 'theta');

rotGain = 1/8.3 * obj.dataTMaze.EXP.aGain * 180/pi;
velGain = -1/53 * obj.dataTMaze.EXP.zGain;

dtBall = median(diff(obj.dataBall.t));
dt2p = median(diff(tData));
dtRatio = (dt2p/dtBall);
% designing an anti-aliasing low-pass filter
N = ceil(dtRatio)*2+1;
[b, a] = fir1(N, 1/dtRatio, 'low');
r = filtfilt(b, a, rotGain/dtBall*obj.dataBall.rotation); % [deg/sec]
r = interp1(obj.dataBall.t, r, tData);
r = r*dtBall; % weird hack to scale properly, for some reason is only needed for the rotation!
v = filtfilt(b, a, velGain/dtBall*obj.dataBall.forward); % [cm/sec]
v = interp1(obj.dataBall.t, v, tData);


nTrials = length(trialIdx);
rCell = cell(nTrials, 1);
vCell = cell(nTrials, 1);
fCell = cell(nTrials, 1);
tCell = cell(nTrials, 1);
for trialNum = 1:nTrials
    iTrial = trialIdx(trialNum);
    tStart = obj.timesTrials.onset(iTrial);
    tEnd = obj.timesTrials.offset(iTrial);
    idx = tData>=tStart & tData<=tEnd;
    rCell{trialNum} = r(idx)';
    vCell{trialNum} = v(idx)';
    fCell{trialNum} = fData(idx);
    tCell{trialNum} = tData(idx)';

    % some plotting for debugging
%     zTrial = -obj.dataTMaze.SESSION.allTrials(iTrial).posdata(130:end, zInd);
%     thTrial = obj.dataTMaze.SESSION.allTrials(iTrial).posdata(130:end, thInd)*180/pi;
%     h1 = subplot(3,1,1);
%     plot(tData(idx)-tStart, (rCell{trialNum}), 'r', tData(idx)-tStart, (vCell{trialNum}), 'b')
%     title(['Trial ', num2str(iTrial)]);
%     legend('rotation', 'velocity');
%     h2 = subplot(3, 1, 2);
%     tPos = [1:length(thTrial)]/60;
%     plot(tPos, thTrial, 'r', tPos, zTrial, 'b');
%     legend('\theta', 'z');
%     h3 = subplot(3,1,3);
%     plot(tData(idx)-tStart, cumsum(rCell{trialNum})*dt2p, 'r', tData(idx)-tStart, cumsum(vCell{trialNum})*dt2p, 'b')
%     hold on;
%     plot(xlim, [30 30], 'k:', xlim, [-30 -30], 'k:');
%     hold off;
%     legend('CUM_{rot}', 'CUM_{vel}');
%     linkaxes([h1, h2, h3], 'x')
%     pause

end

rVector = cell2mat(rCell);
vVector = cell2mat(vCell);
fVector = cell2mat(fCell);
tVector = cell2mat(tCell);

end % buildMotorVectors()
