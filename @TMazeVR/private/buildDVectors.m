function [dThVector, dZVector, fVector, tVector] = buildDVectors(obj, trialIdx, tData, fData)

[~, zInd] = intersect(obj.dataTMaze.SESSION.allTrials(1).pospars, 'Z');
[~, thInd] = intersect(obj.dataTMaze.SESSION.allTrials(1).pospars, 'theta');

% pre-allocating
nTrials = length(trialIdx);
dThCell = cell(nTrials, 1);
dZCell = cell(nTrials, 1);
fCell = cell(nTrials, 1);
tCell = cell(nTrials, 1);

% we will apply slight LP filtering
order = 4;
b = hamming(order+1);
b = b./sum(b);

for trialNum = 1:length(trialIdx)
    iTrial = trialIdx(trialNum);
    tt = obj.timesVRframes(iTrial).t;
    if isempty(tt)
        continue;
    end
    idx = obj.timesVRframes(iTrial).idx(3:end);
    tt = tt(3:end-1);
    theta = obj.dataTMaze.SESSION.allTrials(iTrial).posdata(idx, thInd);
    z = -obj.dataTMaze.SESSION.allTrials(iTrial).posdata(idx, zInd);
    
    ttNew = linspace(tt(1), tt(end), length(tt));
    dt = mean(diff(ttNew)); 
    thetaNew = interp1(tt, theta, ttNew);
    zNew = interp1(tt, z, ttNew);
    
    dThCell{trialNum} = (thetaNew(3:end) - thetaNew(1:end-2))/(2*dt); % [rad/sec]
    dThCell{trialNum} = filtfilt(b, 1, dThCell{trialNum}(:));
    dZCell{trialNum} = (zNew(3:end) - zNew(1:end-2))/(2*dt); % [cm/sec]
    dZCell{trialNum} = filtfilt(b, 1, dZCell{trialNum}(:));
    ttNew = ttNew(2:end-1);
    fCell{trialNum} = interp1(tData, fData, ttNew(:));
    tCell{trialNum} = ttNew(:);

end

dThVector = cell2mat(dThCell);
dZVector = cell2mat(dZCell);
fVector = cell2mat(fCell);
tVector = cell2mat(tCell);

dThVector = dThVector*180/pi; % transform from [radians/sec] to [degrees/sec]

end % buildDVectors()
