function [thMatrix, zMatrix, fMatrix] = buildTraces(obj, trialIdx, tData, fData, nPoints)

if nargin<5 || isempty(nPoints)
    nPoints = 100;
end

[~, zInd] = intersect(obj.dataTMaze.SESSION.allTrials(1).pospars, 'Z');
[~, thInd] = intersect(obj.dataTMaze.SESSION.allTrials(1).pospars, 'theta');


% pre-allocating
nTrials = length(trialIdx);
thMatrix = nan(nTrials, nPoints);
zMatrix = nan(nTrials, nPoints);
fMatrix = nan(nTrials, nPoints);

% nSamplesAccum = 0;
for trialNum = 1:nTrials
    iTrial = trialIdx(trialNum);
    tt = obj.timesVRframes(iTrial).t;
    if isempty(tt)
        continue;
    end
    idx = [obj.timesVRframes(iTrial).idx(2:end); obj.timesVRframes(iTrial).idx(end)+1];
    tt = tt(2:end);
    
    % correction to exclude the long static frame at the end of the trial
    % (when the reward is consumed and the next trial is being prepared)
    tt = tt(1:end-1);
    idx = idx(1:end-1);
    
    tAxis = linspace(tt(1), tt(end), nPoints);
%     nSamplesThisTrial = length(tt);
%     sampleIdx = nSamplesAccum+1:nSamplesAccum+nSamplesThisTrial;
    thMatrix(trialNum, :) = interp1(tt, obj.dataTMaze.SESSION.allTrials(iTrial).posdata(idx, thInd), tAxis);
    zMatrix(trialNum, :) = interp1(tt, -obj.dataTMaze.SESSION.allTrials(iTrial).posdata(idx, zInd), tAxis);
    fMatrix(trialNum, :) = interp1(tData, fData, tAxis);
%     nSamplesAccum = nSamplesAccum + nSamplesThisTrial;
end

thMatrix = thMatrix*180/pi; % transform from radians to degrees

end % buildVectors()
