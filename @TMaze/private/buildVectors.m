function [thVector, zVector, dVector] = buildVectors(obj, trialIdx)

[~, zInd] = intersect(obj.dataTMaze.SESSION.allTrials(1).pospars, 'Z');
[~, thInd] = intersect(obj.dataTMaze.SESSION.allTrials(1).pospars, 'theta');

nSamples = 0;
for trialNum = 1:length(trialIdx)
    iTrial = trialIdx(trialNum);
    if ~isempty(obj.timesVRframes(iTrial).idx)
        nSamples = nSamples + length(obj.timesVRframes(iTrial).idx)-1;
    end
end

% pre-allocating
thVector = nan(nSamples, 1);
zVector = nan(nSamples, 1);
fVector = nan(nSamples, size(fData, 2));
tVector = nan(nSamples, 1);

nSamplesAccum = 0;
for trialNum = 1:length(trialIdx)
    iTrial = trialIdx(trialNum);
    tt = obj.timesVRframes(iTrial).t;
    if isempty(tt)
        continue;
    end
    idx = obj.timesVRframes(iTrial).idx(2:end);
    tt = tt(2:end-1);
    nSamplesThisTrial = length(tt);
    sampleIdx = nSamplesAccum+1:nSamplesAccum+nSamplesThisTrial;
    thVector(sampleIdx) = obj.dataTMaze.SESSION.allTrials(iTrial).posdata(idx, thInd);
    zVector(sampleIdx) = -obj.dataTMaze.SESSION.allTrials(iTrial).posdata(idx, zInd);
    fVector(sampleIdx, :) = interp1(tData, fData, tt);
    tVector(sampleIdx) = tt;
    nSamplesAccum = nSamplesAccum + nSamplesThisTrial;
end

thVector = thVector*180/pi; % transform from radians to degrees

end % buildVectors()
