function [rVector, vVector, fVector, tVector] = buildMotorVectors(obj, trialIdx, tData, fData)

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

snippets = obj.dataTMaze.SESSION.replaySnippets(:);
snippetIdx = [snippets.trialIndex]';
snippets2Use = find(ismember(snippetIdx, trialIdx));
ourSnippets = snippets(snippets2Use); % now we are only left with snippets of interest

nSnippets = length(ourSnippets);
rCell = cell(nSnippets, 1);
vCell = cell(nSnippets, 1);
fCell = cell(nSnippets, 1);
tCell = cell(nSnippets, 1);
nSamplesPerSnippet = size(snippets(1).posdata, 1);
idxStart = nSamplesPerSnippet*(snippets2Use-1)+1;
idxEnd = idxStart+nSamplesPerSnippet-1;
dtVR = median(diff(obj.timesVRframes{1}));

tStart = obj.timesVRframes{1}(idxStart);
tEnd = obj.timesVRframes{1}(idxEnd) + dtVR;

for iSnippet = 1:nSnippets
    
    idx = find(tData>=tStart(iSnippet) & tData<=tEnd(iSnippet)); % using find() is faster here than binary index
    rCell{iSnippet} = r(idx)';
    vCell{iSnippet} = v(idx)';
    fCell{iSnippet} = fData(idx);
    tCell{iSnippet} = tData(idx)';

end

rVector = cell2mat(rCell);
vVector = cell2mat(vCell);
fVector = cell2mat(fCell);
tVector = cell2mat(tCell);

end % buildMotorVectors()
