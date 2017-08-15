function [thVector, zVector, fVector, tVector] = buildVectors(obj, trialIdx, tData, fData, isReplay)

%builds vector for defined trial, if trial is scrambled then unscrambles
%only takes argument of single trial idx not collection of trial idxs (yet)

[~, zInd] = intersect(obj.dataTMaze.SESSION.allTrials(1).pospars, 'Z');
[~, thInd] = intersect(obj.dataTMaze.SESSION.allTrials(1).pospars, 'theta');

if nargin < 5
    isReplay = false;
end 

%calculate number of frames in trial
if ~isReplay
    posdata = obj.dataTMaze.SESSION.allTrials(trialIdx).posdata;
    
    tt = obj.timesVRframes(trialIdx).t;
    tt = tt(2:end-1);
elseif isReplay
    snippets = obj.dataTMaze.SESSION.replaySnippets(:);
    snippetIdx = [snippets.trialIndex]';
    snippets2Use = find(ismember(snippetIdx, trialIdx));
    nSamplesPerSnippet = size(snippets(1).posdata, 1);
    
    posdata = cell2mat({snippets(snippets2Use).posdata}');
    
    tSamples = bsxfun(@plus, [1:nSamplesPerSnippet]', nSamplesPerSnippet*(snippets2Use(:)'-1));
    tSamples = tSamples(:);
    tt = obj.timesVRframes{1}(tSamples);
end

idx = obj.timesVRframes(trialIdx).idx(2:end); %why is first one always 1?

thVector = posdata(idx, thInd);
zVector = -posdata(idx, zInd);

fVector = interp1(tData, fData, tt);
tVector = tt;

thVector = thVector*180/pi; % transform from radians to degrees

end 
