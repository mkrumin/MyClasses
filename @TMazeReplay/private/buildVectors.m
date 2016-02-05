function [thVector, zVector, fVector, tVector] = buildVectors(obj, trialIdx, tData, fData)

[~, zInd] = intersect(obj.dataTMaze.SESSION.allTrials.pospars, 'Z');
[~, thInd] = intersect(obj.dataTMaze.SESSION.allTrials.pospars, 'theta');

snippets = obj.dataTMaze.SESSION.replaySnippets(:);
snippetIdx = [snippets.trialIndex]';
snippets2Use = find(ismember(snippetIdx, trialIdx));

if isempty(snippets2Use)
    thVector = [];
    zVector = [];
    fVector = [];
    tVector = [];
else
    nSamplesPerSnippet = size(snippets(1).posdata, 1);
    posdata = cell2mat({snippets(snippets2Use).posdata}');
    tSamples = bsxfun(@plus, [1:nSamplesPerSnippet]', nSamplesPerSnippet*(snippets2Use(:)'-1));
    tSamples = tSamples(:);
    tt = obj.timesVRframes{1}(tSamples);
    
    thVector = posdata(:, thInd);
    zVector = -posdata(:, zInd);
    fVector = interp1(tData, fData, tt);
    tVector = tt;
    
    thVector = thVector*180/pi; % transform from radians to degrees
end

end % buildVectors()
