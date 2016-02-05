function [dThVector, dZVector, fVector, tVector] = buildDVectors(obj, trialIdx, tData, fData)

[~, zInd] = intersect(obj.dataTMaze.SESSION.allTrials.pospars, 'Z');
[~, thInd] = intersect(obj.dataTMaze.SESSION.allTrials.pospars, 'theta');

snippets = obj.dataTMaze.SESSION.replaySnippets(:);
snippetIdx = [snippets.trialIndex]';
snippets2Use = find(ismember(snippetIdx, trialIdx));
nSnippets = length(snippets2Use);

% we will apply slight LP filtering
order = 4;
b = hamming(order+1);
b = b./sum(b);

if isempty(snippets2Use)
    dThVector = [];
    dZVector = [];
    fVector = [];
    tVector = [];
else
    nSamplesPerSnippet = size(snippets(1).posdata, 1);
    posdata = cell2mat({snippets(snippets2Use).posdata}');
    
    zData = - reshape(posdata(:, zInd), nSamplesPerSnippet, nSnippets);
    thData = reshape(posdata(:, thInd), nSamplesPerSnippet, nSnippets);

    dt = mean(diff(obj.timesVRframes{1}));

    zData = imfilter(zData, b, 'replicate', 'same');
    thData = imfilter(thData, b, 'replicate', 'same');

    dZ = (zData(3:end, :) - zData(1:end-2, :))/(2*dt);
    dTheta = (thData(3:end, :) - thData(1:end-2, :))/(2*dt);
    
    tSamples = bsxfun(@plus, [1:nSamplesPerSnippet]', nSamplesPerSnippet*(snippets2Use(:)'-1));
    tSamples = tSamples(2:end-1, :);
    tSamples = tSamples(:);
    tt = obj.timesVRframes{1}(tSamples);

    dThVector = dTheta(:);
    dZVector = dZ(:);
    fVector = interp1(tData, fData, tt);
    tVector = tt;
    
    dThVector = dThVector*180/pi; % transform from radians to degrees
end

end % buildVectors()
