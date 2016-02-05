function obj = getPosUniformZ(obj, nSamples)

if nargin<2
    nSamples = 100;
end

trials = obj.SESSION.allTrials(1:obj.nTrials);
nTrials = length(trials);
[~, zInd] = intersect(trials(1).pospars, 'Z');
[~, xInd] = intersect(trials(1).pospars, 'X');
[~, thInd] = intersect(trials(1).pospars, 'theta');
[~, spInd] = intersect(trials(1).pospars, 'speed');

%% ---- Z-based and interpolated single trajectories approach ----
% Z coordinate is uniformly sampled (and interpolated if needed) at nSamples points, 
% and we get x(z), theta(z), and speed(z) vectors of length nSamples
zMin = Inf;
zMax = -Inf;
for iTrial = 1:nTrials
    zMin = min(zMin, min(-trials(iTrial).posdata(:, zInd)));
    zMax = max(zMax, max(-trials(iTrial).posdata(:, zInd)));
end

zAxis = linspace(zMin, zMax, nSamples)';
zUniform = repmat(zAxis, 1, nTrials);
xUniform = nan(nSamples, nTrials);
thUniform = nan(nSamples, nTrials);
spUniform = nan(nSamples, nTrials);
for iTrial = 1:nTrials
    zVector = -trials(iTrial).posdata(:, zInd);
    xVector = trials(iTrial).posdata(:, xInd);
    thVector = trials(iTrial).posdata(:, thInd);
    spVector = trials(iTrial).posdata(:, spInd);
    [counts, ind] = histc(zVector, zAxis);
    goodIdx = find(counts~=0);
    emptyIdx = find(counts == 0);
    for iBin = 1:length(goodIdx)
%         zUniform(goodIdx(iBin), iTrial) = mean(zVector(ind == goodIdx(iBin)));
        xUniform(goodIdx(iBin), iTrial) = mean(xVector(ind == goodIdx(iBin)));
        thUniform(goodIdx(iBin), iTrial) = mean(thVector(ind == goodIdx(iBin)));
        spUniform(goodIdx(iBin), iTrial) = mean(spVector(ind == goodIdx(iBin)));
    end
    
    % interpolating for bins without samples;
    % extrapolation points will be NaNs
    if ~isempty(emptyIdx) && length(goodIdx)>1
%         zUniform(emptyIdx, iTrial) = interp1(zAxis(goodIdx), zUniform(goodIdx, iTrial), zAxis(emptyIdx));
        xUniform(emptyIdx, iTrial) = interp1(zAxis(goodIdx), xUniform(goodIdx, iTrial), zAxis(emptyIdx));
        thUniform(emptyIdx, iTrial) = interp1(zAxis(goodIdx), thUniform(goodIdx, iTrial), zAxis(emptyIdx));
        spUniform(emptyIdx, iTrial) = interp1(zAxis(goodIdx), spUniform(goodIdx, iTrial), zAxis(emptyIdx));
        
        xNans = isnan(xUniform(:, iTrial));
        xUniform(xNans, iTrial) = interp1(zAxis(~xNans), xUniform(~xNans, iTrial), zAxis(xNans), 'nearest', 'extrap');
        thNans = isnan(thUniform(:, iTrial));
        thUniform(thNans, iTrial) = interp1(zAxis(~thNans), thUniform(~thNans, iTrial), zAxis(thNans), 'nearest', 'extrap');
        spNans = isnan(spUniform(:, iTrial));
        spUniform(spNans, iTrial) = interp1(zAxis(~spNans), spUniform(~spNans, iTrial), zAxis(spNans), 'nearest', 'extrap');

    end
    
    
end

obj.posUniform.z = zUniform;
obj.posUniform.x = xUniform;
obj.posUniform.th = thUniform;
obj.posUniform.sp = spUniform;

end % getPosUniformZ()


%% ---- dots' cloud approach (doesn't work well) ------
% pool all the data into the same 'cloud' of dots
%     zAll = [];
%     thAll = [];
%     for iTrial = 1:nTrials
%         zStart = trials(iTrial).posdata(1, zInd);
%         % find the last frame of the start position in the maze
%         indStart = find(trials(iTrial).posdata(:, zInd) ~= zStart, 1, 'first') - 1;
%         zVector = trials(iTrial).posdata(indStart:end, zInd);
%         thVector = trials(iTrial).posdata(indStart:end, thInd);
%         zAll = cat(1, zAll, zVector);
%         thAll = cat(1, thAll, thVector);
%     end
%     zAll = -zAll; % make them positive - OpenGL feature to make them negative
%     [zAll, sortedIdx] = sort(zAll);
%     thAll = thAll(sortedIdx);
%
%     zMin = zAll(1);
%     zMax = zAll(end);
%
%     % let's exclude the points of traveling along the back and front walls
%     goodIdx = zAll~=zMin & zAll~=zMax;
%     zAll = zAll(goodIdx);
%     thAll = thAll(goodIdx);
%
%     binEdges = prctile(zAll, 0:1:nSamples);
% %     binEdges = linspace(min(zAll), max(zAll), nSamples+1);
%     nBins = length(binEdges)-1;
%
%     zMean = nan(nBins, 1);
%     thMean = nan(nBins, 1);
%     for iBin = 1:nBins
%         idx = zAll>binEdges(iBin) & zAll<=binEdges(iBin+1);
%         zVector = zAll(idx);
%         thVector = thAll(idx);
%         zMean(iBin) = median(zVector);
%         thMean(iBin) = median(thVector);
%     end

