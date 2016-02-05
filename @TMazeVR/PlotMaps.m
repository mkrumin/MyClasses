function PlotMaps(obj, iPlane, iROI)

fData = obj.data2p{iPlane}.F(:, iROI);
fData = fData - min(fData);
fData = fData/max(fData)*63;
contrasts = unique(obj.dataTMaze.contrastSequence);
%     contrasts = [50 25 12 6 0 -50 -25 -12 -6];
contrastGroups = {contrasts(contrasts<0);...
    0;...
    contrasts(contrasts>0);
    contrasts};
nContrastGroups = length(contrastGroups);
contrastGroupLabels = {'StimL All';
    '0% contrast';
    'StimR All';
    'All Trials'};

figHandle = gcf;
clf(figHandle, 'reset');
set(figHandle, 'Name', sprintf('ROI # %d processing...', iROI));
drawnow;

trialIdx = 1:length(obj.dataTMaze.contrastSequence);
[thVector, zVector, ~] = buildVectors(obj, trialIdx, obj.times2p{iPlane}, fData);

zMax = max(abs(zVector));
zMin = min(abs(zVector));
thMax = max(abs(thVector));

nBinsTheta = 2*round(thMax/5)+1; % defining ~5deg large bins
nBinsZ = round(zMax-zMin)/5+1; %defining ~5cm large bins

zEdges = linspace(zMin-eps(zMin), zMax+eps(zMax), nBinsZ+1);
thEdges = linspace(-thMax-eps(thMax), thMax+eps(thMax), nBinsTheta+1);
dZ = median(diff(zEdges));
dTh = median(diff(thEdges));
thStd = 5;
zStd = 5;
hGauss = fspecial('gaussian', [ceil(7*zStd/dZ) 1], zStd/dZ)*fspecial('gaussian', [1 ceil(7*thStd/dTh)], thStd/dTh);
for iContrastGroup = 1:nContrastGroups
    % subdividing the trials
    trialIdx = find(ismember(obj.dataTMaze.contrastSequence, contrastGroups{iContrastGroup}));
    if isempty(trialIdx)
        continue;
    end
    [thVector, zVector, fVector] = buildVectors(obj, trialIdx, obj.times2p{iPlane}, fData);
    
    subplot(1, nContrastGroups, iContrastGroup);
    plot(thVector, zVector, '.k', 'MarkerSize', 4);
%     xlim([-thMax, thMax]+[-0.1 0.1]);
%     ylim([zMin-0.1 zMax+0.1]);
    title(contrastGroupLabels{iContrastGroup});
    hold on;

    [occMap{iContrastGroup}, occMapCenters] = hist3([zVector, thVector], 'Edges', {zEdges, thEdges});
    occMap{iContrastGroup} = occMap{iContrastGroup}(1:nBinsZ, 1:nBinsTheta);
    occMapCenters{1} = occMapCenters{1}(1:nBinsZ);
    occMapCenters{2} = occMapCenters{2}(1:nBinsTheta);
    occMap{iContrastGroup} = imfilter(occMap{iContrastGroup}, hGauss, 'replicate');
    fMap{iContrastGroup} = buildAccumMap([zVector, thVector], fVector, {zEdges, thEdges});
    fMap{iContrastGroup} = imfilter(fMap{iContrastGroup}, hGauss, 'replicate');
    finalMap{iContrastGroup} = fMap{iContrastGroup}./occMap{iContrastGroup};
end

minF = Inf;
maxF = -Inf;
for iContrastGroup = 1:nContrastGroups
    minF = min(minF, min(finalMap{iContrastGroup}(:)));
    maxF = max(maxF, max(finalMap{iContrastGroup}(:)));
end

for iContrastGroup = 1:nContrastGroups

%     subplot(4, nContrastGroups, iContrastGroup+nContrastGroups);
%     imagesc(occMapCenters{2}, occMapCenters{1}, log(1+occMap));
%     set(gca, 'YDir', 'Normal');
%     subplot(4, nContrastGroups, iContrastGroup+2*nContrastGroups);
%     imagesc(occMapCenters{2}, occMapCenters{1}, log(1+fMap));
%     set(gca, 'YDir', 'Normal');

    subplot(1, nContrastGroups, iContrastGroup);
    hold on;
    imagesc(occMapCenters{2}, occMapCenters{1}, finalMap{iContrastGroup}, 'AlphaData', 0.8);
    axis xy tight
    set(gca, 'CLim', [minF maxF]);
    box off;

%     subplot(2, nContrastGroups, iContrastGroup+nContrastGroups);
%     imagesc(occMapCenters{2}, occMapCenters{1}, finalMap{iContrastGroup}, 'AlphaData', 1);
%     set(gca, 'YDir', 'Normal');
%     set(gca, 'CLim', [minF maxF]);
%     box off;
end

set(figHandle, 'Name', sprintf('Cell # %d drawing...', iROI));
drawnow;

set(figHandle, 'Name', sprintf('Cell # %d done.', iROI));
end % PlotMaps()


%=========================================================================
function accMap = buildAccumMap(coords, signal, binEdges)

[nSamples, nDims] = size(coords);
indices = cell(nDims, 1);
nBins = nan(1, nDims);
% first, lets build indices matrices
for iDim = 1:nDims
    nBins(iDim) = length(binEdges{iDim})-1;
    indices{iDim} = nan(nSamples, nBins(iDim));
    for iBin = 1:nBins(iDim)
        indices{iDim}(:, iBin) = ...
            coords(:, iDim)>=binEdges{iDim}(iBin) & ...
            coords(:, iDim)<binEdges{iDim}(iBin+1);
    end
end
accMap = nan(nBins);

% trying to be independent of the nDims (not to hard-code the
% dimensionality)
str = '[sub(1)';
for iDim = 2:nDims
    str = sprintf('%s, sub(%d)', str, iDim);
end
str = [str, '] = ind2sub(nBins, iElement);'];

nElements = prod(nBins);
for iElement = 1:nElements
    eval(str);
    ind = ones(nSamples, 1);
    for iDim=1:length(sub);
        ind = ind & indices{iDim}(:, sub(iDim));
    end
    accMap(iElement) = sum(signal(ind));
end

end % buildAccumMap()