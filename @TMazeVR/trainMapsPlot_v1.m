function trainMap(obj, iPlane, iROI, options)

if nargin<4
    options = struct;
end
options = fillOptions(obj, options);

nGroups = length(options.contrastGroups);
nChunks = options.cvFactor;

fData = obj.data2p{iPlane}.F(:, iROI);
tData = obj.times2p{iPlane};

figHandle = gcf;
clf(figHandle, 'reset');
set(figHandle, 'Name', sprintf('ROI # %d processing...', iROI));
drawnow;
set(figHandle, 'Units', 'normalized', 'MenuBar', 'none', 'OuterPosition', [0 0 1 1]);

trialIdx = 1:length(obj.dataTMaze.contrastSequence);
[thAll, zAll, ~] = buildVectors(obj, trialIdx, tData, fData);

zMax = max(abs(zAll));
zMin = min(abs(zAll));
thMax = max(abs(thAll));
clear zAll thAll;

% now, a complicated way to define bin edges
thEdges = buildEdges(thMax, options.dTheta);
zEdges = buildEdges([zMin, zMax], options.dZ);

% fMap = cell(nContrastGroups, options.cvFactor);
thVector = cell(nGroups, nChunks);
zVector = cell(nGroups, nChunks);
fVector = cell(nGroups, nChunks);
traces = cell(nGroups, 1);
for iGroup = 1:nGroups
    % picking the correct trials
    trialIdx = find(ismember(obj.dataTMaze.contrastSequence, options.contrastGroups{iGroup}));
    if isempty(trialIdx)
        continue;
    end
    nTrials = length(trialIdx);
    nTrialsPerSet = ceil(nTrials/nChunks);
    randIdx = nan(nChunks, nTrialsPerSet);
    randIdx(1:nTrials) = randperm(nTrials);
    for iChunk = 1:nChunks
        % extracting part of the data
        idx = randIdx(iChunk, :);
        idx = idx (~isnan(idx));
        [thVector{iGroup, iChunk}, zVector{iGroup, iChunk}, fVector{iGroup, iChunk}] = buildVectors(obj, trialIdx(idx), tData, fData);
    end
    
    traces{iGroup} = getTimeScaledETA(tData, fData, ...
         [obj.timesTmazeOnset(trialIdx), obj.timesTmazeOffset(trialIdx)]);

end

optStd = nan(nGroups, 2);
errVals = nan(nGroups, 1);
coords = cell(nGroups, nChunks);
errValMatrix = cell(nGroups, 1);
stdGridValues = cell(nGroups, 1);
% estimating the optimal filter STDs
validGroups = 1:nGroups;
for iGroup = 1:nGroups
    for iChunk = 1:nChunks
        coords{iGroup, iChunk} = [zVector{iGroup, iChunk}, thVector{iGroup, iChunk}];
    end
    if isempty(coords{iGroup, 1})
        validGroups = setdiff(validGroups, iGroup);
        continue;
    end
    [optStd(iGroup, :), errVals(iGroup), errValMatrix{iGroup}, stdGridValues{iGroup}] = ...
        estOptimalStd(coords(iGroup, :), fVector(iGroup, :)', {zEdges, thEdges});
end


% calculating the optimal maps (chunk-by-chunk)
% optStd = abs(optStd); % HACK! Unnecessary for the grid search
fGroupMap = cell(nGroups, 1);
rawGroupMap = cell(nGroups, 1);
for iGroup = validGroups
    hGauss = ndGaussian(optStd(iGroup, :));
    [fGroupMap{iGroup}, binCentres] = estMap([cell2mat(zVector(iGroup, :)'), cell2mat(thVector(iGroup, :)')], ...
        cell2mat(fVector(iGroup, :)'), {zEdges, thEdges}, hGauss);
    rawGroupMap{iGroup} = estMap([cell2mat(zVector(iGroup, :)'), cell2mat(thVector(iGroup, :)')], ...
        cell2mat(fVector(iGroup, :)'), {zEdges, thEdges}, 1, 0);
end

% calculating the range of the colormaps
minF = Inf;
maxF = -Inf;
for i = validGroups
    minF = min(minF, min(fGroupMap{i}(:)));
    maxF = max(maxF, max(fGroupMap{i}(:)));
end
minRawF = Inf;
maxRawF = -Inf;
for i = validGroups
    minRawF = min(minRawF, min(rawGroupMap{i}(:)));
    maxRawF = max(maxRawF, max(rawGroupMap{i}(:)));
end

nRows = 5;
nColumns = nGroups + 2;
dZ = mean(diff(binCentres{1}));
dTh = mean(diff(binCentres{2}));
for iGroup = validGroups%1:nGroups

    subplot(nRows, nColumns, iGroup);
    set(gca, 'Position', get(gca, 'Position'));
    plot(cell2mat(thVector(iGroup, :)'), cell2mat(zVector(iGroup, :)'), '.k', 'MarkerSize', 4);
    hold on;
    imagesc(binCentres{2}, binCentres{1}, rawGroupMap{iGroup}, 'AlphaData', 0.8);
    axis xy tight off
    set(gca, 'CLim', [minRawF maxRawF]);
    box off;
    title(options.contrastGroupLabels{iGroup});
    if iGroup == 1
        tmp = xlim;
        text(tmp(1)-0.2*diff(tmp), mean(ylim), '''Raw'' Data', 'Rotation', 90, ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
    end

    subplot(nRows, nColumns, nColumns+iGroup);
    set(gca, 'Position', get(gca, 'Position'));
    surf(stdGridValues{iGroup}{2}*dTh, stdGridValues{iGroup}{1}*dZ, errValMatrix{iGroup});
    view([0 90]);
    axis tight off;
    box off;
    title(sprintf('err: %.2f', errVals(iGroup)));
    if iGroup == 1
        tmp = xlim;
        text(tmp(1)-0.2*diff(tmp), mean(ylim), 'Optimization Matrix', 'Rotation', 90, ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
    end

    subplot(nRows, nColumns, 2*nColumns+iGroup);
    set(gca, 'Position', get(gca, 'Position'));
    plot(cell2mat(thVector(iGroup, :)'), cell2mat(zVector(iGroup, :)'), '.k', 'MarkerSize', 4);
    hold on;
    imagesc(binCentres{2}, binCentres{1}, fGroupMap{iGroup}, 'AlphaData', 0.8);
    axis xy tight off
    set(gca, 'CLim', [minF maxF]);
    box off;
    if iGroup == 1
        tmp = xlim;
        text(tmp(1)-0.2*diff(tmp), mean(ylim), 'Optimal Map', 'Rotation', 90, ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
    end
    
    subplot(nRows, nColumns, 3*nColumns+iGroup);
    set(gca, 'Position', get(gca, 'Position'));
    imagesc(traces{iGroup}.t, 1:size(traces{iGroup}.all, 1), traces{iGroup}.all);
    axis tight off
    set(gca, 'CLim', [minRawF maxRawF]);
    box off;
    if iGroup == 1
        tmp = xlim;
        text(tmp(1)-0.2*diff(tmp), mean(ylim), 'Trials', 'Rotation', 90, ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
    end

end

subplot(nRows, nColumns, [nColumns-1, nColumns, 2*nColumns-1, 2*nColumns]);
obj.ShowROI(iPlane, iROI);

subplot(nRows, nColumns, [3*nColumns-1, 3*nColumns, 4*nColumns-1, 4*nColumns]);
obj.dataTMaze.showPC(gca);

subplot(nRows, nColumns, ((nRows-1)*nColumns + 1):nRows*nColumns)
% set(gca, 'Position', get(gca, 'Position'));
plot(tData, fData);
xlabel('time [sec]');
ylabel('F');
axis tight
box off
% axis off;

set(figHandle, 'Name', sprintf('Cell # %d drawing...', iROI));
drawnow;

set(figHandle, 'Name', sprintf('Cell # %d done.', iROI));
end % PlotMaps()

%==========================================================================

function optOut = fillOptions(obj, optIn)

optOut = optIn;
if ~isfield(optIn, 'contrastGroups')
    contrasts = unique(obj.dataTMaze.contrastSequence);
    optOut.contrastGroups = {contrasts(contrasts<0);...
        0;...
        contrasts(contrasts>0);
        contrasts};
    optOut.contrastGroupLabels = {'Stim L';
        'Stim 0%';
        'Stim R';
        'All'};
end

if ~isfield(optIn, 'dZ')
    optOut.dZ = 3; % [cm]
end

if ~isfield(optIn, 'dTheta')
    optOut.dTheta = 2; % [deg]
end

if ~isfield(optIn, 'cvFactor')
    optOut.cvFactor = 5; %
end

if ~isfield(optIn, 'cvMethod')
    optOut.cvMethod = 'random'; % 'random' or 'block'
end

end %fillOptions();

%==========================================================================

function [stdOut, errVal, errValMatrix, gridValues] = estOptimalStd(coords, signal, binEdges)

coords = coords(:);
signal = signal(:);
nChunks = length(coords);
nDims = size(coords{1}, 2);

for iChunk = 1:nChunks
    trChunks = setdiff(1:nChunks, iChunk);
    trCoords{iChunk} = cell2mat(coords(trChunks));
    trSignal{iChunk} = cell2mat(signal(trChunks));
    testCoords{iChunk} = coords{iChunk};
    testSignal{iChunk} = signal{iChunk};
end

% precalculate the occupancy and signal maps (unsmoothed - filter will be
% our optimization parameter),
% also linear indices of the testCoordinates
for iChunk = 1:nChunks
    [occMaps{iChunk}, binCentres] = buildOccupMap(trCoords{iChunk}, binEdges);
    fMaps{iChunk} = buildAccumMap(trCoords{iChunk}, trSignal{iChunk}, binEdges);
    ind = nan(size(testCoords{iChunk}));
    [~, ind(:, 1)] = histc(testCoords{iChunk}(:,1), binEdges{1});
    [~, ind(:, 2)] = histc(testCoords{iChunk}(:,2), binEdges{2});
    testXYLlinIdx{iChunk} = sub2ind(size(occMaps{iChunk}), ind(:,1), ind(:,2));
end

% let's run a grid search for the optimal parameters

sz = size(occMaps{1});
for iDim = 1:nDims
    minValue = 1;
    maxValue = sz(iDim);
    pp = 3;
    gridValues{iDim} = linspace(minValue^(1/pp), maxValue^(1/pp), 15).^pp;
end

if nDims == 2
    for mInd = 1:length(gridValues{1})
        for nInd = 1:length(gridValues{2})
            x = [gridValues{1}(mInd), gridValues{2}(nInd)];
            errValMatrix(mInd, nInd) = mapError(x, occMaps, fMaps, testXYLlinIdx, testSignal);
        end
    end
else
    fprint('For now, grid search is implemented only for 2D grids\n');
end

ind = find(errValMatrix == min(errValMatrix(:)));
[mOpt, nOpt] = ind2sub(size(errValMatrix), ind);
stdOut = [gridValues{1}(mOpt), gridValues{2}(nOpt)];
errVal = errValMatrix(ind);

% % now run the optimization
% initStd = 5*ones(1, nDims);
% [stdOut, errVal] = fminsearch(@(x) mapError(x, occMaps, fMaps, testXYLlinIdx, testSignal), initStd);

end % estOptimalStd()

%==========================================================================

function errValue = mapError(x, occMaps, fMaps, testXYIdx, testF)

x(x<0.2) = 0.2; % HACK! This value is in pixels, not real units

hGauss = ndGaussian(x);

nChunks = length(occMaps);
chunkErrors = zeros(nChunks, 1);
for iChunk = 1:nChunks
    predictionMap = filterAndDivideMaps(occMaps{iChunk}, fMaps{iChunk}, hGauss);
    
    predictedF = predictionMap(testXYIdx{iChunk});
    
    %     chunkErrors(iChunk) = var(predictedF - testF{iChunk})/var(testF{iChunk});
    meanFTraining = sum(fMaps{iChunk}(:))/sum(occMaps{iChunk}(:));
    chunkErrors(iChunk) = sum((predictedF - testF{iChunk}).^2)/sum((testF{iChunk}-meanFTraining).^2);
end
errValue = median(chunkErrors);

end % mapError()
