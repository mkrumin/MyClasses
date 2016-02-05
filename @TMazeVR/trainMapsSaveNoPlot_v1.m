function obj = trainMapsSaveNoPlot_v1(obj, iPlane, iROI, options)

if nargin<4
    options = struct;
end
options = fillOptions(obj, options);

nGroups = length(options.contrastGroups);
nChunks = options.cvFactor;

fData = obj.data2p{iPlane}.F(:, iROI);
tData = obj.times2p{iPlane};

% figHandle = gcf;
% clf(figHandle, 'reset');
% set(figHandle, 'Name', sprintf('ROI # %d processing...', iROI));
% drawnow;
% set(figHandle, 'Units', 'normalized', 'MenuBar', 'none', 'OuterPosition', [0 0 1 1]);

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

obj.trainingData{iPlane}(iROI).optStd = optStd;
obj.trainingData{iPlane}(iROI).errVals = errVals;
obj.trainingData{iPlane}(iROI).errValMatrix = errValMatrix;
obj.trainingData{iPlane}(iROI).stdGridValues = stdGridValues;
obj.trainingData{iPlane}(iROI).zEdges = zEdges;
obj.trainingData{iPlane}(iROI).thEdges = thEdges;
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
