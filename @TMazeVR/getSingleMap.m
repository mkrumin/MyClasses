function [mapOut, binsOut] = getSingleMap(obj, iPlane, iROI, options)

% vSfN of TrainMaps - based on v6, data and modeled traces overlayed, only
% scaled model presented in the traces

if nargin<4
    options = struct;
end
options = fillOptions(obj, options);

nChunks = options.cvFactor;

fData = obj.data2p{iPlane}.F(:, iROI);
tData = obj.times2p{iPlane};

trialIdx = 1:length(obj.dataTMaze.contrastSequence);
validTrials = find(~isnan(obj.timesTmazeOnset));
trialIdx = intersect(trialIdx, validTrials);
[thAll, zAll, fAll] = buildVectors(obj, trialIdx, tData, fData);

zMax = max(abs(zAll));
zMin = min(abs(zAll));
thMax = max(abs(thAll));
% clear zAll thAll;

% now, a complicated way to define bin edges
thEdges = buildEdges(thMax, options.dTheta);
zEdges = buildEdges([zMin, zMax], options.dZ);

% first, estimate the optimal filter parameters

nTrials = length(trialIdx);
nTrialsPerSet = ceil(nTrials/nChunks);
randIdx = nan(nChunks, nTrialsPerSet);
randIdx(1:nTrials) = randperm(nTrials);
coords = cell(nChunks, 1);
fVector = cell(nChunks, 1);
zVector = cell(nChunks, 1);
thVector = cell(nChunks, 1);
for iChunk = 1:nChunks
    % extracting part of the data
    idx = randIdx(iChunk, :);
    idx = idx (~isnan(idx));
    [thVector{iChunk}, zVector{iChunk}, fVector{iChunk}] = buildVectors(obj, trialIdx(idx), tData, fData);
    coords{iChunk} = [zVector{iChunk}, thVector{iChunk}];
end

[optStd, errVals, ~, ~] = estOptimalStd(coords, fVector, {zEdges, thEdges});


hGauss = ndGaussian(optStd);
[theMap, binCentres] = estMap([zAll, thAll], fAll, {zEdges, thEdges}, hGauss);

mapOut = theMap;
binsOut = binCentres;

end
% return;
% 
% [thMatrix, zMatrix, fMatrix] = buildTraces(obj, 1:length(obj.dataTMaze.contrastSequence), tData, fData);

%% ==========================================================================

function optOut = fillOptions(obj, optIn)

optOut = optIn;
if ~isfield(optIn, 'contrastGroups')
    contrasts = unique(obj.dataTMaze.contrastSequence);
    optOut.contrastGroups = mat2cell(contrasts(:), ones(size(contrasts)), 1);
    optOut.contrastGroupLabels = cell(size(optOut.contrastGroups));
    for iLabel = 1:length(optOut.contrastGroups)
        optOut.contrastGroupLabels{iLabel} = [num2str(contrasts(iLabel)), '%'];
    end
    optOut.contrastGroups{end+1} = contrasts;
    optOut.contrastGroupLabels{end+1} = 'All';
end

if ~isfield(optIn, 'decisionGroups')
    optOut.decisionGroups = {'L'; 'R'; 'LR'};
    optOut.decisionGroupLabels = {'Go L'; 'Go R'; 'All'};
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

%=====================================================
function imOut = setNanColor(im, cmap, nanColor)
nColors = length(cmap);
nanColor = nanColor(:)';
cmapFull = [cmap; nanColor];
im = im*(nColors-1);
im(isnan(im)) = nColors;
idx = floor(im);
imOut = cmapFull(idx(:)+1, :);
imOut = reshape(imOut, [size(im), 3]);
end