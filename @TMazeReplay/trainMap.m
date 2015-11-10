function obj = trainMap(obj, iPlane, iROI, options)

if nargin<4
    options = struct;
end
options = fillOptions(obj, options);

fData = obj.data2p{iPlane}.F(:, iROI);
tData = obj.times2p{iPlane};

% Fs = 1/median(diff(tData));
% % notch filter at 2 Hz
% wo = 2/(Fs/2);  bw = wo/5;
% [b,a] = iirnotch(wo,bw);
% % fvtool(b,a);
% fData = filtfilt(b,a,fDataRaw);


F0 = prctile(fData, 20);
fData = (fData-F0)/F0;


tData = tData - options.delay;

nTrials = obj.nTrials;
[thAll, zAll, ~] = buildVectors(obj, 1:nTrials, tData, fData);

zMin = min(abs(zAll));
zMax = max(abs(zAll));
thMax = max(abs(thAll));
clear zAll thAll;

% now, a complicated way to define bin edges
thEdges = buildEdges(thMax, options.dTheta);
zEdges = buildEdges([zMin, zMax], options.dZ);

% fMap = cell(nContrastGroups, options.cvFactor);
thVector = cell(nTrials, 1);
zVector = cell(nTrials, 1);
fVector = cell(nTrials, 1);
% traces = cell(nGroups, 1);
for iTrial = 1:nTrials
    [thVector{iTrial}, zVector{iTrial}, fVector{iTrial}] = buildVectors(obj, iTrial, tData, fData);
end

coords = cell(nTrials, 1);
% estimating the optimal filter STDs
for iTrial = 1:nTrials
    coords{iTrial} = [zVector{iTrial}, thVector{iTrial}];
end
% [optStd, errVals, errValMatrix, stdGridValues] = ...
%         estOptimalStd(coords, fVector, {zEdges, thEdges});
[optStd, errVals, errValMatrix, stdGridValues] = ...
        estOptimalStdFaster(coords, fVector, {zEdges, thEdges});

% obj.trainingData{iPlane}(iROI).delay = options.delay;
obj.trainingData{iPlane}(iROI).optStd = optStd;
obj.trainingData{iPlane}(iROI).errVals = errVals;
obj.trainingData{iPlane}(iROI).errValMatrix = errValMatrix;
obj.trainingData{iPlane}(iROI).stdGridValues = stdGridValues;
obj.trainingData{iPlane}(iROI).zEdges = zEdges;
obj.trainingData{iPlane}(iROI).thEdges = thEdges;

[theMap, binCentres] = estMap(cell2mat(coords), cell2mat(fVector), {zEdges, thEdges}, optStd);

obj.trainingData{iPlane}(iROI).map = theMap;
obj.trainingData{iPlane}(iROI).binCentres = binCentres;
obj.trainingData{iPlane}(iROI).options = options;

% plotting the map, useful for debugging
subplot(1, 2, 1)
imagesc(binCentres{2}, binCentres{1}, estMap(cell2mat(coords), cell2mat(fVector), {zEdges, thEdges}, [1 1]));
title(sprintf('iPlane %d, iROI %d', iPlane, iROI), 'FontWeight', 'Normal');
colorbar;
clim = caxis;
axis xy equal tight;

subplot(1, 2, 2)
imagesc(binCentres{2}, binCentres{1}, theMap);
colorbar;
% title(sprintf('iPlane %d, iROI %d, err = %4.3f, delay = %4.2f', iPlane, iROI, errVals, options.delay), 'FontWeight', 'Normal');
title(sprintf('err = %4.3f, delay = %4.2f', errVals, options.delay), 'FontWeight', 'Normal');
axis xy equal tight;
caxis(clim);

end % TrainMap()

%==========================================================================

function optOut = fillOptions(obj, optIn)

optOut = optIn;

if ~isfield(optIn, 'delay')
    % delay>0 means we will be looking at the neural activity relative to
    % positions in the past: tNew = tActual - delay;
    optOut.delay = 0; % [sec] 
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
nTrials = length(coords);
nDims = size(coords{1}, 2);

% precalculate the occupancy and signal maps (unsmoothed - filter will be
% our optimization parameter),
% also linear indices of the testCoordinates
occMaps = zeros(length(binEdges{1})-1, length(binEdges{2})-1, nTrials);
fMaps = zeros(length(binEdges{1})-1, length(binEdges{2})-1, nTrials);
for iTrial = 1:nTrials
    if isempty(coords{iTrial})
        continue;
    end

    [occMaps(:,:,iTrial), binCentres] = buildOccupMap(coords{iTrial}, binEdges);
    fMaps(:,:,iTrial) = buildAccumMap(coords{iTrial}, signal{iTrial}, binEdges);
    ind = nan(size(coords{iTrial}));
    [~, ind(:, 1)] = histc(coords{iTrial}(:,1), binEdges{1});
    [~, ind(:, 2)] = histc(coords{iTrial}(:,2), binEdges{2});
    testXYLinearIdx{iTrial} = sub2ind(size(occMaps(:,:,1)), ind(:,1), ind(:,2));
end

occMapsXOR = zeros(size(occMaps));
fMapsXOR = zeros(size(fMaps));
for iTrial = 1:nTrials
    idx = setdiff(1:nTrials, iTrial);
    occMapsXOR(:,:,iTrial) = sum(occMaps(:,:,idx), 3);
    fMapsXOR(:,:,iTrial) = sum(fMaps(:,:,idx), 3);
end

% let's run a grid search for the optimal parameters

sz = size(occMaps);
gridValues = cell(nDims, 1);
for iDim = 1:nDims
    minValue = 1;
    maxValue = sz(iDim);
    pp = 3;
    gridValues{iDim} = linspace(minValue^(1/pp), maxValue^(1/pp), 15).^pp;
end

if nDims == 2
    errValMatrix = nan(length(gridValues{1}), length(gridValues{2}));
    for mInd = 1:length(gridValues{1})
        for nInd = 1:length(gridValues{2})
            x = [gridValues{1}(mInd), gridValues{2}(nInd)];
            errValMatrix(mInd, nInd) = mapError(x, occMapsXOR, fMapsXOR, testXYLinearIdx, signal);
        end
    end
else
    fprint('For now, grid search is implemented only for 2D grids\n');
end

[errVal, ind] = min(errValMatrix(:));
[mOpt, nOpt] = ind2sub(size(errValMatrix), ind);
stdOut = [gridValues{1}(mOpt), gridValues{2}(nOpt)];


% % now run the optimization
% initStd = 5*ones(1, nDims);
% [stdOut, errVal] = fminsearch(@(x) mapError(x, occMaps, fMaps, testXYLlinIdx, testSignal), initStd);

end % estOptimalStd()

%==========================================================================

function errValue = mapError(x, occMaps, fMaps, testXYIdx, testF)

x(x<0.2) = 0.2; % HACK! This value is in pixels, not real units

hGauss{1} = ndGaussian(x(1));
hGauss{2} = ndGaussian(x(2))';

nTrials = size(occMaps, 3);
chunkErrors = nan(nTrials, 1);
for iTrial = 1:nTrials
    if isempty(testF{iTrial})
        % there is no data for this trial
        continue;
    end
    occM = occMaps(:,:,iTrial);
    fM = fMaps(:,:,iTrial);
    predictionMap = filterAndDivideMaps(occM, fM, hGauss);
    
    predictedF = predictionMap(testXYIdx{iTrial});
    
    %     chunkErrors(iChunk) = var(predictedF - testF{iChunk})/var(testF{iChunk});
    meanFTraining = sum(fM(:))/sum(occM(:));
    % the idea here is to check how much better the map is relative to just
    % using mean F map.
    chunkErrors(iTrial) = sum((predictedF - testF{iTrial}).^2)/sum((meanFTraining - testF{iTrial}).^2);
end
% errValue = nanmedian(chunkErrors);
errValue = nanmean(chunkErrors); % mean() works better, maps are smoother and nicer, "outliers have less effect"

end % mapError()

%=========================================================================
function [stdOut, errVal, errValMatrix, gridValues] = estOptimalStdFaster(coords, signal, binEdges)

coords = coords(:);
signal = signal(:);
nTrials = length(coords);
nDims = size(coords{1}, 2);

% precalculate the occupancy and signal maps (unsmoothed - filter will be
% our optimization parameter),
% also linear indices of the testCoordinates
nRows = length(binEdges{1})-1;
nColumns = length(binEdges{2})-1;
occMaps = zeros(nRows, nColumns, nTrials);
fMaps = zeros(nRows, nColumns, nTrials);
for iTrial = 1:nTrials
    if isempty(coords{iTrial})
        continue;
    end

    [occMaps(:,:,iTrial), binCentres] = buildOccupMap(coords{iTrial}, binEdges);
    fMaps(:,:,iTrial) = buildAccumMap(coords{iTrial}, signal{iTrial}, binEdges);
    ind = nan(size(coords{iTrial}));
    [~, ind(:, 1)] = histc(coords{iTrial}(:,1), binEdges{1});
    [~, ind(:, 2)] = histc(coords{iTrial}(:,2), binEdges{2});
    testXYLinearIdx{iTrial} = sub2ind(size(occMaps(:,:,1)), ind(:,1), ind(:,2));
end

occMapsXOR = zeros(size(occMaps));
fMapsXOR = zeros(size(fMaps));
for iTrial = 1:nTrials
    idx = setdiff(1:nTrials, iTrial);
    occMapsXOR(:,:,iTrial) = sum(occMaps(:,:,idx), 3);
    fMapsXOR(:,:,iTrial) = sum(fMaps(:,:,idx), 3);
end

% let's run a grid search for the optimal parameters

sz = size(occMaps);
gridValues = cell(nDims, 1);
for iDim = 1:nDims
    minValue = 1;
    maxValue = sz(iDim);
    pp = 3;
    gridValues{iDim} = linspace(minValue^(1/pp), maxValue^(1/pp), 15).^pp;
end

if nDims == 2
    errValMatrix = nan(length(gridValues{1}), length(gridValues{2}));
    for nInd = 1:length(gridValues{2})
        % prefilter along the 2nd dimension
        x2 = gridValues{2}(nInd);
        hGauss2 = ndGaussian(x2);
        
        meanVal = squeeze(sum(sum(fMapsXOR), 2)./sum(sum(occMapsXOR), 2));
        occMapsXORSpecial = permute(occMapsXOR, [2, 1, 3]);
        occMapsXORSpecial = reshape(occMapsXORSpecial, nColumns, nRows*nTrials);
        occMapsXORFilt = conv2(hGauss2, 1, occMapsXORSpecial, 'same');
        occMapsXORFilt = reshape(occMapsXORFilt, nColumns, nRows, nTrials);
        occMapsXORFilt = permute(occMapsXORFilt, [2, 1, 3]);
        fMapsXORSpecial = permute(fMapsXOR, [2, 1, 3]);
        fMapsXORSpecial = reshape(fMapsXORSpecial, nColumns, nRows*nTrials);
        fMapsXORFilt = conv2(hGauss2, 1, fMapsXORSpecial, 'same');
        fMapsXORFilt = reshape(fMapsXORFilt, nColumns, nRows, nTrials);
        fMapsXORFilt = permute(fMapsXORFilt, [2, 1, 3]);

        for mInd = 1:length(gridValues{1})
            x1 = gridValues{1}(mInd);
            errValMatrix(mInd, nInd) = mapError1D(x1, occMapsXORFilt, fMapsXORFilt, meanVal, testXYLinearIdx, signal);
        end
    end
else
    fprint('For now, grid search is implemented only for 2D grids\n');
end

[errVal, ind] = min(errValMatrix(:));
[mOpt, nOpt] = ind2sub(size(errValMatrix), ind);
stdOut = [gridValues{1}(mOpt), gridValues{2}(nOpt)];

% % now run the optimization
% initStd = 5*ones(1, nDims);
% [stdOut, errVal] = fminsearch(@(x) mapError(x, occMaps, fMaps, testXYLlinIdx, testSignal), initStd);

end % estOptimalStdFaster()

%==========================================================================

function errValue = mapError1D(x1, occMaps, fMaps, meanF, testXYIdx, testF)

x1(x1<0.2) = 0.2; % HACK! This value is in pixels, not real units

hGauss = ndGaussian(x1);

[nRows, nColumns, nTrials] = size(occMaps);
chunkErrors = nan(nTrials, 1);

epsilon = 0.01;
meanFMatrix = repmat(permute(meanF, [2, 3, 1]), nRows, nColumns, 1);
% allPredictionMaps = (imfilter(fMaps, hGauss)+epsilon*meanFMatrix)./(imfilter(occMaps, hGauss)+epsilon);
fMapsSpecial = reshape(fMaps, nRows, []);
occMapsSpecial = reshape(occMaps, nRows, []);
meanFSpecial = reshape(meanFMatrix, nRows, []);
allPredictionMaps = (conv2(hGauss, 1, fMapsSpecial, 'same')+epsilon*meanFSpecial)./(conv2(hGauss, 1, occMapsSpecial, 'same')+epsilon);
allPredictionMaps = reshape(allPredictionMaps, nRows, nColumns, nTrials);

for iTrial = 1:nTrials
    if isempty(testF{iTrial})
        % there is no data for this trial
        continue;
    end
%     occM = occMaps(:,:,iTrial);
%     fM = fMaps(:,:,iTrial);
%     predictionMap = filterAndDivideMaps1D(occM, fM, meanF(iTrial), hGauss);
    predictionMap = allPredictionMaps(:,:,iTrial);
    
    predictedF = predictionMap(testXYIdx{iTrial});
    
    %     chunkErrors(iChunk) = var(predictedF - testF{iChunk})/var(testF{iChunk});
    meanFTraining = meanF(iTrial);
    % the idea here is to check how much better the map is relative to just
    % using mean F map.
    chunkErrors(iTrial) = sum((predictedF - testF{iTrial}).^2)/sum((meanFTraining - testF{iTrial}).^2);
end
% errValue = nanmedian(chunkErrors);
errValue = nanmean(chunkErrors); % mean() works better, maps are smoother and nicer, "outliers have less effect"

end % mapErrorFaster()
