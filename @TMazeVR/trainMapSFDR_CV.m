function [zthErr, zthdErr, zrErr] = trainMapSFDR_CV(obj, iPlane, iROI, options)

if nargin<4
    options = struct;
end
options = fillOptions(obj, options);

fData = obj.data2p{iPlane}.F(:, iROI);

F0 = prctile(fData, 20);
fData = (fData-F0)/F0;

tData = obj.times2p{iPlane};
tData = tData - options.delay;

trialIdx = 1:length(obj.dataTMaze.contrastSequence);
nTrials = length(trialIdx);
thAll = cell(nTrials, 1);
zAll = cell(nTrials, 1);
rotAll = cell(nTrials, 1);
fAll = cell(nTrials, 1);
validTrials = true(nTrials, 1);
for iTrial = 1:nTrials
    [thAll{iTrial}, zAll{iTrial}, fAll{iTrial}, tAll{iTrial}] = buildVectors(obj, iTrial, tData, fData);
    [rotAll{iTrial}, ~, fRot{iTrial}, tRot{iTrial}] = buildMotorVectors(obj, iTrial, tData, fData);
    if ~isempty(tAll{iTrial})
        thAll{iTrial} = interp1(tAll{iTrial}, thAll{iTrial}, tRot{iTrial});
        zAll{iTrial} = interp1(tAll{iTrial}, zAll{iTrial}, tRot{iTrial});
        fAll{iTrial} = fRot{iTrial};
        tAll{iTrial} = tRot{iTrial};
    else
        rotAll{iTrial} = [];
    end
    validTrials(iTrial) = ~isempty(fAll{iTrial});
end

zMin = min(abs(cell2mat(zAll)));
zMax = max(abs(cell2mat(zAll)));
thMax = max(abs(cell2mat(thAll)));
rotMin = min(cell2mat(rotAll));
rotMax = max(cell2mat(rotAll));
% clear zAll thAll;

% now, a complicated way to define bin edges
thEdges = buildEdges(thMax, options.dTheta);
zEdges = buildEdges([zMin, zMax], options.dZ);
rotEdges = buildEdges([rotMin, rotMax], options.dTheta);

nChunks = options.cvFactor;
validLTrials = find(validTrials & ismember(obj.dataTMaze.report, 'L')');
validRTrials = find(validTrials & ismember(obj.dataTMaze.report, 'R')');
scrambledLIdx = validLTrials(randperm(length(validLTrials)));
scrambledRIdx = validRTrials(randperm(length(validRTrials)));
matrixL = nan(nChunks, ceil(length(scrambledLIdx)/nChunks));
matrixR = nan(nChunks, ceil(length(scrambledRIdx)/nChunks));
matrixL(1:numel(scrambledLIdx)) = scrambledLIdx(:);
matrixR(1:numel(scrambledRIdx)) = scrambledRIdx(:);
thVector = cell(nChunks, 1);
zVector = cell(nChunks, 1);
rotVector = cell(nChunks, 1);
fVector = cell(nChunks, 1);
zR = cell(nChunks, 1);
thR = cell(nChunks, 1);
rotR = cell(nChunks, 1);
fR = cell(nChunks, 1);
zL = cell(nChunks, 1);
thL = cell(nChunks, 1);
rotL = cell(nChunks, 1);
fL = cell(nChunks, 1);
for iChunk = 1:nChunks
    idxL = matrixL(iChunk, :)';
    idxL = idxL(~isnan(idxL));
    idxR = matrixR(iChunk, :)';
    idxR = idxR(~isnan(idxR));
    trialIdx = [idxL; idxR];
    thVector{iChunk} = cell2mat(thAll(trialIdx));
    zVector{iChunk} = cell2mat(zAll(trialIdx));
    rotVector{iChunk} = cell2mat(rotAll(trialIdx));
    fVector{iChunk} = cell2mat(fAll(trialIdx));
    zR{iChunk} = cell2mat(zAll(idxR));
    thR{iChunk} = cell2mat(thAll(idxR));
    rotR{iChunk} = cell2mat(rotAll(idxR));
    fR{iChunk} = cell2mat(fAll(idxR));
    zL{iChunk} = cell2mat(zAll(idxL));
    thL{iChunk} = cell2mat(thAll(idxL));
    rotL{iChunk} = cell2mat(rotAll(idxL));
    fL{iChunk} = cell2mat(fAll(idxL));
end

coords = cell(nChunks, 1);
coordsZR = cell(nChunks, 1);
coordsR = cell(nChunks, 1);
coordsL = cell(nChunks, 1);
% estimating the optimal filter STDs
for iChunk = 1:nChunks
    coords{iChunk} = [zVector{iChunk}, thVector{iChunk}];
    coordsZR{iChunk} = [zVector{iChunk}, rotVector{iChunk}];
    coordsR{iChunk} = [zR{iChunk}, thR{iChunk}];
    coordsL{iChunk} = [zL{iChunk}, thL{iChunk}];
end
[optStd, errVals, errValMatrix, stdGridValues] = ...
    estOptimalStd(coords, fVector, {zEdges, thEdges}, options.errPower);
[optStdZR, errValsZR, ~, ~] = ...
    estOptimalStd(coordsZR, fVector, {zEdges, rotEdges}, options.errPower);
[optStdR3D, errValR3D, ~, ~] = ...
    estOptimalStd(coordsR, fR, {zEdges, thEdges}, options.errPower);
[optStdL3D, errValL3D, ~, ~] = ...
    estOptimalStd(coordsL, fL, {zEdges, thEdges}, options.errPower);

% hFilter = ndGaussian(optStd);
% [zthMap, binCntrs] = estMap(cell2mat(coords), cell2mat(fVector), {zEdges, thEdges}, hFilter);
nR = length(validRTrials);
nL = length(validLTrials);
zthErr = errVals;
zthdErr = (errValL3D*nL + errValR3D*nR)/(nL + nR);
zrErr = errValsZR;
% z_d_Err = errVals15D;

end % TrainMaps()

%==========================================================================
function optOut = fillOptions(obj, optIn)

optOut = optIn;

if ~isfield(optIn, 'dZ')
    optOut.dZ = 3; % [cm]
end

if ~isfield(optIn, 'dTheta')
    optOut.dTheta = 2; % [deg]
end

if ~isfield(optIn, 'cvFactor')
    optOut.cvFactor = 10; %
end

if ~isfield(optIn, 'delay')
    optOut.delay = 0; %
end


end %fillOptions();

%==========================================================================
function [stdOut, errVal, errValMatrix, gridValues] = estOptimalStd(coords, signal, binEdges, errPow)

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
testXYLinearIdx = cell(nTrials, 1);
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
            errValMatrix(mInd, nInd) = mapError(x1, occMapsXORFilt, fMapsXORFilt, meanVal, testXYLinearIdx, signal, errPow);
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
function errValue = mapError(x1, occMaps, fMaps, meanF, testXYIdx, testF, errPow)

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
    switch errPow
        case 1
            chunkErrors(iTrial) = sum(abs(predictedF - testF{iTrial}))/sum(abs(meanFTraining - testF{iTrial}));
        case 2
            chunkErrors(iTrial) = sum((predictedF - testF{iTrial}).^2)/sum((meanFTraining - testF{iTrial}).^2);
        otherwise
            chunkErrors(iTrial) = sum(abs(predictedF - testF{iTrial}).^errPow)/sum(abs(meanFTraining - testF{iTrial}).^errPow);
    end
end
% errValue = nanmedian(chunkErrors);
errValue = nanmean(chunkErrors); % mean() works better, maps are smoother and nicer, "outliers have less effect"

end % mapErrorFaster()

%==========================================================================
function [stdOut, errVal, errValMatrix, gridValues] = estOptimalStd1D(coords, signal, binEdges, errPow)

coords = coords(:);
signal = signal(:);
nTrials = length(coords);
nDims = size(coords{1}, 2);

% precalculate the occupancy and signal maps (unsmoothed - filter will be
% our optimization parameter),
% also linear indices of the testCoordinates
nRows = length(binEdges{1})-1;
occMaps = zeros(nRows, 1, nTrials);
fMaps = zeros(nRows, 1, nTrials);
testXYLinearIdx = cell(nTrials, 1);
for iTrial = 1:nTrials
    if isempty(coords{iTrial})
        continue;
    end
    
    [occMaps(:,:,iTrial), binCentres] = buildOccupMap(coords{iTrial}, binEdges);
    fMaps(:,:,iTrial) = buildAccumMap(coords{iTrial}, signal{iTrial}, binEdges);
    %     ind = nan(size(coords{iTrial}));
    [~, ind] = histc(coords{iTrial}, binEdges{1});
    %     [~, ind(:, 2)] = histc(coords{iTrial}(:,2), binEdges{2});
    testXYLinearIdx{iTrial} = ind;
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

switch nDims
    case 1
        errValMatrix = nan(length(gridValues{1}), 1);
        
        meanVal = squeeze(sum(sum(fMapsXOR), 2)./sum(sum(occMapsXOR), 2));
        
        for mInd = 1:length(gridValues{1})
            x1 = gridValues{1}(mInd);
            errValMatrix(mInd) = mapError1D(x1, occMapsXOR, fMapsXOR, meanVal, testXYLinearIdx, signal, errPow);
        end
    otherwise
        fprintf('For now, grid search is implemented only for 2D grids\n');
end

ind = find(errValMatrix == min(errValMatrix(:)));
mOpt = ind;
stdOut = [gridValues{1}(mOpt)];
errVal = errValMatrix(ind);

% % now run the optimization
% initStd = 5*ones(1, nDims);
% [stdOut, errVal] = fminsearch(@(x) mapError(x, occMaps, fMaps, testXYLlinIdx, testSignal), initStd);

end % estOptimalStd()

%==========================================================================
function errValue = mapError1D(x1, occMaps, fMaps, meanF, testXYIdx, testF, errPow)

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
    switch errPow
        case 1
            chunkErrors(iTrial) = sum(abs(predictedF - testF{iTrial}))/sum(abs(meanFTraining - testF{iTrial}));
        case 2
            chunkErrors(iTrial) = sum((predictedF - testF{iTrial}).^2)/sum((meanFTraining - testF{iTrial}).^2);
        otherwise
            chunkErrors(iTrial) = sum(abs(predictedF - testF{iTrial}).^errPow)/sum(abs(meanFTraining - testF{iTrial}).^errPow);
    end
end
% errValue = nanmedian(chunkErrors);
errValue = nanmean(chunkErrors); % mean() works better, maps are smoother and nicer, "outliers have less effect"

end % mapErrorFaster()

%==========================================================================
function [stdOut, errVal, errValMatrix, gridValues] = estOptimalStd25D(coords2D, signal2D, binEdges, errPow)

% coords = coords(:);
% signal = signal(:);
[nTrials, nSubsets] = size(coords2D);
coords = cell(nTrials, 1);
signal = cell(nTrials, 1);
coordsXOR2D = cell(nTrials, nSubsets);
signalXOR2D = cell(nTrials, nSubsets);
for iTrial = 1:nTrials
    coords{iTrial, 1} = cat(1, coords2D{iTrial, 1}, coords2D{iTrial, 2});
    signal{iTrial, 1} = cat(1, signal2D{iTrial, 1}, signal2D{iTrial, 2});
    idx = setdiff(1:nTrials, iTrial);
    for iSubset = 1:nSubsets
        coordsXOR2D{iTrial, iSubset} = cell2mat(coords2D(idx, iSubset));
        signalXOR2D{iTrial, iSubset} = cell2mat(signal2D(idx, iSubset));
    end
end
nDims = size(coords{1}, 2);

% precalculate the occupancy and signal maps (unsmoothed - filter will be
% our optimization parameter),
% also linear indices of the testCoordinates
nRows = length(binEdges{1})-1;
nColumns = length(binEdges{2})-1;
occMaps = zeros(nRows, nColumns, nTrials);
fMaps = zeros(nRows, nColumns, nTrials);
testXYLinearIdx = cell(nTrials, nSubsets);
for iTrial = 1:nTrials
    if isempty(coords{iTrial})
        continue;
    end
    
    [occMaps(:,:,iTrial), binCentres] = buildOccupMap(coords{iTrial}, binEdges);
    fMaps(:,:,iTrial) = buildAccumMap(coords{iTrial}, signal{iTrial}, binEdges);
    for iSubset = 1:nSubsets
        ind = nan(size(coords2D{iTrial, iSubset}));
        [~, ind(:, 1)] = histc(coords2D{iTrial, iSubset}(:,1), binEdges{1});
        [~, ind(:, 2)] = histc(coords2D{iTrial, iSubset}(:,2), binEdges{2});
        testXYLinearIdx{iTrial, iSubset} = sub2ind([nRows, nColumns], ind(:,1), ind(:,2));
        indTrain = nan(size(coordsXOR2D{iTrial, iSubset}));
        [~, indTrain(:, 1)] = histc(coordsXOR2D{iTrial, iSubset}(:,1), binEdges{1});
        [~, indTrain(:, 2)] = histc(coordsXOR2D{iTrial, iSubset}(:,2), binEdges{2});
        trainXYLinearIdx{iTrial, iSubset} = sub2ind([nRows, nColumns], indTrain(:,1), indTrain(:,2));
    end
end
% these are training sets (XOR because they exclude trial #1 data)
occMapsXOR = zeros(size(occMaps));
fMapsXOR = zeros(size(fMaps));
for iTrial = 1:nTrials
    idx = setdiff(1:nTrials, iTrial);
    occMapsXOR(:,:,iTrial) = sum(occMaps(:,:,idx), 3);
    fMapsXOR(:,:,iTrial) = sum(fMaps(:,:,idx), 3);
    meanVal(iTrial, 1) = mean(cell2mat(signal2D(idx, 1)));
    meanVal(iTrial, 2) = mean(cell2mat(signal2D(idx, 2)));
    meanVal(iTrial, 3) = mean(cell2mat(signal(idx)));
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
            errValMatrix(mInd, nInd) = ...
                mapError25D(x1, occMapsXORFilt, fMapsXORFilt, meanVal, trainXYLinearIdx, signalXOR2D, testXYLinearIdx, signal2D, errPow);
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
function errValue = mapError25D(x1, occMaps, fMaps, meanF, trainXYIdx, trainF, testXYIdx, testF, errPow)

x1(x1<0.2) = 0.2; % HACK! This value is in pixels, not real units

hGauss = ndGaussian(x1);

[nRows, nColumns, nTrials] = size(occMaps);
chunkErrors = nan(nTrials, 1);

epsilon = 0.01;
meanFMatrix = repmat(permute(meanF(:,3), [2, 3, 1]), nRows, nColumns, 1);
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
    
    predictionMap = allPredictionMaps(:,:,iTrial);
    
    predictedF{1} = predictionMap(testXYIdx{iTrial, 1});
    predictedF{2} = predictionMap(testXYIdx{iTrial, 2});
    trainPredictedF{1} = predictionMap(trainXYIdx{iTrial, 1});
    trainPredictedF{2} = predictionMap(trainXYIdx{iTrial, 2});
    pp = getLinearScaling(trainPredictedF{1}, trainF{iTrial, 1});
    predictedF{1} = predictedF{1}*pp(1) + pp(2);
    pp = getLinearScaling(trainPredictedF{2}, trainF{iTrial, 2});
    predictedF{2} = predictedF{2}*pp(1) + pp(2);
    
    meanFTraining = meanF(iTrial);
    % the idea here is to check how much better the map is relative to just
    % using mean F map.
    switch errPow
        case 1
            chunkErrors(iTrial) = sum(abs(predictedF - testF{iTrial}))/sum(abs(meanFTraining - testF{iTrial}));
        case 2
            err1 = sum((predictedF{1}-testF{iTrial, 1}).^2);
            base1 = sum((testF{iTrial, 1}-meanF(iTrial, 1)).^2);
            err2 = sum((predictedF{2}-testF{iTrial, 2}).^2);
            base2 = sum((testF{iTrial, 2}-meanF(iTrial, 2)).^2);
            len1 = length(predictedF{1});
            len2 = length(predictedF{2});
            
            chunkErrors(iTrial) = (err1/base1*len1+err2/base2*len2)/(len1+len2);
            %             sum((predictedF - testF{iTrial}).^2)/sum((meanFTraining - testF{iTrial}).^2);
        otherwise
            chunkErrors(iTrial) = sum(abs(predictedF - testF{iTrial}).^errPow)/sum(abs(meanFTraining - testF{iTrial}).^errPow);
    end
end
% errValue = nanmedian(chunkErrors);
errValue = nanmean(chunkErrors); % mean() works better, maps are smoother and nicer, "outliers have less effect"

end % mapErrorFaster()

%============================================================================
function  out = getLinearScaling(xx, yy)

mY = mean(yy);
mX = mean(xx);
% sc = pinv(xx-mX)*(yy-mY);
sc = ((xx-mX)'*(yy-mY))/((xx-mX)'*(xx-mX));
out(1) = sc;
out(2) =  mY-sc*mX;
out = out(:);

end
