function obj = trainMotorMap_CVFR(obj, iPlane, iROI, options)

if nargin<4
    options = struct;
end
options = fillOptions(obj, options);

fData = obj.data2p{iPlane}.FR(:, iROI);

tData = obj.times2p{iPlane};
tData = tData - options.delay;

nTrials = obj.nTrials; % this is a method of TMazeReplay object
trialIdx = 1:nTrials;
[rAll, vAll, ~, ~] = buildMotorVectors(obj, trialIdx, tData, fData);

rMin = min(rAll);
rMax = max(rAll);
vMin = min(vAll);
vMax = max(vAll);
clear rAll vAll;

% now, a complicated way to define bin edges
rEdges = buildEdges([rMin, rMax], options.dR);
vEdges = buildEdges([vMin, vMax], options.dV);

nChunks = options.cvFactor;
scrambledIdx = randperm(nTrials);
rVector = cell(nChunks, 1);
vVector = cell(nChunks, 1);
fVector = cell(nChunks, 1);
endIdx = round([1:nChunks]/nChunks*nTrials);
startIdx = [1, endIdx(1:end-1)+1];
for iChunk = 1:nChunks
    trialIdx = scrambledIdx(startIdx(iChunk):endIdx(iChunk));
    [rVector{iChunk}, vVector{iChunk}, fVector{iChunk}] = buildMotorVectors(obj, trialIdx, tData, fData);
end

coords = cell(nChunks, 1);
% estimating the optimal filter STDs
for iChunk = 1:nChunks
    coords{iChunk} = [vVector{iChunk}, rVector{iChunk}];
end
[optStd, errVals, errValMatrix, stdGridValues] = ...
        estOptimalStd(coords, fVector, {vEdges, rEdges}, options.errPower);

obj.trainingData{iPlane}(iROI).optStd = optStd;
obj.trainingData{iPlane}(iROI).errVals = errVals;
obj.trainingData{iPlane}(iROI).errValMatrix = errValMatrix;
obj.trainingData{iPlane}(iROI).stdGridValues = stdGridValues;
obj.trainingData{iPlane}(iROI).vEdges = vEdges;
obj.trainingData{iPlane}(iROI).rEdges = rEdges;

[theMap, binCentres] = estMap(cell2mat(coords), cell2mat(fVector), {vEdges, rEdges}, optStd);
[rawMap, ~] = estMap(cell2mat(coords), cell2mat(fVector), {vEdges, rEdges}, [0.2 0.2]);

obj.trainingData{iPlane}(iROI).options = options;
obj.trainingData{iPlane}(iROI).rvMap = theMap;
obj.trainingData{iPlane}(iROI).rvBinCentres = binCentres;
obj.trainingData{iPlane}(iROI).rawDataMap = rawMap;

end % TrainMotorMaps_CV()

%==========================================================================
function optOut = fillOptions(obj, optIn)

optOut = optIn;

if ~isfield(optIn, 'dV') % V - velocity
    optOut.dV = 2; % [cm/sec]
end

if ~isfield(optIn, 'dR') % R - rotation
    optOut.dR = 2; % [deg/sec]
end

if ~isfield(optIn, 'cvFactor')
    optOut.cvFactor = 10; %
end

if ~isfield(optIn, 'delay')
    optOut.delay = 0; %
end

if ~isfield(optIn, 'errPower')
    optOut.errPower = 2; %
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
            errValMatrix(mInd, nInd) = mapError1D(x1, occMapsXORFilt, fMapsXORFilt, meanVal, testXYLinearIdx, signal, errPow);
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

end % estOptimalStdFaster()

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
