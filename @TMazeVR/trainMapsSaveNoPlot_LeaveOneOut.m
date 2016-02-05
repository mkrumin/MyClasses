function obj = trainMapsSaveNoPlot_LeaveOneOut(obj, iPlane, iROI, options)

if nargin<4
    options = struct;
end
options = fillOptions(obj, options);

fData = obj.data2p{iPlane}.F(:, iROI);

F0 = prctile(fData, 10);
fData = (fData-F0)/F0;

tData = obj.times2p{iPlane};
tData = tData - options.delay;

trialIdx = 1:length(obj.dataTMaze.contrastSequence);
nTrials = length(trialIdx);
[thAll, zAll, ~] = buildVectors(obj, trialIdx, tData, fData);

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
[optStd, errVals, errValMatrix, stdGridValues] = ...
        estOptimalStd(coords, fVector, {zEdges, thEdges});

% obj.trainingData{iPlane}(iROI).delay = options.delay;
obj.trainingData{iPlane}(iROI).optStd = optStd;
obj.trainingData{iPlane}(iROI).errVals = errVals;
obj.trainingData{iPlane}(iROI).errValMatrix = errValMatrix;
obj.trainingData{iPlane}(iROI).stdGridValues = stdGridValues;
obj.trainingData{iPlane}(iROI).zEdges = zEdges;
obj.trainingData{iPlane}(iROI).thEdges = thEdges;

hFilter = ndGaussian(optStd);
[theMap, binCentres] = estMap(cell2mat(coords), cell2mat(fVector), {zEdges, thEdges}, hFilter);

imagesc(binCentres{2}, binCentres{1}, theMap);
% title(sprintf('iPlane %d, iROI %d, err = %4.3f, delay = %4.2f', iPlane, iROI, errVals, options.delay), 'FontWeight', 'Normal');
title(sprintf('err = %4.3f, delay = %4.2f', errVals, options.delay), 'FontWeight', 'Normal');
colorbar;
axis xy equal tight;

obj.trainingData{iPlane}(iROI).options = options;
obj.trainingData{iPlane}(iROI).map = theMap;
obj.trainingData{iPlane}(iROI).binCentres = binCentres;



end % TrainMaps()

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

if ~isfield(optIn, 'delay')
    optOut.delay = 0; %
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
for iTrial = 1:nTrials
    if isempty(coords{iTrial})
        continue;
    end

    [occMaps(:,:,iTrial), binCentres] = buildOccupMap(coords{iTrial}, binEdges);
    fMaps(:,:,iTrial) = buildAccumMap(coords{iTrial}, signal{iTrial}, binEdges);
    ind = nan(size(coords{iTrial}));
    [~, ind(:, 1)] = histc(coords{iTrial}(:,1), binEdges{1});
    [~, ind(:, 2)] = histc(coords{iTrial}(:,2), binEdges{2});
    testXYLlinIdx{iTrial} = sub2ind(size(occMaps(:,:,1)), ind(:,1), ind(:,2));
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
            errValMatrix(mInd, nInd) = mapError(x, occMapsXOR, fMapsXOR, testXYLlinIdx, signal);
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

nTrials = length(occMaps);
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
