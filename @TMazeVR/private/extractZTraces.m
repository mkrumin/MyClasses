function [dataTraces, modelTraces, extras] = extractZTraces(obj, res, options)

extras.report = obj.dataTMaze.report;
extras.outcome = obj.dataTMaze.outcome;
extras.contrasts = obj.dataTMaze.contrastSequence;
extras.mapData = struct('map', [], 'coords', {}, 'CoM', struct);
nTrials = obj.dataTMaze.nTrials;

thVector = cell(nTrials, 1);
zVector = cell(nTrials, 1);
fVector = cell(nTrials, max(obj.Planes));
zStd = cell(1, max(obj.Planes));
nSamples = zeros(nTrials, 1);

getCoords = true; % this is a flag of the first pass (the first plane is analyzed)
for iPlane = obj.Planes
    cellIdx = ismember(cell2mat(obj.data2p{iPlane}.ROI.CellClasses), options.cellClasses2Use);
    fData = obj.data2p{iPlane}.F(:, cellIdx);
    F0 = prctile(fData, 20);
    fData = bsxfun(@rdivide, bsxfun(@minus, fData, F0), F0);
    tData = obj.times2p{iPlane}';

    if getCoords
        for iTrial = 1:nTrials
            [thVector{iTrial}, zVector{iTrial}, fVector{iTrial, iPlane}, ~] =...
                buildVectors(obj, iTrial, tData, fData);
            nSamples(iTrial) = length(zVector{iTrial});
        end
        getCoords = false;
    else
        [~, ~, fVectorTmp, ~] = buildVectors(obj, 1:nTrials, tData, fData);
        fVector(:, iPlane) = mat2cell(fVectorTmp, nSamples, size(fVectorTmp, 2));
    end
    
    fModelTmp{iPlane} = nan(sum(nSamples), sum(cellIdx));
    cellNumbers = find(cellIdx);
    allTh = cell2mat(thVector(:));
    allZ = cell2mat(zVector(:));
    for iCell = 1:length(cellNumbers)
        map = res{iPlane}(cellNumbers(iCell)).zThetaMap;
        zCoords = res{iPlane}(cellNumbers(iCell)).zThetaBinCentres{1};
        thCoords = res{iPlane}(cellNumbers(iCell)).zThetaBinCentres{2};
        allTh = max(min(allTh, max(thCoords)), min(thCoords));
        fModelTmp{iPlane}(:, iCell) = interp2(thCoords, zCoords, map, allTh, allZ, 'linear');
        % getting the smoothing parameter (translated into [cm])
        zStd{iPlane}(1, iCell) = res{iPlane}(cellNumbers(iCell)).optStd(1)*...
            res{iPlane}(cellNumbers(iCell)).options.dZ;
        mapData.map = map;
        mapData.coords = res{iPlane}(cellNumbers(iCell)).zThetaBinCentres;
        thrMap = map - prctile(map(:), 80);
        thrMap = max(thrMap, 0);
        mapData.CoM = centerOfMass(thrMap, mapData.coords);
        extras.mapData(end+1) = mapData;
    end
end

fVector = cell2mat(fVector);
nCells = size(fVector, 2);
fVector = mat2cell(fVector, nSamples, nCells);
fVectorModel = cell2mat(fModelTmp(:)');
fVectorModel = mat2cell(fVectorModel, nSamples, nCells);
zStd = cell2mat(zStd(:)');

extras.zEdges = linspace(min(cell2mat(zVector)), max(cell2mat(zVector)), options.nBins+1);
dataTraces = binData(zVector, fVector, extras.zEdges);
modelTraces = binData(zVector, fVectorModel, extras.zEdges);

if options.scaleByContrast
    modelTraces = scaleMaps(obj, dataTraces, modelTraces);
end

%%  smooth the data in z with the same filter as the model (the
% cross-validated Gaussian filter)
% tic
dataTracesDenan = nan(size(dataTraces));
nanPositions = isnan(dataTraces);
zAxis = 1:options.nBins;
for iTrial = 1:nTrials
    legal = ~nanPositions(:, 1, iTrial);
    if sum(legal)>1
        dataTracesDenan(:, :, iTrial) = ...
            interp1(zAxis(legal), dataTraces(legal, :, iTrial), zAxis, 'previous', 'extrap');%, dataTraces(find(legal, 1, 'last'), iCell, iTrial));
    elseif sum(legal) == 1
        dataTracesDenan(:, :, iTrial) = ...
            repmat(dataTraces(legal, :, iTrial), options.nBins, 1);
    end
end
for iCell = 1:nCells
    hGauss = ndGaussian(zStd(iCell));
    dataTraces(:,iCell, :) = imfilter(dataTracesDenan(:, iCell, :), hGauss(:), 'replicate', 'same');
end
dataTraces(nanPositions) = NaN;
% toc

% ===================================================
function dataOut = binData(x, dataIn, binEdges)

binEdges(end) = Inf;
nBins = length(binEdges)-1;
nCells = size(dataIn{1}, 2);
nTrials = length(x);
dataOut = nan(nBins, nCells, nTrials);

for iTrial = 1:nTrials
    [~, ~, idx] = histcounts(x{iTrial}, binEdges);
    for iBin = 1:nBins
        tmp = nanmean(dataIn{iTrial}(idx == iBin, :), 1);
        if ~isempty(tmp)
            dataOut(iBin, :, iTrial) = tmp;
        end
    end
end
