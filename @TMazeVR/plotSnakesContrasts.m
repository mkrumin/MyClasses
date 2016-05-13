function plotSnakesContrasts(data, results)

separateLRColormaps = false;
climPrctiles = [50 99];
margin = 0.0;

nDatasets = length(data);
for iDataset = 1:nDatasets
    fprintf('iDataset %d/%d\n', iDataset, nDatasets);
    obj = data(iDataset);
    res = results(iDataset).res;
    % take only legal trials
    
    idxAll = find(ismember(obj.dataTMaze.report, 'LR'));
    nTrials = length(idxAll);
    % divide trials into groups (e.g. going L/R)
    
    idxL = obj.dataTMaze.report(idxAll) == 'L';
    idxR = obj.dataTMaze.report(idxAll) == 'R';
    idxC = obj.dataTMaze.outcome(idxAll) == 'C';
    idxW = obj.dataTMaze.outcome(idxAll) == 'W';

% taking only correct trials    
%     idxL = idxL & idxC;
%     idxR = idxR & idxC;
    
    % get traces with uniform z
    
    % [thVector, zVector, fVector, tVector] = buildVectors(obj, trialIdx, tData, fData)
    
    thVector = cell(nTrials, 1);
    zVector = cell(nTrials, 1);
    fVector = cell(nTrials, max(obj.Planes));
    zStd = cell(1, max(obj.Planes));
    nSamples = zeros(nTrials, 1);
    % tVector = cell(0);
    zMin = [];
    zMax = [];
    
    getCoords = true; % this is a flag of the first pass (the first plane is analyzed)
    for iPlane = obj.Planes
        fprintf('iPlane %d\n', iPlane);
        cellIdx = cell2mat(obj.data2p{iPlane}.ROI.CellClasses) == 's';
        %     cellIdx = ismember(cell2mat(obj.data2p{iPlane}.ROI.CellClasses), 'sa');
        fData = obj.data2p{iPlane}.F(:, cellIdx);
        F0 = prctile(fData, 20);
        fData = bsxfun(@minus, fData, F0);
        fData = bsxfun(@rdivide, fData, F0);
        
        tData = obj.times2p{iPlane}';
        if getCoords
            for iTrial = 1:nTrials
                [thVector{iTrial}, zVector{iTrial}, fVector{iTrial, iPlane}, ~] =...
                    buildVectors(obj, idxAll(iTrial), tData, fData);
                if ~isempty(zVector{iTrial})
                    zMin(iTrial) = min(zVector{iTrial});
                    zMax(iTrial) = max(zVector{iTrial});
                else
                    zMin(iTrial) = NaN;
                    zMax(iTrial) = NaN;
                end
                nSamples(iTrial) = length(zVector{iTrial});
            end
            getCoords = false;
            
            zMinLimit = max(zMin);
            zMaxLimit = min(zMax);
            
            nBins = 50;
            zEdges = linspace(zMinLimit, zMaxLimit, nBins+1);
            zAxis = (zEdges(1:end-1)+zEdges(2:end))/2;
            dZ = diff(zAxis(1:2));
            zEdges(end) = Inf;
        else
            [~, ~, fVectorTmp, ~] =...
                buildVectors(obj, idxAll, tData, fData);
            fVector(:, iPlane) = mat2cell(fVectorTmp, nSamples, size(fVectorTmp, 2));
        end
        %             res{iPlane}
        
        fModelTmp{iPlane} = nan(sum(nSamples), sum(cellIdx));
        cellNumbers = find(cellIdx);
        allTh = cell2mat(thVector(:));
        allZ = cell2mat(zVector(:));
        for iCell = 1:length(cellNumbers)
            map = res{iPlane}(cellNumbers(iCell)).zThetaMap;
            zCoords = res{iPlane}(cellNumbers(iCell)).zThetaBinCentres{1};
            thCoords = res{iPlane}(cellNumbers(iCell)).zThetaBinCentres{2};
            allTh(allTh>max(thCoords)) = max(thCoords);
            allTh(allTh<min(thCoords)) = min(thCoords);
            fModelTmp{iPlane}(:, iCell) = interp2(thCoords, zCoords, map, allTh, allZ, 'linear');
        end
        
 
        % get the smoothing parameter used in training
        for iCell = 1:length(cellNumbers)
            % the value in pixels is calculated in [cm]
            zStd{iPlane}(1, iCell) = res{iPlane}(cellNumbers(iCell)).optStd(1)*...
                res{iPlane}(cellNumbers(iCell)).options.dZ;
        end
        
    end
    
    fModelTmp = fModelTmp(:)'; % make sure it is a row array of cells
    fModel = cell(length(nSamples), size(fModelTmp, 2));
    for iPlane = obj.Planes
        fModel(:, iPlane) = mat2cell(fModelTmp{1, iPlane}, nSamples, size(fModelTmp{1, iPlane}, 2));
    end
    
    for iTrial = 1:nTrials;
        fVector{iTrial} = cell2mat(fVector(iTrial, :));
        fModel{iTrial} = cell2mat(fModel(iTrial, :));
    end
    fVector = fVector(:,1);
    fModel = fModel(:,1);
    zStd = cell2mat(zStd);
    
    nCells = size(fVector{1}, 2);
    
    fTraces = nan(nBins, nCells, nTrials);
    fTracesModel = nan(nBins, nCells, nTrials);
    for iTrial = 1:nTrials
        [~, ~, idx] = histcounts(zVector{iTrial}, zEdges);
        for iBin = 1:nBins
            tmp = mean(fVector{iTrial}(idx == iBin, :), 1);
            if ~isempty(tmp)
                fTraces(iBin, :, iTrial) = tmp;
            end
            tmp = mean(fModel{iTrial}(idx == iBin, :), 1);
            if ~isempty(tmp)
                fTracesModel(iBin, :, iTrial) = tmp;
            end
        end
    end
    
    % smooth the data in z with the same filter as the model (the
    % cross-validated Gaussian filter)
    for iCell = 1:nCells
        hGauss = ndGaussian(zStd(iCell));
        fTraces(:,iCell, :) = imfilter(fTraces(:, iCell, :), hGauss(:), 'replicate', 'same');
    end
    
    meanAll{iDataset} = nanmean(fTraces, 3);
    meanL{iDataset} = nanmean(fTraces(:,:,idxL), 3);
    meanR{iDataset} = nanmean(fTraces(:,:,idxR), 3);
    
    meanAllModel{iDataset} = nanmean(fTracesModel, 3);
    meanLModel{iDataset} = nanmean(fTracesModel(:,:,idxL), 3);
    meanRModel{iDataset} = nanmean(fTracesModel(:,:,idxR), 3);
    
    cSequence = obj.dataTMaze.contrastSequence(idxAll);
    cValues = unique(cSequence);
    [~, cIndices] = ismember(cSequence, cValues);
    nContrasts = length(cValues);
    idxLC = cell(nContrasts, 1);
    idxRC = cell(nContrasts, 1);
    for iContrast = 1:nContrasts
        idxLC{iContrast} = idxL' & (cIndices == iContrast);
        idxRC{iContrast} = idxR' & (cIndices == iContrast);
        meanLC{iContrast, iDataset} = nanmean(fTraces(:,:,idxLC{iContrast}), 3);
        meanRC{iContrast, iDataset} = nanmean(fTraces(:,:,idxRC{iContrast}), 3);
        meanLCModel{iContrast, iDataset} = nanmean(fTracesModel(:,:,idxLC{iContrast}), 3);
        meanRCModel{iContrast, iDataset} = nanmean(fTracesModel(:,:,idxRC{iContrast}), 3);
    end

end % iDataset

% average across trials to find mean response
% use this to normalize the traces and to get latency;

%% plotting
fTracesMean = cell2mat(meanAll);
meanL = cell2mat(meanL);
meanR = cell2mat(meanR);
nCells = size(fTracesMean, 2);
fTracesMeanModel = cell2mat(meanAllModel);
meanLModel = cell2mat(meanLModel);
meanRModel = cell2mat(meanRModel);
tmpBinning = nBins*ones(nContrasts, 1);
meanLC = cell2mat(meanLC);
meanRC = cell2mat(meanRC);
meanLCModel = cell2mat(meanLCModel);
meanRCModel = cell2mat(meanRCModel);
meanLC = mat2cell(meanLC, tmpBinning, nCells);
meanRC = mat2cell(meanRC, tmpBinning, nCells);
meanLCModel = mat2cell(meanLCModel, tmpBinning, nCells);
meanRCModel = mat2cell(meanRCModel, tmpBinning, nCells);

% h = mean(diff(zAxis));
% p = 1/(1+h^3/0.6);
% [fTracesMean, p] = csaps(zAxis, fTracesMean', p, zAxis);
% % [sp, fTracesMean, p] = spaps(zAxis, fTracesMean', 10.97);
% fTracesMean = fTracesMean';
% [fTracesMeanModel, p] = csaps(zAxis, fTracesMeanModel', p, zAxis);
% fTracesMeanModel = fTracesMeanModel';

minValues = nanmin(fTracesMean);
fTracesMean = bsxfun(@minus, fTracesMean, minValues);
maxValues = nanmax(fTracesMean);
fTracesMean = bsxfun(@rdivide, fTracesMean, maxValues);

% minValues = nanmin(fTracesMeanModel);
fTracesMeanModel = bsxfun(@minus, fTracesMeanModel, minValues);
% maxValues = nanmax(fTracesMeanModel);
fTracesMeanModel = bsxfun(@rdivide, fTracesMeanModel, maxValues);

[~, latency] = max(fTracesMean);
[~, order] = sort(latency);

figure
subplot(4, 4, [1 2 5 6])
% imagesc(zAxis, 1:nCells, fTracesMean(:, order)');
imagesc(zAxis, 1:nCells, fTracesMean(:, order)');
caxis(prctile(fTracesMean(:), climPrctiles))

axis xy
% xlabel('z [cm]')
ylabel('All Cells');
title('Data (All trials)');
colorbar

subplot(4, 4, [3 4 7 8])
% plot(zAxis, fTracesMean(:, order))
imagesc(zAxis, 1:nCells, fTracesMeanModel(:, order)');
caxis(prctile(fTracesMean(:), climPrctiles))

axis xy
% xlabel('z [cm]')
% ylabel('Cell #');
title('SF Model (All trials)');
colorbar


% meanL = csaps(zAxis, meanL', p, zAxis)';
% meanR = csaps(zAxis, meanR', p, zAxis)';

% meanLModel = csaps(zAxis, meanLModel', p, zAxis)';
% meanRModel = csaps(zAxis, meanRModel', p, zAxis)';

% divide cells into two groups, based on the maximum F
prefR = find(max(meanR)>(max(meanL)+margin));
prefL = find(max(meanL)>(max(meanR)+margin));
% OR, divide cells into two groups based on the mean F
% prefR = find(mean(meanR)-mean(meanL)>margin);
% prefL = find(mean(meanR)-mean(meanL)<-margin);

minValues(prefR) = nanmin(meanR(:, prefR));
minValues(prefL) = nanmin(meanL(:, prefL));
meanL = bsxfun(@minus, meanL, minValues);
meanR = bsxfun(@minus, meanR, minValues);
maxValues(prefR) = nanmax(meanR(:, prefR));
maxValues(prefL) = nanmax(meanL(:, prefL));
meanL = bsxfun(@rdivide, meanL, maxValues);
meanR = bsxfun(@rdivide, meanR, maxValues);

% minValues(prefR) = nanmin(meanRModel(:, prefR));
% minValues(prefL) = nanmin(meanLModel(:, prefL));
meanLModel = bsxfun(@minus, meanLModel, minValues);
meanRModel = bsxfun(@minus, meanRModel, minValues);
% maxValues(prefR) = nanmax(meanRModel(:, prefR));
% maxValues(prefL) = nanmax(meanLModel(:, prefL));
meanLModel = bsxfun(@rdivide, meanLModel, maxValues);
meanRModel = bsxfun(@rdivide, meanRModel, maxValues);


nCellsL = numel(prefL);
nCellsR = numel(prefR);

% plot the snakes, while ordering the cells according to the peak latency
% in the overall mean response

% clims = [min([meanL(:); meanR(:)]'), max([meanL(:); meanR(:)]')];
clims = prctile([meanL(:); meanR(:)]', climPrctiles);
if separateLRColormaps
    climsL = prctile(reshape([meanL(:, prefL), meanR(:, prefL)], 1, []), climPrctiles);
    climsR = prctile(reshape([meanL(:, prefR), meanR(:, prefR)], 1, []), climPrctiles);
else
    climsR = clims;
    climsL = clims;
end

% subplot(4, 4, 9)
% imagesc(zAxis, 1:nCellsL, meanL(:, order(ismember(order, prefL)))');
% caxis(climsL);
% title('Left trials');
% ylabel('Cells Prefer Left');
% axis xy
% colorbar;
% 
% subplot(4, 4, 10)
% imagesc(zAxis, 1:nCellsL, meanR(:, order(ismember(order, prefL)))');
% caxis(climsL);
% title('Right trials');
% axis xy
% colorbar
% 
% subplot(4, 4, 13)
% imagesc(zAxis, 1:nCellsR, meanL(:, order(ismember(order, prefR)))');
% caxis(climsR);
% ylabel('Cells prefer Right')
% axis xy
% colorbar;
% 
% subplot(4, 4, 14)
% imagesc(zAxis, 1:nCellsR, meanR(:, order(ismember(order, prefR)))');
% caxis(climsR);
% axis xy
% colorbar

% And now order the cells 'within group': left-preferring cells by the left
% trial response latencies, right preferring cells - by right trial
% response latencies

[~, latency] = max(meanL(:, prefL));
[~, orderPrefL] = sort(latency);

[~, latency] = max(meanR(:, prefR));
[~, orderPrefR] = sort(latency);

% control ordering - L cells by their R response latency, and vice-versa

% [~, latency] = max(meanR(:, prefL));
% [~, orderPrefL] = sort(latency);
%
% [~, latency] = max(meanL(:, prefR));
% [~, orderPrefR] = sort(latency);

subplot(4, 4, 9)
imagesc(zAxis, 1:nCellsL, meanL(:, prefL(orderPrefL))');
caxis(climsL);
title('Left trials');
ylabel('Cells Prefer Left');
axis xy
colorbar;

subplot(4, 4, 10)
imagesc(zAxis, 1:nCellsL, meanR(:, prefL(orderPrefL))');
caxis(climsL);
title('Right trials');
axis xy
colorbar

subplot(4, 4, 13)
imagesc(zAxis, 1:nCellsR, meanL(:, prefR(orderPrefR))');
caxis(climsR);
ylabel('Cells prefer Right')
xlabel('z [cm]');
axis xy
colorbar;

subplot(4, 4, 14)
imagesc(zAxis, 1:nCellsR, meanR(:, prefR(orderPrefR))');
caxis(climsR);
xlabel('z [cm]');
axis xy
colorbar

% plotting the SF model-predicted 'snakes'

subplot(4, 4, 11)
imagesc(zAxis, 1:nCellsL, meanLModel(:, prefL(orderPrefL))');
caxis(climsL);
title('Left trials');
ylabel('Cells Prefer Left');
axis xy
colorbar;

subplot(4, 4, 12)
imagesc(zAxis, 1:nCellsL, meanRModel(:, prefL(orderPrefL))');
caxis(climsL);
title('Right trials');
axis xy
colorbar

subplot(4, 4, 15)
imagesc(zAxis, 1:nCellsR, meanLModel(:, prefR(orderPrefR))');
caxis(climsR);
ylabel('Cells prefer Right')
xlabel('z [cm]');
axis xy
colorbar;

subplot(4, 4, 16)
imagesc(zAxis, 1:nCellsR, meanRModel(:, prefR(orderPrefR))');
caxis(climsR);
xlabel('z [cm]');
axis xy
colorbar

%%
for iContrast = 1:nContrasts
    meanLC{iContrast} = bsxfun(@minus, meanLC{iContrast}, minValues);
    meanLC{iContrast} = bsxfun(@rdivide, meanLC{iContrast}, maxValues);
    meanRC{iContrast} = bsxfun(@minus, meanRC{iContrast}, minValues);
    meanRC{iContrast} = bsxfun(@rdivide, meanRC{iContrast}, maxValues);
    meanLCModel{iContrast} = bsxfun(@minus, meanLCModel{iContrast}, minValues);
    meanLCModel{iContrast} = bsxfun(@rdivide, meanLCModel{iContrast}, maxValues);
    meanRCModel{iContrast} = bsxfun(@minus, meanRCModel{iContrast}, minValues);
    meanRCModel{iContrast} = bsxfun(@rdivide, meanRCModel{iContrast}, maxValues);
end

%%
doColorbar = false;
for iContrast = 1:nContrasts
    figure('Name', sprintf('%d %% Contrast', cValues(iContrast)));
    
    subplot(2, 4, 1)
    imagesc(zAxis, 1:nCellsL, meanLC{iContrast}(:, prefL(orderPrefL))');
    caxis(climsL);
    title('Left trials');
    ylabel('Cells Prefer Left');
    axis xy
    if doColorbar
    colorbar;
    end
    
    subplot(2, 4, 2)
    imagesc(zAxis, 1:nCellsL, meanRC{iContrast}(:, prefL(orderPrefL))');
    caxis(climsL);
    title('Right trials');
    axis xy
    if doColorbar
    colorbar;
    end
    
    subplot(2, 4, 5)
    imagesc(zAxis, 1:nCellsR, meanLC{iContrast}(:, prefR(orderPrefR))');
    caxis(climsR);
    ylabel('Cells prefer Right')
    xlabel('z [cm]');
    axis xy
    if doColorbar
    colorbar;
    end
    
    subplot(2, 4, 6)
    imagesc(zAxis, 1:nCellsR, meanRC{iContrast}(:, prefR(orderPrefR))');
    caxis(climsR);
    xlabel('z [cm]');
    axis xy
    if doColorbar
    colorbar;
    end

    subplot(2, 4, 3)
    imagesc(zAxis, 1:nCellsL, meanLCModel{iContrast}(:, prefL(orderPrefL))');
    caxis(climsL);
    title('Left trials');
    ylabel('Cells Prefer Left');
    axis xy
    if doColorbar
    colorbar;
    end
    
    subplot(2, 4, 4)
    imagesc(zAxis, 1:nCellsL, meanRCModel{iContrast}(:, prefL(orderPrefL))');
    caxis(climsL);
    title('Right trials');
    axis xy
    if doColorbar
    colorbar;
    end
    
    subplot(2, 4, 7)
    imagesc(zAxis, 1:nCellsR, meanLCModel{iContrast}(:, prefR(orderPrefR))');
    caxis(climsR);
    ylabel('Cells prefer Right')
    xlabel('z [cm]');
    axis xy
    if doColorbar
    colorbar;
    end
    
    subplot(2, 4, 8)
    imagesc(zAxis, 1:nCellsR, meanRCModel{iContrast}(:, prefR(orderPrefR))');
    caxis(climsR);
    xlabel('z [cm]');
    axis xy
    if doColorbar
    colorbar;
    end
end
%%
% keyboard;
