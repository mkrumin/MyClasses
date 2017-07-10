function plotSnakes(data, results)

scaleByContrast = true;
separateLRColormaps = false;
climPrctiles = [50 99];
margin = 0.0;
cellClasses2Use = 'sad';
nBins = 50;

%% organizing the data
nDatasets = length(data);
for iDataset = 1:nDatasets
    fprintf('iDataset %d/%d\n', iDataset, nDatasets);
    obj = data(iDataset);
    res = results(iDataset).res;
    
    % take only finished trials
    idxAll = find(ismember(obj.dataTMaze.report, 'LR'));
    nTrials = length(idxAll);
    
    
    % get traces with uniform z
    
    thVector = cell(nTrials, 1);
    zVector = cell(nTrials, 1);
    zStd = cell(1, max(obj.Planes));
    nSamples = zeros(nTrials, 1);
    % tVector = cell(0);
    zMin = nan(nTrials, 1);
    zMax = nan(nTrials, 1);
    
    
    % get the trajectories information
    
    for iTrial = 1:nTrials
        [thVector{iTrial}, zVector{iTrial}, ~, ~] =...
            buildVectors(obj, idxAll(iTrial));%, tData, fData);
        if ~isempty(zVector{iTrial})
            zMin(iTrial) = min(zVector{iTrial});
            zMax(iTrial) = max(zVector{iTrial});
        end
        nSamples(iTrial) = length(zVector{iTrial});
    end

    % exclude empty trials
    validTrials = find(nSamples>0);
    nTrials = length(validTrials);
    nSamples = nSamples(validTrials);
    idxAll = idxAll(validTrials);
    zMin = zMin(validTrials);
    zMax = zMax(validTrials);
    thVector = thVector(validTrials);
    zVector = zVector(validTrials);
    
    zMinLimit = max(zMin);
    zMaxLimit = min(zMax);
    
    dZ = (zMaxLimit-zMinLimit)/(nBins+1);
    zEdges = buildEdges([zMinLimit, zMaxLimit], dZ);
    zAxis = (zEdges(1:end-1)+zEdges(2:end))/2;
    % the next line will collapse all the data from the final part of the
    % corridor (where the 'junction' is) into a single bin
    zEdges(end) = Inf;

    fVector = cell(nTrials, max(obj.Planes));
    tVector = cell(nTrials, max(obj.Planes));
    preF = cell(nTrials, max(obj.Planes));
    postF = cell(nTrials, max(obj.Planes));
    preT = cell(nTrials, max(obj.Planes));
    postT = cell(nTrials, max(obj.Planes));

    for iPlane = obj.Planes
        fprintf('iPlane %d\n', iPlane);
        cellIdx = ismember(cell2mat(obj.data2p{iPlane}.ROI.CellClasses), cellClasses2Use);
        fData = obj.data2p{iPlane}.F(:, cellIdx);
        F0 = prctile(fData, 20);
        fData = bsxfun(@minus, fData, F0);
        fData = bsxfun(@rdivide, fData, F0);
        
        tData = obj.times2p{iPlane}';
        
        [~, ~, fVectorTmp, tVectorTmp] =...
            buildVectors(obj, idxAll, tData, fData);
        [preF(:, iPlane), postF(:, iPlane), preT(:, iPlane), postT(:, iPlane)] = ...
            buildVectorsPrePost(obj, idxAll, tData, fData);
        fVector(:, iPlane) = mat2cell(fVectorTmp, nSamples, size(fVectorTmp, 2));
        tVector(:, iPlane) = mat2cell(tVectorTmp, nSamples, 1);
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
        
        
        % get the smoothing parameter used during model cross-validation
        for iCell = 1:length(cellNumbers)
            % the value in [pixels], translated into [cm]
            zStd{iPlane}(1, iCell) = res{iPlane}(cellNumbers(iCell)).optStd(1)*...
                res{iPlane}(cellNumbers(iCell)).options.dZ;
            % and then translated into bins (current resolution);
            zStd{iPlane}(1, iCell) = zStd{iPlane}(1, iCell)/dZ;
        end
        
    end
    
    fModelTmp = fModelTmp(:)'; % make sure it is a row array of cells
    fModel = cell(length(nSamples), size(fModelTmp, 2));
    for iPlane = obj.Planes
        fModel(:, iPlane) = mat2cell(fModelTmp{1, iPlane}, nSamples, size(fModelTmp{1, iPlane}, 2));
    end
    
    for iTrial = 1:nTrials
        fVector{iTrial} = cell2mat(fVector(iTrial, :));
        fModel{iTrial} = cell2mat(fModel(iTrial, :));
    end
    fVector = fVector(:,1);
    tVector = tVector(:,1);
    fModel = fModel(:,1);
    zStd = cell2mat(zStd);
    
    nCells = size(fVector{1}, 2);
    
    fTraces = nan(nBins, nCells, nTrials);
    fTracesT = nan(nBins, nCells, nTrials);
    fTracesModel = nan(nBins, nCells, nTrials);
    fTracesModelT = nan(nBins, nCells, nTrials);
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
        if nSamples(iTrial)>0
            fTracesT(:, :, iTrial) = resample(fVector{iTrial}, nBins, nSamples(iTrial));
            fTracesModelT(:, :, iTrial) = resample(fModel{iTrial}, nBins, nSamples(iTrial));
        end
    end
    
    % remove NaNs from fTraces and fTracesModel
    for iTrial = 1:nTrials
        nanIdx = find(isnan(fTraces(:, 1, iTrial)));
        if isempty(nanIdx)
            continue;
        end
        notNanIdx = find(~isnan(fTraces(:, 1, iTrial)));
        % use 'pchip' method to allow extrapolation
        fTraces(nanIdx, :, iTrial) = interp1(notNanIdx, ...
            fTraces(notNanIdx, :, iTrial), nanIdx, 'pchip', 'extrap');
        fTracesModel(nanIdx, :, iTrial) = interp1(notNanIdx, ...
            fTracesModel(notNanIdx, :, iTrial), nanIdx, 'pchip', 'extrap');
    end
    
    % smooth the data in z with the same filter as the model (the
    % cross-validated Gaussian filter)
    for iCell = 1:nCells
        hGauss = ndGaussian(zStd(iCell));
        fTraces(:,iCell, :) = imfilter(fTraces(:, iCell, :), hGauss(:), 'replicate', 'same');
    end
    
%     if scaleByContrast
%         cSequence = obj.dataTMaze.contrastSequence(idxAll);
%         cValues = unique(cSequence);
%         [~, cIndices] = ismember(cSequence, cValues);
%         fTracesModelScaled = scaleMaps(fTraces, fTracesModel, cIndices);
%     end

    % divide trials into groups (e.g. going L/R)
    % be careful - these are indices within idxAll, already excluding the unfinished
    % trials
    idxL = obj.dataTMaze.report(idxAll) == 'L';
    idxR = obj.dataTMaze.report(idxAll) == 'R';
    idxC = obj.dataTMaze.outcome(idxAll) == 'C';
    idxW = obj.dataTMaze.outcome(idxAll) == 'W';
    % only correct trials
    idxCL = idxL & idxC;
    idxCR = idxR & idxC;

    meanAll{iDataset} = nanmean(fTraces, 3);
    meanL{iDataset} = nanmean(fTraces(:,:,idxL), 3);
    meanR{iDataset} = nanmean(fTraces(:,:,idxR), 3);

    stdAll{iDataset} = std(fTraces, [], 3, 'omitnan');
    stdL{iDataset} = std(fTraces(:,:,idxL), [], 3, 'omitnan');
    stdR{iDataset} = std(fTraces(:,:,idxR), [], 3, 'omitnan');

    meanAllModel{iDataset} = nanmean(fTracesModel, 3);
    meanLModel{iDataset} = nanmean(fTracesModel(:,:,idxL), 3);
    meanRModel{iDataset} = nanmean(fTracesModel(:,:,idxR), 3);

    meanAllT{iDataset} = nanmean(fTracesT, 3);
    meanLT{iDataset} = nanmean(fTracesT(:,:,idxL), 3);
    meanRT{iDataset} = nanmean(fTracesT(:,:,idxR), 3);
    
    meanAllModelT{iDataset} = nanmean(fTracesModelT, 3);
    meanLModelT{iDataset} = nanmean(fTracesModelT(:,:,idxL), 3);
    meanRModelT{iDataset} = nanmean(fTracesModelT(:,:,idxR), 3);

    
%     meanAll{iDataset} = nanmedian(fTraces, 3);
%     meanL{iDataset} = nanmedian(fTraces(:,:,idxL), 3);
%     meanR{iDataset} = nanmedian(fTraces(:,:,idxR), 3);
%     
%     meanAllModel{iDataset} = nanmedian(fTracesModel, 3);
%     meanLModel{iDataset} = nanmedian(fTracesModel(:,:,idxL), 3);
%     meanRModel{iDataset} = nanmedian(fTracesModel(:,:,idxR), 3);
% 
%     meanAllT{iDataset} = nanmedian(fTracesT, 3);
%     meanLT{iDataset} = nanmedian(fTracesT(:,:,idxL), 3);
%     meanRT{iDataset} = nanmedian(fTracesT(:,:,idxR), 3);
%     
%     meanAllModelT{iDataset} = nanmedian(fTracesModelT, 3);
%     meanLModelT{iDataset} = nanmedian(fTracesModelT(:,:,idxL), 3);
%     meanRModelT{iDataset} = nanmedian(fTracesModelT(:,:,idxR), 3);

end % iDataset

% average across trials to find mean response
% use this to normalize the traces and to get latency;

%% figure out what cells are responsive to the task
[nTr, nPl] = size(preF);
meanPreF = cell(nTr, nPl);
meanPostF = cell(nTr, nPl);
for iPl = 1:nPl
    for iTr = 1:nTr
        meanPreF{iTr, iPl} = mean(preF{iTr, iPl});
        meanPostF{iTr, iPl} = mean(postF{iTr, iPl});
    end
end
meanPreF = cell2mat(meanPreF);
meanPostF = cell2mat(meanPostF);

[hResponsive] = ttest(meanPreF, meanPostF, 'Tail', 'left', 'alpha', 0.05);
respIdx = find(hResponsive);

%% preparing for plotting
meanAll = cell2mat(meanAll);
meanL = cell2mat(meanL);
meanR = cell2mat(meanR);
meanAllT = cell2mat(meanAllT);
meanLT = cell2mat(meanLT);
meanRT = cell2mat(meanRT);

meanAllModel = cell2mat(meanAllModel);
meanLModel = cell2mat(meanLModel);
meanRModel = cell2mat(meanRModel);
meanAllModelT = cell2mat(meanAllModelT);
meanLModelT = cell2mat(meanLModelT);
meanRModelT = cell2mat(meanRModelT);

meanAll = meanAll(:, respIdx);
meanL = meanL(:, respIdx);
meanR = meanR(:, respIdx);
meanAllT = meanAllT(:, respIdx);
meanLT = meanLT(:, respIdx);
meanRT = meanRT(:, respIdx);

meanAllModel = meanAllModel(:, respIdx);
meanLModel = meanLModel(:, respIdx);
meanRModel = meanRModel(:, respIdx);
meanAllModelT = meanAllModelT(:, respIdx);
meanLModelT = meanLModelT(:, respIdx);
meanRModelT = meanRModelT(:, respIdx);

nCells = size(meanAll, 2);

% h = mean(diff(zAxis));
% p = 1/(1+h^3/0.6);
% [fTracesMean, p] = csaps(zAxis, fTracesMean', p, zAxis);
% % [sp, fTracesMean, p] = spaps(zAxis, fTracesMean', 10.97);
% fTracesMean = fTracesMean';
% [fTracesMeanModel, p] = csaps(zAxis, fTracesMeanModel', p, zAxis);
% fTracesMeanModel = fTracesMeanModel';

% minValues = min(meanAll);
% meanAll = bsxfun(@minus, meanAll, minValues);
% maxValues = max(meanAll);
% meanAll = bsxfun(@rdivide, meanAll, maxValues);
% 
% meanAllModel = bsxfun(@minus, meanAllModel, minValues);
% meanAllModel = bsxfun(@rdivide, meanAllModel, maxValues);
% 
% minValues = min(meanAllT);
% meanAllT = bsxfun(@minus, meanAllT, minValues);
% maxValues = max(meanAllT);
% meanAllT = bsxfun(@rdivide, meanAllT, maxValues);
% 
% meanAllModelT = bsxfun(@minus, meanAllModelT, minValues);
% meanAllModelT = bsxfun(@rdivide, meanAllModelT, maxValues);


[~, latency] = max(meanAll);
[orderedLatency, order] = sort(latency);
% idxResponsive = ismember(order, find(hResponsive));
% order = order(idxResponsive);
% orderedLatency = orderedLatency(idxResponsive);
orderProper = order(orderedLatency>1);
[~, orderInv] = sort(order);

[~, latencyT] = max(meanAllT);
[orderedLatencyT, orderT] = sort(latencyT);
% idxResponsive = ismember(orderT, find(hResponsive));
% orderT = orderT(idxResponsive);
% orderedLatencyT = orderedLatencyT(idxResponsive);
orderProperT = orderT(orderedLatencyT>2);
[~, orderInvT] = sort(orderT);


%% plotting

figure
subplot(2, 2, 1)
imagesc(zAxis, 1:nCells, meanAll(:, order)');
imagesc(zAxis, 1:nCells, meanAll(:, orderProper)');
caxis(prctile(meanAll(:), climPrctiles))
axis xy
xlabel('z [cm]')
ylabel('All Cells');
title('Data (All trials)');
colorbar

subplot(2, 2, 2)
imagesc(zAxis, 1:nCells, meanAllModel(:, order)');
imagesc(zAxis, 1:nCells, meanAllModel(:, orderProper)');
caxis(prctile(meanAll(:), climPrctiles))
axis xy
xlabel('z [cm]')
% ylabel('Cell #');
title('SF Model (All trials)');
colorbar

subplot(2, 2, 3)
imagesc((1:nBins)/nBins, 1:nCells, meanAllT(:, orderT)');
imagesc((1:nBins)/nBins, 1:nCells, meanAllT(:, orderProperT)');
caxis(prctile(meanAllT(:), climPrctiles))
axis xy
xlabel('Normalized time [a.u.]')
ylabel('All Cells');
title('Data (All trials)');
colorbar

subplot(2, 2, 4)
imagesc((1:nBins)/nBins, 1:nCells, meanAllModelT(:, orderT)');
imagesc((1:nBins)/nBins, 1:nCells, meanAllModelT(:, orderProperT)');
caxis(prctile(meanAllT(:), climPrctiles))
axis xy
xlabel('Normalized Time [a.u.]')
% ylabel('Cell #');
title('SF Model (All trials)');
colorbar

%%
tmpRes = [res{:}];
for iCell = 1 :nCells
% iCell = 57;
set(gcf, 'Name', sprintf('iCell = %d', iCell))
nRows = 3;
nColumns = 3;
traces = squeeze(fTraces(:, iCell, :));
minF = min(traces(:));
maxF = max(traces(:));
% figure
subplot(nRows, nColumns, 1)
plot(1:nBins, traces(:, idxL))
hold on;
plot(1:nBins, mean(traces(:, idxL), 2), 'm', 'LineWidth', 3)
plot(1:nBins, median(traces(:, idxL), 2), 'm:', 'LineWidth', 3)
hold off
ylim([minF, maxF]);
xlim([1, nBins]);
title('Left Trials');
ylabel('Data');

subplot(nRows, nColumns, 2);
plot(1:nBins, traces(:, idxR));
hold on;
plot(1:nBins, mean(traces(:, idxR), 2), 'm', 'LineWidth', 3)
plot(1:nBins, median(traces(:, idxR), 2), 'm:', 'LineWidth', 3)
hold off
ylim([minF, maxF]);
xlim([1, nBins]);
title('Right Trials');

subplot(nRows, nColumns, 3);
cla;
meanL = mean(traces(:, idxL), 2);
meanR = mean(traces(:, idxR), 2);
semL = std(traces(:, idxL), [], 2)/sqrt(sum(idxL));
semR = std(traces(:, idxR), [], 2)/sqrt(sum(idxR));
alpha = 0.01;
[hRight, pVal, ci] = ttest2(traces(:, idxL), traces(:, idxR), ...
    'alpha', alpha, 'Dim', 2, 'Tail', 'left', 'Vartype', 'unequal');
[hLeft, pVal, ci] = ttest2(traces(:, idxL), traces(:, idxR), ...
    'alpha', alpha, 'Dim', 2, 'Tail', 'right', 'Vartype', 'unequal');

plot(1:nBins, meanL, 'r', 'LineWidth', 1)
hold on;
plot(1:nBins, meanR, 'b', 'LineWidth', 1)
xx = [1:nBins, nBins:-1:1]';
yy = [meanL+semL; flipud(meanL-semL)];
p = patch(xx, yy, 'r', 'FaceAlpha', 0.2, 'LineStyle', 'none');
xx = [1:nBins, nBins:-1:1]';
yy = [meanR+semR; flipud(meanR-semR)];
p = patch(xx, yy, 'b', 'FaceAlpha', 0.2, 'LineStyle', 'none');
scatter(find(hLeft), repmat(minF, size(find(hLeft))), [], 'r');
scatter(find(hRight), repmat(minF, size(find(hRight))), [], 'b');
hold off

ylim([minF, maxF]);
xlim([1, nBins]);
legend('Left trials', 'Right trials')


tracesM = squeeze(fTracesModel(:, iCell, :));

subplot(nRows, nColumns, 4)
plot(1:nBins, tracesM(:, idxL))
hold on;
plot(1:nBins, mean(tracesM(:, idxL), 2), 'c', 'LineWidth', 3)
plot(1:nBins, median(tracesM(:, idxL), 2), 'c:', 'LineWidth', 3)
hold off;
ylim([minF, maxF]);
xlim([1, nBins]);
ylabel('Model');

subplot(nRows, nColumns, 5);
plot(1:nBins, tracesM(:, idxR));
hold on;
plot(1:nBins, mean(tracesM(:, idxR), 2), 'c', 'LineWidth', 3)
plot(1:nBins, median(tracesM(:, idxR), 2), 'c:', 'LineWidth', 3)
hold off
ylim([minF, maxF]);
xlim([1, nBins]);

subplot(nRows, nColumns, 6)
cla
meanL = mean(tracesM(:, idxL), 2);
meanR = mean(tracesM(:, idxR), 2);
semL = std(tracesM(:, idxL), [], 2)/sqrt(sum(idxL));
semR = std(tracesM(:, idxR), [], 2)/sqrt(sum(idxR));

plot(1:nBins, meanL, 'r', 'LineWidth', 1)
hold on;
plot(1:nBins, meanR, 'b', 'LineWidth', 1)
xx = [1:nBins, nBins:-1:1]';
yy = [meanL+semL; flipud(meanL-semL)];
p = patch(xx, yy, 'r', 'FaceAlpha', 0.2, 'LineStyle', 'none');
xx = [1:nBins, nBins:-1:1]';
yy = [meanR+semR; flipud(meanR-semR)];
p = patch(xx, yy, 'b', 'FaceAlpha', 0.2, 'LineStyle', 'none');
hold off
ylim([minF, maxF]);
xlim([1, nBins]);
legend('Left trials', 'Right trials')

subplot(nRows, nColumns, 7)
plot(1:nBins, mean(traces(:, idxL), 2), 'm', 'LineWidth', 3)
hold on;
plot(1:nBins, median(traces(:, idxL), 2), 'm:', 'LineWidth', 3)
plot(1:nBins, mean(tracesM(:, idxL), 2), 'c', 'LineWidth', 3)
plot(1:nBins, median(tracesM(:, idxL), 2), 'c:', 'LineWidth', 3)
hold off;
ylim([minF, maxF]);
xlim([1, nBins]);
ylabel('Comparison');
legend('mean Data', 'median', 'mean Model', 'median');

subplot(nRows, nColumns, 8)
plot(1:nBins, mean(traces(:, idxR), 2), 'm', 'LineWidth', 3)
hold on;
plot(1:nBins, median(traces(:, idxR), 2), 'm:', 'LineWidth', 3)
plot(1:nBins, mean(tracesM(:, idxR), 2), 'c', 'LineWidth', 3)
plot(1:nBins, median(tracesM(:, idxR), 2), 'c:', 'LineWidth', 3)
hold off;
ylim([minF, maxF]);
xlim([1, nBins]);

subplot(nRows, nColumns, 9)
th = tmpRes(iCell).zThetaBinCentres{2};
zz = tmpRes(iCell).zThetaBinCentres{1};
imagesc(th, zz, tmpRes(iCell).zThetaMap)
axis xy equal tight

pause
end


%%


% divide cells into two groups, based on the maximum F
prefR = find(max(meanR)>(max(meanL)+margin));
prefL = find(max(meanL)>(max(meanR)+margin));
% OR, divide cells into two groups based on the mean F
% prefR = find(mean(meanR)-mean(meanL)>margin);
% prefL = find(mean(meanR)-mean(meanL)<-margin);

% minValues(prefR) = nanmin(meanR(:, prefR));
% minValues(prefL) = nanmin(meanL(:, prefL));
% meanL = bsxfun(@minus, meanL, minValues);
% meanR = bsxfun(@minus, meanR, minValues);
% maxValues(prefR) = nanmax(meanR(:, prefR));
% maxValues(prefL) = nanmax(meanL(:, prefL));
% meanL = bsxfun(@rdivide, meanL, maxValues);
% meanR = bsxfun(@rdivide, meanR, maxValues);
% 
% % minValues(prefR) = nanmin(meanRModel(:, prefR));
% % minValues(prefL) = nanmin(meanLModel(:, prefL));
% meanLModel = bsxfun(@minus, meanLModel, minValues);
% meanRModel = bsxfun(@minus, meanRModel, minValues);
% % maxValues(prefR) = nanmax(meanRModel(:, prefR));
% % maxValues(prefL) = nanmax(meanLModel(:, prefL));
% meanLModel = bsxfun(@rdivide, meanLModel, maxValues);
% meanRModel = bsxfun(@rdivide, meanRModel, maxValues);


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

figure

subplot(2, 4, 1)
imagesc(zAxis, 1:nCellsL, meanL(:, prefL(orderPrefL))');
caxis(climsL);
title('Left trials');
ylabel('Cells Prefer Left');
axis xy
colorbar;

subplot(2, 4, 2)
imagesc(zAxis, 1:nCellsL, meanR(:, prefL(orderPrefL))');
caxis(climsL);
title('Right trials');
axis xy
colorbar

subplot(2, 4, 5)
imagesc(zAxis, 1:nCellsR, meanL(:, prefR(orderPrefR))');
caxis(climsR);
ylabel('Cells prefer Right')
xlabel('z [cm]');
axis xy
colorbar;

subplot(2, 4, 6)
imagesc(zAxis, 1:nCellsR, meanR(:, prefR(orderPrefR))');
caxis(climsR);
xlabel('z [cm]');
axis xy
colorbar

% plotting the SF model-predicted 'snakes'

subplot(2, 4, 3)
imagesc(zAxis, 1:nCellsL, meanLModel(:, prefL(orderPrefL))');
caxis(climsL);
title('Left trials');
ylabel('Cells Prefer Left');
axis xy
colorbar;

subplot(2, 4, 4)
imagesc(zAxis, 1:nCellsL, meanRModel(:, prefL(orderPrefL))');
caxis(climsL);
title('Right trials');
axis xy
colorbar

subplot(2, 4, 7)
imagesc(zAxis, 1:nCellsR, meanLModel(:, prefR(orderPrefR))');
caxis(climsR);
ylabel('Cells prefer Right')
xlabel('z [cm]');
axis xy
colorbar;

subplot(2, 4, 8)
imagesc(zAxis, 1:nCellsR, meanRModel(:, prefR(orderPrefR))');
caxis(climsR);
xlabel('z [cm]');
axis xy
colorbar

%%

figure
iCell = 1;
subplot(1, 2, 1)
plot(zAxis, squeeze(fTraces(:, iCell, idxL)));
title('Went Left Trials');
xlabel('z [cm]');
ylabel('\DeltaF/F_0')
subplot(1, 2, 2)
plot(zAxis, squeeze(fTraces(:, iCell, idxR)));
title('Went Right Trials')
xlabel('z [cm]');

% keyboard;
