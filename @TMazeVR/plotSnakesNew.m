function plotSnakesNew(data, results)

%% setting some parameters
separateLRColormaps = false;
climPrctiles = [50 99];
margin = 0.1;
correctOnly = true;

options.nBins = 50;
options.cellClasses2Use = 's'; % 'sad' for all ROIs, 'sa' to exclude dendrites
options.scaleByContrast = false;


%% extract and re-organize data
nDatasets = length(data);
for iDataset = 1:nDatasets
    fprintf('iDataset %d/%d\n', iDataset, nDatasets);
    obj = data(iDataset);
    res = results(iDataset).res;
    
    [fTraces, fTracesModel, extras] = extractZTraces(obj, res, options);
    
    idxL = extras.report == 'L';
    idxR = extras.report == 'R';
    if correctOnly
        idxL = idxL & extras.outcome == 'C';
        idxR = idxR & extras.outcome == 'C';
    end
    meanAll{iDataset} = nanmean(fTraces, 3);
    meanL{iDataset} = nanmean(fTraces(:,:,idxL), 3);
    meanR{iDataset} = nanmean(fTraces(:,:,idxR), 3);
    
    meanAllModel{iDataset} = nanmean(fTracesModel, 3);
    meanLModel{iDataset} = nanmean(fTracesModel(:,:,idxL), 3);
    meanRModel{iDataset} = nanmean(fTracesModel(:,:,idxR), 3);
    
    cValues = unique(abs(extras.contrasts));
    nContrasts = length(cValues);
    for iContrast = 1:nContrasts
        idxC = (abs(extras.contrasts) == cValues(iContrast))';
        meanLC{iContrast, iDataset} = nanmean(fTraces(:, :, idxL & idxC), 3);
        meanRC{iContrast, iDataset} = nanmean(fTraces(:, :, idxR & idxC), 3);
        meanLCModel{iContrast, iDataset} = nanmean(fTracesModel(:, :, idxL & idxC), 3);
        meanRCModel{iContrast, iDataset} = nanmean(fTracesModel(:, :, idxR & idxC), 3);
    end
    
    tmp = [extras.mapData.CoM];
    sfCenters{iDataset} = [tmp.center];
    
end % iDataset

% average across trials to find mean response
% use this to normalize the traces and to get latency;

%% calculate mean responses, for L and R trials, and also separate by contrast
fTracesMean = cell2mat(meanAll);
meanL = cell2mat(meanL);
meanR = cell2mat(meanR);
nCells = size(fTracesMean, 2);
fTracesMeanModel = cell2mat(meanAllModel);
meanLModel = cell2mat(meanLModel);
meanRModel = cell2mat(meanRModel);

meanLC = cell2mat(meanLC);
meanRC = cell2mat(meanRC);
meanLCModel = cell2mat(meanLCModel);
meanRCModel = cell2mat(meanRCModel);

zAxis = (extras.zEdges(1:end-1)+extras.zEdges(2:end))/2;
sfCenters = cell2mat(sfCenters(:)');

minValues = nanmin(fTracesMean);
fTracesMean = bsxfun(@minus, fTracesMean, minValues);
maxValues = nanmax(fTracesMean);
fTracesMean = bsxfun(@rdivide, fTracesMean, maxValues);

fTracesMeanModel = bsxfun(@minus, fTracesMeanModel, minValues);
fTracesMeanModel = bsxfun(@rdivide, fTracesMeanModel, maxValues);

[~, latency] = max(fTracesMean);
[~, order] = sort(latency);

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

meanLModel = bsxfun(@minus, meanLModel, minValues);
meanRModel = bsxfun(@minus, meanRModel, minValues);
meanLModel = bsxfun(@rdivide, meanLModel, maxValues);
meanRModel = bsxfun(@rdivide, meanRModel, maxValues);

meanLC = bsxfun(@minus, meanLC, minValues);
meanRC = bsxfun(@minus, meanRC, minValues);
meanLC = bsxfun(@rdivide, meanLC, maxValues);
meanRC = bsxfun(@rdivide, meanRC, maxValues);
meanLCModel = bsxfun(@minus, meanLCModel, minValues);
meanRCModel = bsxfun(@minus, meanRCModel, minValues);
meanLCModel = bsxfun(@rdivide, meanLCModel, maxValues);
meanRCModel = bsxfun(@rdivide, meanRCModel, maxValues);

meanLC = mat2cell(meanLC, options.nBins*ones(nContrasts, 1), nCells);
meanRC = mat2cell(meanRC, options.nBins*ones(nContrasts, 1), nCells);
meanLCModel = mat2cell(meanLCModel, options.nBins*ones(nContrasts, 1), nCells);
meanRCModel = mat2cell(meanRCModel, options.nBins*ones(nContrasts, 1), nCells);

nCellsL = numel(prefL);
nCellsR = numel(prefR);

%% plot the snakes, while ordering the cells according to the peak latency
% in the overall mean response

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

%% SF centers plot, all contrasts on one plot
figure
plot(sfCenters(2, prefL), sfCenters(1, prefL), 'ko', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
hold on;
plot(sfCenters(2, prefR), sfCenters(1, prefR), 'ko', 'MarkerSize', 8, 'MarkerFaceColor', 'b');

xlabel('\theta [deg]')
ylabel('z [cm]')
title('Position of SF centers');
% histogram(sfCenters(2, prefL), 20)
% histogram(sfCenters(2, prefR), 20)
leg = legend('L_{pref} Cells', 'R_{pref} Cells');
set(leg, 'Location', 'bestoutside');

plot([0 0], ylim, 'k:');
axis equal
xlim([-30 30]);
ylim([5 105]);

%% plotting contrast-wise SF position maps
figure
respLGlobal = max(meanL);
respRGlobal = max(meanR);
for iContrast = 1:nContrasts
%     h(iContrast) =subplot(2, 3, iContrast);
    h(iContrast) = axes('OuterPosition', [(iContrast-1)/nContrasts, 0, 1/nContrasts, 1]);
    respL = max(meanLC{iContrast});
    respR = max(meanRC{iContrast});
    R = (respR+respL);
    min(R);
    R = R+6;
    score = (respR-respL)./(respR+respL);
%     score = (respR-respL)./(respRGlobal+respLGlobal);
    maxVal = 1.0;
    minVal = 0.2;
    score(abs(score)>maxVal) = sign(score(abs(score)>maxVal));
    score(abs(score)<=minVal) = 0;
    score = (abs(score)-minVal)/(maxVal-minVal).*sign(score);
    idx = find(score);
%     length(idx)
    red = 1-(score+abs(score))/2;
    green = 1-abs(score);
    blue = 1+(score-abs(score))/2;
%     red = respL./(respR+respL);
%     blue = respR./(respR+respL);
    colorVector = [red(:), green(:), blue(:)];
    scatter(sfCenters(2, idx), sfCenters(1, idx), max(10*R(idx), 1), colorVector(idx, :), 'filled', 'MarkerEdgeColor', 0.5*[1 1 1]);
%     scatter(sfCenters(2, idx), sfCenters(1, idx), 10*R(idx), colorVector(idx, :), 'filled');
%     scatter(sfCenters(2, idx), sfCenters(1, idx), 10*R(idx), colorVector(idx, :), 'filled');
    title(['\pm' sprintf('%d%%', cValues(iContrast))])
    xlim([-30 30]);
    ylim([5 105]);
    
    nBinsTh = 9;
    nBinsZ = 15;
    thEdges = linspace(-30, 30, nBinsTh+1);
    zEdges =  linspace(5, 105, nBinsZ+1);
    thAxis = (thEdges(1:end-1)+thEdges(2:end))/2;
    zAxis = (zEdges(1:end-1)+zEdges(2:end))/2;
    [~, ~, thBins] = histcounts(sfCenters(2, :), thEdges);
    [~, ~, zBins] = histcounts(sfCenters(1, :), zEdges);
    
    idxR = score>0;
    idxL = score<0;
    cumR  = sparse(zBins(idxR), thBins(idxR), score(idxR));
    cumL  = sparse(zBins(idxL), thBins(idxL), score(idxL));
    cumR = cumR/sum(cumR(:))*50;
    cumL = cumL/sum(cumL(:))*50;
%     figure
hold on;
%     plot(thAxis, sum(cumR, 1), 'c', thAxis, sum(-cumL, 1), 'm', 'LineWidth', 3)
%     plot(sum(cumR, 2)-30, zAxis, 'c:', sum(-cumL, 2)-30, zAxis, 'm:', 'LineWidth', 3)
    baseValue = min(ylim);
    b(1) = bar(thAxis, -[sum(cumL, 1)]+baseValue, 'BaseValue', baseValue, ...
        'FaceColor', [1 0.5 0.5], 'EdgeColor', 'r', 'LineWidth', 2);
    b(2) = bar(thAxis, -[sum(cumR, 1)]+baseValue, 'BaseValue', baseValue, ...
        'FaceColor', 'none', 'EdgeColor', 'b', 'LineWidth', 2);
    baseValue = -30;
    b(3) = barh(zAxis, -2*[sum(cumL, 2)]+baseValue, 'BaseValue', baseValue, ...
        'FaceColor', [1 0.5 0.5], 'EdgeColor', 'r', 'LineWidth', 2);
    b(4) = barh(zAxis, -2*[sum(cumR, 2)]+baseValue, 'BaseValue', baseValue, ...
        'FaceColor', 'none', 'EdgeColor', 'b', 'LineWidth', 2);
    set(b, 'ShowBaseline', 'off')
    axis tight off
    
    zCenterR = sum(sfCenters(1, idxR).*score(idxR))/sum(score(idxR));
    zSpreadR = sqrt(sum((sfCenters(1, idxR)-zCenterR).^2.*score(idxR))/sum(score(idxR)));
    thCenterR = sum(sfCenters(2, idxR).*score(idxR))/sum(score(idxR));
    thSpreadR = sqrt(sum((sfCenters(2, idxR)-thCenterR).^2.*score(idxR))/sum(score(idxR)));

    zCenterL = sum(sfCenters(1, idxL).*score(idxL))/sum(score(idxL));
    zSpreadL = sqrt(sum((sfCenters(1, idxL)-zCenterL).^2.*score(idxL))/sum(score(idxL)));
    thCenterL = sum(sfCenters(2, idxL).*score(idxL))/sum(score(idxL));
    thSpreadL = sqrt(sum((sfCenters(2, idxL)-thCenterL).^2.*score(idxL))/sum(score(idxL)));
    
    plot(thCenterL + [-1 1]*thSpreadL, [zCenterL, zCenterL], 'r', 'LineWidth', 3);
    plot([thCenterL, thCenterL], zCenterL + [-1 1]*zSpreadL , 'r', 'LineWidth', 3);
    plot(thCenterR + [-1 1]*thSpreadR, [zCenterR, zCenterR], 'b', 'LineWidth', 3);
    plot([thCenterR, thCenterR], zCenterR + [-1 1]*zSpreadR , 'b', 'LineWidth', 3);

end

% equalize the axes limits
for iH = 1:length(h)
    axes(h(iH));
    xx(iH,:) = xlim;
    yy(iH,:) = ylim;
end
xx = [min(xx(:, 1)), max(xx(:, 2))];
yy = [min(yy(:, 1)), max(yy(:, 2))];
for iH = 1:length(h)
    axes(h(iH));
    xlim(xx);
    ylim(yy);
end

%% plotting snakes, data separated by contrast

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
