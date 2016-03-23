function plotSnakes(data)

separateColormap = false;
climPrctiles = [50 99];

nDatasets = length(data);
for iDataset = 1:nDatasets
    obj = data(iDataset);
% take only legal trials

idxAll = find(ismember(obj.dataTMaze.report, 'LR'));
nTrials = length(idxAll);
% divide trials into groups (e.g. going L/R)

idxL = obj.dataTMaze.report(idxAll) == 'L';
idxR = obj.dataTMaze.report(idxAll) == 'R';
idxC = obj.dataTMaze.outcome(idxAll) == 'C';
idxW = obj.dataTMaze.outcome(idxAll) == 'W';

% get traces with uniform z

% [thVector, zVector, fVector, tVector] = buildVectors(obj, trialIdx, tData, fData)

thVector = cell(0);
zVector = cell(0);
fVector = cell(0);
% tVector = cell(0);
zMin = [];
zMax = [];

getCoords = true; % this is a flag of the first pass (the first plane is analyzed)
for iPlane = obj.Planes
    cellIdx = cell2mat(obj.data2p{iPlane}.ROI.CellClasses) == 's';
%     cellIdx = ismember(cell2mat(obj.data2p{iPlane}.ROI.CellClasses), 'sa');
    fData = obj.data2p{iPlane}.F(:, cellIdx);
    tData = obj.times2p{iPlane}';
    for iTrial = 1:nTrials
        if getCoords
            [thVector{iTrial}, zVector{iTrial}, fVector{iTrial, iPlane}, ~] =...
                buildVectors(obj, idxAll(iTrial), tData, fData);
            if ~isempty(zVector{iTrial})
                zMin(iTrial) = min(zVector{iTrial});
                zMax(iTrial) = max(zVector{iTrial});
            else
                zMin(iTrial) = NaN;
                zMax(iTrial) = NaN;
            end
        else
            [~, ~, fVector{iTrial, iPlane}, ~] =...
                buildVectors(obj, idxAll(iTrial), tData, fData);
        end
    end
    getCoords = false;
end

for iTrial = 1:nTrials
    fVector{iTrial} = cell2mat(fVector(iTrial, :));
end
fVector = fVector(:,1);

zMinLimit = max(zMin);
zMaxLimit = min(zMax);

nBins = 50;
zEdges = linspace(zMinLimit, zMaxLimit, nBins+1);
zAxis = (zEdges(1:end-1)+zEdges(2:end))/2;
zEdges(end) = Inf;

nCells = size(fVector{1}, 2);

fTraces = nan(nBins, nCells, nTrials);
for iTrial = 1:nTrials
    [~, ~, idx] = histcounts(zVector{iTrial}, zEdges);
    for iBin = 1:nBins
        tmp = mean(fVector{iTrial}(idx == iBin, :), 1);
        if ~isempty(tmp)
        fTraces(iBin, :, iTrial) = tmp;
        end
    end
end

meanAll{iDataset} = nanmean(fTraces, 3);
meanL{iDataset} = nanmean(fTraces(:,:,idxL), 3);
meanR{iDataset} = nanmean(fTraces(:,:,idxR), 3);

end % iDataset

% average across trials to find mean response
% use this to normalize the traces and to get latency;

%% plotting
fTracesMean = cell2mat(meanAll);
meanL = cell2mat(meanL);
meanR = cell2mat(meanR);
nCells = size(fTracesMean, 2);

h = mean(diff(zAxis));
p = 1/(1+h^3/0.6);
[fTracesMean, p] = csaps(zAxis, fTracesMean', p, zAxis);
% [sp, fTracesMean, p] = spaps(zAxis, fTracesMean', 10.97);
fTracesMean = fTracesMean';

minValues = nanmin(fTracesMean);
fTracesMean = bsxfun(@minus, fTracesMean, minValues);
maxValues = nanmax(fTracesMean);
fTracesMean = bsxfun(@rdivide, fTracesMean, maxValues);

[~, latency] = max(fTracesMean);
[~, order] = sort(latency);
[~, orderInv] = sort(order);

figure
subplot(2, 2, 1)
% imagesc(zAxis, 1:nCells, fTracesMean(:, order)');
imagesc(zAxis, 1:nCells, fTracesMean(:, order)');
caxis(prctile(fTracesMean(:), climPrctiles))

axis xy
xlabel('z [cm]')
ylabel('Cell #');
colorbar

subplot(2, 2, 2)
plot(zAxis, fTracesMean(:, order))

meanL = csaps(zAxis, meanL', p, zAxis)';
meanR = csaps(zAxis, meanR', p, zAxis)';

% divide cells into two groups, based on the maximum F
prefR = find(max(meanR)>(max(meanL)+0.0));
prefL = find(max(meanL)>(max(meanR)+0.0));
% OR, divide cells into two groups based on the mean F
% prefR = find(mean(meanR)-mean(meanL)>0.1);
% prefL = find(mean(meanR)-mean(meanL)<-0.1);

minValues(prefR) = nanmin(meanR(:, prefR));
minValues(prefL) = nanmin(meanL(:, prefL));
meanL = bsxfun(@minus, meanL, minValues);
meanR = bsxfun(@minus, meanR, minValues);
maxValues(prefR) = nanmax(meanR(:, prefR));
maxValues(prefL) = nanmax(meanL(:, prefL));
meanL = bsxfun(@rdivide, meanL, maxValues);
meanR = bsxfun(@rdivide, meanR, maxValues);


nCellsL = numel(prefL);
nCellsR = numel(prefR);

% plot the snakes, while ordering the cells according to the peak latency
% in the overall mean response

clims = minmax([meanL(:); meanR(:)]');
clims = prctile([meanL(:); meanR(:)]', climPrctiles);
if separateColormap
    climsL = prctile(reshape([meanL(:, prefL), meanR(:, prefL)], 1, []), climPrctiles);
    climsR = prctile(reshape([meanL(:, prefR), meanR(:, prefR)], 1, []), climPrctiles);
else
    climsR = clims;
    climsL = clims;
end

subplot(4, 4, 9)
imagesc(zAxis, 1:nCellsL, meanL(:, order(ismember(order, prefL)))');
caxis(climsL);
title('Left trials');
ylabel('Cells Prefer Left');
axis xy
colorbar;

subplot(4, 4, 10)
imagesc(zAxis, 1:nCellsL, meanR(:, order(ismember(order, prefL)))');
caxis(climsL);
title('Right trials');
axis xy
colorbar

subplot(4, 4, 13)
imagesc(zAxis, 1:nCellsR, meanL(:, order(ismember(order, prefR)))');
caxis(climsR);
ylabel('Cells prefer Right')
axis xy
colorbar;

subplot(4, 4, 14)
imagesc(zAxis, 1:nCellsR, meanR(:, order(ismember(order, prefR)))');
caxis(climsR);
axis xy
colorbar

% And now order the cells 'within group': left-preferring cells by the left
% trial response latencies, right preferring cells - by right trial
% response latencies

[~, latency] = max(meanL(:, prefL));
[~, orderPrefL] = sort(latency);

[~, latency] = max(meanR(:, prefR));
[~, orderPrefR] = sort(latency);

subplot(4, 4, 11)
imagesc(zAxis, 1:nCellsL, meanL(:, prefL(orderPrefL))');
caxis(climsL);
title('Left trials');
ylabel('Cells Prefer Left');
axis xy
colorbar;

subplot(4, 4, 12)
imagesc(zAxis, 1:nCellsL, meanR(:, prefL(orderPrefL))');
caxis(climsL);
title('Right trials');
axis xy
colorbar

subplot(4, 4, 15)
imagesc(zAxis, 1:nCellsR, meanL(:, prefR(orderPrefR))');
caxis(climsR);
ylabel('Cells prefer Right')
axis xy
colorbar;

subplot(4, 4, 16)
imagesc(zAxis, 1:nCellsR, meanR(:, prefR(orderPrefR))');
caxis(climsR);
axis xy
colorbar

%%
% keyboard;
