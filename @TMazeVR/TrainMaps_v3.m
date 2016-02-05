function TrainMaps_v3(obj, iPlane, iROI, options)

if nargin<4
    options = struct;
end
options = fillOptions(obj, options);

nContrastGroups = length(options.contrastGroups);
nDecisionGroups = length(options.decisionGroups);
nChunks = options.cvFactor;

fData = obj.data2p{iPlane}.F(:, iROI);
tData = obj.times2p{iPlane};

figHandle = gcf;
clf(figHandle, 'reset');
set(figHandle, 'Name', sprintf('ROI # %d processing...', iROI));
drawnow;
set(figHandle, 'Units', 'normalized', 'MenuBar', 'none')%, 'OuterPosition', [0 0 1 1]);

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
traces = getTimeScaledETA(tData, fData, ...
    [obj.timesTmazeOnset(trialIdx), obj.timesTmazeOffset(trialIdx)]);
[thMatrix, zMatrix, fMatrix] = buildTraces(obj, 1:length(obj.dataTMaze.contrastSequence), tData, fData);
% subplot(1, 2, 1)
% plot(thAll, zAll, '.k', 'MarkerSize', 4);
% hold on;
% imagesc(binCentres{2}, binCentres{1}, theMap, 'AlphaData', 0.8);
% axis xy tight off
% subplot(1, 2, 2)
% imagesc(traces.t, 1:size(traces.all, 1), traces.all);
% axis tight off
% colormap jet(256)


% fMap = cell(nContrastGroups, options.cvFactor);
thVector = cell(nDecisionGroups, nContrastGroups);
zVector = cell(nDecisionGroups, nContrastGroups);
fVector = cell(nDecisionGroups, nContrastGroups);
traces = cell(nDecisionGroups, nContrastGroups);
% figure
% colormap jet(256)
for iD = 1:nDecisionGroups
    for iC = 1:nContrastGroups
        % picking the correct trials
        trialIdx = intersect(find(ismember(obj.dataTMaze.contrastSequence, options.contrastGroups{iC})), ...
            find(ismember(obj.dataTMaze.report, options.decisionGroups{iD})));
        trialIdx = intersect(trialIdx, validTrials);
        if isempty(trialIdx)
            continue;
        end
        
%         traces{iD, iC} = getTimeScaledETA(tData, fData, ...
%             [obj.timesTmazeOnset(trialIdx), obj.timesTmazeOffset(trialIdx)]);
        
        [thVector{iD, iC}, zVector{iD, iC}, fVector{iD, iC}] = ...
            buildVectors(obj, trialIdx, tData, fData);

%         subInd = sub2ind([nContrastGroups, nDecisionGroups], iC, iD);
%         subplot(nDecisionGroups, nContrastGroups, subInd)
%         imagesc(traces{iD, iC}.all);
%         caxis([min(fAll), max(fAll)]);
%         title(options.contrastGroupLabels{iC});
%         ylabel(options.decisionGroupLabels{iD});
    end
end

X = binCentres{2};
X(1) = X(1) - 0.001;
X(end) = X(end) + 0.001;
Y = binCentres{1};
Y(1) = Y(1) - 0.001;
Y(end) = Y(end) + 0.001;

% baseline = median(theMap(:));
% baseline = mode(theMap(:));
% baseline = median(fAll);
baseline = 0;
alpha = nan(nDecisionGroups, nContrastGroups);
beta = zeros(nDecisionGroups, nContrastGroups);
err = nan(nDecisionGroups, nContrastGroups);
rho = nan(nDecisionGroups, nContrastGroups);
for iC = 1:nContrastGroups
    for iD = 1:nDecisionGroups
        if isempty(thVector{iD, iC})
            continue;
        end
        fMap = interp2(X, Y, theMap-baseline, thVector{iD, iC}, zVector{iD, iC}, 'linear');
        fReal = fVector{iD, iC}-baseline;%median(fAll);
        % estimation of the optimal scaling parameter
%         alpha(iD, iC) = sum(fReal.*fMap)/sum(fMap.^2);
        tmp = pinv([fMap, ones(size(fMap))])*fReal;
        alpha(iD, iC) = tmp(1);
        beta(iD, iC) = tmp(2);
        fEst = alpha(iD, iC)*fMap+beta(iD, iC);
        err(iD, iC) = mean((fReal-fEst).^2)/var(fReal);
        cc = corrcoef(fReal(:), fEst(:));
        rho(iD, iC) = cc(2);
        
        fGlobalModel{iD, iC} = fMap;
        fBestModel{iD, iC} = fEst;
        fActualData{iD, iC} = fReal;
    end
end

% figure
colormap jet(256)
nRows = nDecisionGroups*2;
nColumns = nContrastGroups+2;
maxF = max(fMatrix(:));
minF = min(fMatrix(:));
for iD = 1:nDecisionGroups
    for iC = 1:nContrastGroups
        % picking the correct trials
        trialIdx = intersect(find(ismember(obj.dataTMaze.contrastSequence, options.contrastGroups{iC})), ...
            find(ismember(obj.dataTMaze.report, options.decisionGroups{iD})));
        trialIdx = intersect(trialIdx, validTrials);
        if isempty(trialIdx)
            continue;
        end
        
        fMap = interp2(X, Y, theMap, thMatrix(trialIdx, :), zMatrix(trialIdx, :), 'linear');
        fEst = alpha(iD, iC)*(fMap-baseline)+beta(iD, iC)+baseline;
        
        subplot(nRows, nColumns, sub2ind([nColumns, nRows], iC, 2*iD-1))
        set(gca, 'Position', get(gca, 'Position'));
        imagesc(fMatrix(trialIdx, :));
        caxis([minF, maxF]);
        if iD==1
            title(options.contrastGroupLabels{iC});
        end
        axis off

        subplot(nRows, nColumns, sub2ind([nColumns, nRows], iC, 2*iD))
        set(gca, 'Position', get(gca, 'Position'));
        imagesc(fEst);
        caxis([minF, maxF]);
%         title({sprintf('\\alpha=%4.2f \\epsilon=%4.2f', alpha(iD, iC), err(iD, iC));...
%             sprintf('\\rho=%4.2f', rho(iD, iC))}, 'FontSize', 8);
        axis off
    end
    pos = get(gca, 'Position');
    annotation('textbox', [0.01, pos(2)+pos(4)/2, 0.1, 0.1], 'String', options.decisionGroupLabels{iD}, 'FontSize', 20);
    annotation('textbox', [0.09, pos(2)+pos(4), 0.1, 0.1], 'String', 'Data', 'FontSize', 12, 'LineStyle', 'none');
    annotation('textbox', [0.09, pos(2)-0.03, 0.1, 0.1], 'String', 'Model', 'FontSize', 12, 'LineStyle', 'none');
    if iD<nDecisionGroups
        annotation('line', [0 1], [pos(2) pos(2)]-0.02);
    end
    
end
% pos = get(gca, 'Position');
% annotation('line', [pos(1)+pos(3) pos(1)+pos(3)]+0.01, [0 1]);

cIdx = find(~ismember(options.contrastGroupLabels, 'All'));
allCInd = find(ismember(options.contrastGroupLabels, 'All'));
dIdx = find(~ismember(options.decisionGroupLabels, 'All'));
allDInd = find(ismember(options.decisionGroupLabels, 'All'));

% error of the original cross-validated map
fModel = cell2mat(reshape(fGlobalModel(allDInd, allCInd), [], 1));
fReal = cell2mat(reshape(fActualData(allDInd, allCInd), [], 1));
errOneMap = mean((fModel-fReal).^2)/var(fReal);
% error of the SCALED original cross-validated map
fModel = cell2mat(reshape(fBestModel(allDInd, allCInd), [], 1));
fReal = cell2mat(reshape(fActualData(allDInd, allCInd), [], 1));
errOneMapScaled = mean((fModel-fReal).^2)/var(fReal);
% error of the model, which takes into account different contrasts
fModel = cell2mat(reshape(fBestModel(allDInd, cIdx), [], 1));
fReal = cell2mat(reshape(fActualData(allDInd, cIdx), [], 1));
errC = mean((fModel-fReal).^2)/var(fReal);
% error of the model, which takes into account different contrasts AND
% different behaviours
fModel = cell2mat(reshape(fBestModel(dIdx, cIdx), [], 1));
fReal = cell2mat(reshape(fActualData(dIdx, cIdx), [], 1));
errCnD = mean((fModel-fReal).^2)/var(fReal);


subIdx = [sub2ind([nColumns, nRows], nColumns-1, 1), sub2ind([nColumns, nRows], nColumns, 2)];
subplot(nRows, nColumns, subIdx);
plot(thAll, zAll, 'k.', 'MarkerSize', 4);
hold on;
imagesc(binCentres{2}, binCentres{1}, theMap, 'AlphaData', 0.8);
axis tight xy
caxis([minF, maxF]);
title(sprintf('iPlane = %d, iROI = %d', iPlane, iROI));

subIdx = [sub2ind([nColumns, nRows], nColumns-1, 3), sub2ind([nColumns, nRows], nColumns, 4)];
subplot(nRows, nColumns, subIdx);
obj.dataTMaze.showPC(gca);
delete(get(gca, 'YLabel'));
ch = get(gcf, 'Children');
for iCh=1:length(ch)
    if isequal(get(ch(iCh), 'Type'), 'legend')
        delete(ch(iCh));
        break;
    end
end

subIdx = [sub2ind([nColumns, nRows], nColumns-1, nRows-1), sub2ind([nColumns, nRows], nColumns, nRows)];
subplot(nRows, nColumns, subIdx);
linecolors = 'rbk';
idx = find(~ismember(options.contrastGroupLabels, 'All'));
allInd = find(ismember(options.contrastGroupLabels, 'All'));
for iD = 1:nDecisionGroups
%     plot(idx, alpha(iD, idx), linecolors(iD), 'LineWidth', 3);
%     semilogy(idx, alpha(iD, idx), linecolors(iD), 'LineWidth', 3);
    plot(idx, log10(alpha(iD, idx)), [linecolors(iD) 'o-'], 'LineWidth', 3);

    hold on;
end
for iD = 1:nDecisionGroups
%     plot(allInd, alpha(iD, allInd), 'o', 'Color', linecolors(iD));
%     semilogy(allInd, alpha(iD, allInd), 'o', 'Color', linecolors(iD));
    plot(allInd, log10(alpha(iD, allInd)), 'o', 'Color', linecolors(iD));
end
axis tight
set(gca, 'Position', get(gca, 'Position'));
legend(options.decisionGroupLabels, 'Location', 'Best');
set(gca, 'XTick', 1:nContrastGroups, 'XTickLabel', options.contrastGroupLabels, 'XTickLabelRotation', 90)
iC = find(ismember(options.contrastGroupLabels, '0%'));
if ~isempty(iC)
    plot([iC, iC], ylim, 'k:');
end
plot(xlim, log10([1 1]), 'k:');
title('Scaling (log_{10}(\alpha))');

annotation('textbox', [0.91, 0.8, 0.1, 0.1], 'String', sprintf('\\epsilon_{Global} = %5.3f', errOneMap), 'FontSize', 14);
annotation('textbox', [0.91, 0.75, 0.1, 0.1], 'String', sprintf('\\epsilon_{Scaled} = %5.3f', errOneMapScaled), 'FontSize', 14);
annotation('textbox', [0.91, 0.7, 0.1, 0.1], 'String', sprintf('\\epsilon_{C} = %5.3f', errC), 'FontSize', 14);
annotation('textbox', [0.91, 0.65, 0.1, 0.1], 'String', sprintf('\\epsilon_{CandD} = %5.3f', errCnD), 'FontSize', 14);


end % TrainMaps_v3()

%==========================================================================

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
