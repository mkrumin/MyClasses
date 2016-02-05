function TrainMaps_v4(obj, iPlane, iROI, options)

% v4 of TrainMaps - the main difference from v3 is in the way the plots are
% made - prepared for Matteo's talk in Florence, July 2015
% in addition, predicted traces are calculated and some unnecessary
% calculations were taken out

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
% set(figHandle, 'Units', 'normalized', 'MenuBar', 'none')%, 'OuterPosition', [0 0 1 1]);
set(figHandle, 'Units', 'normalized');%, 'OuterPosition', [0 0 1 1]);

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
% traces = getTimeScaledETA(tData, fData, ...
%     [obj.timesTmazeOnset(trialIdx), obj.timesTmazeOffset(trialIdx)]);
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
% traces = cell(nDecisionGroups, nContrastGroups);
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

% now estimate the model-based trace

fDataModel = nan(size(fData));
fGlobalModel = nan(size(fData));
cTrace = nan(size(fData));
trialTrace = nan(size(fData));
for iValidTrial = 1:length(validTrials)
    iTrial = validTrials(iValidTrial);
    for iC = 1:length(options.contrastGroups)
        if ismember(obj.dataTMaze.contrastSequence(iTrial), options.contrastGroups{iC})
            break;
        end
    end
    [thValues, zValues, ~, tValues] = buildVectors(obj, iTrial, tData, fData);
    snippet = interp2(X, Y, theMap-baseline, thValues, zValues, 'linear');
    snippetGlobal = snippet*alpha(end, end)+beta(end, end);
    snippet = snippet*alpha(end, iC)+beta(end, iC);
    idx = find(tData>=tValues(1) & tData<=tValues(end));
    snippet = interp1(tValues, snippet, tData(idx), 'linear');
    snippetGlobal = interp1(tValues, snippetGlobal, tData(idx), 'linear');
    fDataModel(idx) = snippet;
    fGlobalModel(idx) = snippetGlobal;
    cTrace(idx) = obj.dataTMaze.contrastSequence(iTrial);
    trialTrace(idx) = iTrial;
end

% figure
cCmap = colormap(gcf, sprintf('jet(%d)', nContrastGroups-1));
% nColors = nContrastGroups-1;
% cCmap = [(0:nColors-1)/(nColors-1); 0.1*ones(1, nColors); (nColors-1:-1:0)/(nColors-1)]';
% cCmap = [(0:nColors-1)/(nColors-1); tripuls((-nColors+1)/2:(nColors-1)/2, nColors-1)/2; (nColors-1:-1:0)/(nColors-1)]';
colormap jet(256)
nRows = 4;
nColumns = nContrastGroups+1;

subplot(nRows, nColumns, [nColumns-1, nColumns]);
obj.ShowROI(iPlane, iROI);

hData = subplot(nRows, nColumns, [1:nColumns-2]);
plot(tData, fData, 'b');
hold on;
xlim([tData(1), tData(end)]);
yrange = [min([fData(:); fDataModel(:); fGlobalModel(:)]), max([fData(:); fDataModel(:); fGlobalModel(:)])];
ylim(yrange);
ylabel('F_{data}');


hModel = subplot(nRows, nColumns, [nColumns+1:2*nColumns-2]);
plot(tData, fDataModel, 'b', tData, fGlobalModel, 'r');
hold on;
xlim([tData(1), tData(end)]);
ylim(yrange);
ylabel('F_{model}');

pos = get(gca, 'Position');
hleg = legend('scaled', 'not scaled');
set(hleg, 'Location', 'NorthEast');
set(gca, 'Position', pos);

linkaxes([hData, hModel], 'xy');

for iC = 1:nContrastGroups-1
    trials = find(ismember(obj.dataTMaze.contrastSequence, options.contrastGroups{iC}));
    idx = find(ismember(trialTrace, trials));
    ySig = nan(size(tData));
    ySig(idx) = yrange(1);
    axes(hData);
    plot(tData, ySig, 'LineWidth', 3, 'Color', cCmap(iC, :));
    axes(hModel);
    plot(tData, ySig, 'LineWidth', 3, 'Color', cCmap(iC, :));
end
% zoom off;
% zoom xon;

maxF = max(fMatrix(:));
minF = min(fMatrix(:));
for iC = 1:nContrastGroups-1
    % picking the correct trials
    trialIdx = intersect(find(ismember(obj.dataTMaze.contrastSequence, options.contrastGroups{iC})), ...
        find(obj.dataTMaze.report=='L'));
    trialIdx = intersect(trialIdx, validTrials);
    if isempty(trialIdx)
        fRealL = [];
        fEstL = [];
    else
        fRealL = fMatrix(trialIdx, :);
        fMap = interp2(X, Y, theMap, thMatrix(trialIdx, :), zMatrix(trialIdx, :), 'linear');
        fEstL = alpha(end, iC)*(fMap-baseline)+beta(end, iC)+baseline;
    end
    
    trialIdx = intersect(find(ismember(obj.dataTMaze.contrastSequence, options.contrastGroups{iC})), ...
        find(obj.dataTMaze.report=='R'));
    trialIdx = intersect(trialIdx, validTrials);
    if isempty(trialIdx)
        fRealR = [];
        fEstR = [];
    else
        fRealR = fMatrix(trialIdx, :);
        fMap = interp2(X, Y, theMap, thMatrix(trialIdx, :), zMatrix(trialIdx, :), 'linear');
        fEstR = alpha(end, iC)*(fMap-baseline)+beta(end, iC)+baseline;
    end
    
    nTrials = size(fRealR, 1) + size(fRealL, 1);
    nanLine = nan(ceil(nTrials/10), size(thMatrix, 2));
    
    fReal = cat(1, fRealL, nanLine, fRealR);
    fEst = cat(1, fEstL, nanLine, fEstR);
    fReal = (fReal-minF)/(maxF-minF);
    fEst = (fEst-minF)/(maxF-minF);
    fReal(fReal<0) = 0;
    fReal(fReal>1) = 1;
    fEst(fEst<0) = 0;
    fEst(fEst>1) = 1;
    
    cmap = colormap;
    
    fReal = setNanColor(fReal, cmap, [1 1 1]);
    fEst = setNanColor(fEst, cmap, [1 1 1]);
    
    subplot(nRows, nColumns, sub2ind([nColumns, nRows], iC, 3))
    set(gca, 'Position', get(gca, 'Position'));
    image(fReal);
    axis off
    title(options.contrastGroupLabels{iC}, 'Color', cCmap(iC, :));
    if iC==1
        ylabel('F_{data}');
    end
    
    subplot(nRows, nColumns, sub2ind([nColumns, nRows], iC, 4))
    set(gca, 'Position', get(gca, 'Position'));
    image(fEst);
    axis off
    if iC==1
        ylabel('F_{model}');
    end
end


subIdx = [sub2ind([nColumns, nRows], nColumns-1, 2), sub2ind([nColumns, nRows], nColumns, 2)];
subplot(nRows, nColumns, subIdx);
plot(thAll, zAll, 'k.', 'MarkerSize', 4);
hold on;
imagesc(binCentres{2}, binCentres{1}, theMap, 'AlphaData', 0.8);
axis equal tight xy
caxis([minF, maxF]);
% title(sprintf('iPlane = %d, iROI = %d', iPlane, iROI));

subIdx = [sub2ind([nColumns, nRows], nColumns-1, nRows-1), sub2ind([nColumns, nRows], nColumns, nRows-1)];
subplot(nRows, nColumns, subIdx);
idx = find(~ismember(options.contrastGroupLabels, 'All'));
scaling = alpha(end, idx);
scaling(scaling<0) = nan;
plot(idx, log10(scaling), 'go-', 'LineWidth', 3);
% plot(idx, alpha(end, idx), 'ko-', 'LineWidth', 3);
hold on;

axis tight square
set(gca, 'Position', get(gca, 'Position'));
% legend(options.decisionGroupLabels, 'Location', 'Best');
set(gca, 'XTick', 1:nContrastGroups, 'XTickLabel', options.contrastGroupLabels, 'XTickLabelRotation', 90)
ylim([min(-1, min(log10(scaling))), max(1, max(log10(scaling)))]);
iC = find(ismember(options.contrastGroupLabels, '0%'));
if ~isempty(iC)
    plot([iC, iC], ylim, 'b:');
end
plot(xlim, log10([1 1]), 'b:');
title('Scaling (log_{10}(\alpha))');
% plot(xlim, [1 1], 'k:');
% title('Scaling factor');


subIdx = [sub2ind([nColumns, nRows], nColumns-1, nRows), sub2ind([nColumns, nRows], nColumns, nRows)];
subplot(nRows, nColumns, subIdx);
idx = ~isnan(fDataModel);
plot(fData, fDataModel, '.b', fData, fGlobalModel, '.r');
hold on;
axis equal
xlim([min([fData(idx); fDataModel(idx)]), max([fData(idx); fDataModel(idx)])]);
ylim(xlim);
box off
plot(xlim, ylim, 'k', 'LineWidth', 2);
tmp = corrcoef(fData(idx), fDataModel(idx));
rho = tmp(2);
tmp = corrcoef(fData(idx), fGlobalModel(idx));
rho2 = tmp(2);

pos = get(gca, 'Position');

hleg = legend('scaled', 'not scaled');
set(hleg, 'Location', 'EastOutside');

set(gca, 'Position', pos);

text(max(fData), min(yrange), ['\rho = ', num2str(rho, 2)],...
    'HorizontalAlignment', 'Right', 'VerticalAlignment', 'Bottom', 'FontSize', 14, 'Color', 'b')
text(min(fData), max(yrange), ['\rho = ', num2str(rho2, 2)],...
    'HorizontalAlignment', 'Left', 'VerticalAlignment', 'Top', 'FontSize', 14, 'Color', 'r')
xlabel('F_{Data}');
ylabel('F_{Model}');

set(gcf, 'Color', [0.7 0.7 0.7])

end % TrainMaps_v4()

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

%=====================================================
function imOut = setNanColor(im, cmap, nanColor)
nColors = length(cmap);
nanColor = nanColor(:)';
cmapFull = [cmap; nanColor];
im = im*(nColors-1);
im(isnan(im)) = nColors;
idx = floor(im);
imOut = cmapFull(idx(:)+1, :);
imOut = reshape(imOut, [size(im), 3]);
end