function rocAnalysis_DataVSModel(obj, res)

%% Getting the trial indices

% keyboard
trialsAll = 1:length(obj.dataTMaze.contrastSequence);

cc = unique(obj.dataTMaze.contrastSequence);
nContrasts = length(cc);
for iCC = 1:nContrasts
    trialsGoR{iCC} = find(obj.dataTMaze.contrastSequence' == cc(iCC) & obj.dataTMaze.report == 'R');
    trialsGoL{iCC} = find(obj.dataTMaze.contrastSequence' == cc(iCC) & obj.dataTMaze.report == 'L');
end

[trialsGoR, groupLabels] = groupTrials(trialsGoR, cc);
trialsGoL = groupTrials(trialsGoL, cc);

nGroups = length(trialsGoR);
nPoints = 100;
cellIdx = [];
planeIdx = [];

%% Cutting F traces (trial-wise) for all the ROIs
fprintf('Analyzing experiment %s\n', obj.expRef);
for iPlane = obj.Planes
    fprintf('Plane %d: ', iPlane);
    tData = obj.times2p{iPlane}';
    nCells = obj.nROIs(iPlane);
    fDataMatrix{iPlane} = nan(length(trialsAll), nPoints, nCells);
    fModelMatrix{iPlane} = nan(length(trialsAll), nPoints, nCells);
    nChars = 0;
    for iCell = 1:nCells
        fprintf(repmat('\b', 1, nChars));
        nChars = fprintf('ROI %d/%d', iCell, nCells);
        roiType = obj.data2p{iPlane}.ROI.CellClasses{iCell};
        if roiType ~= 's'
            % only process 'soma' ROIs
            continue;
        end
        
        fData = obj.data2p{iPlane}.F(:, iCell);
        if any(isnan(fData))
            continue;
        end
        
        cellIdx = cat(1, cellIdx, iCell);
        planeIdx = cat(1, planeIdx, iPlane);
        
        [~, ~, fDataMatrix{iPlane}(:,:,iCell)] = ...
            buildTraces(obj, trialsAll, tData, fData, nPoints);
        fModelMatrix{iPlane}(:,:,iCell) = estTraces(obj, res, trialsAll, iPlane, iCell, nPoints);
        
    end
    fprintf('\n');
end

fDataMatrix = cell2mat(reshape(fDataMatrix, 1, 1, []));
fModelMatrix = cell2mat(reshape(fModelMatrix, 1, 1, []));
% keyboard;

%% Performing the decision-style (ROC) analysis

nCells = size(fDataMatrix, 3);
fStd = [3];
hGauss = ndGaussian(fStd);

rocRLData = nan(nCells, nGroups);
rocRLModel = nan(nCells, nGroups);
pValRLData = nan(nCells, nGroups);
pValRLModel = nan(nCells, nGroups);

nChars = 0;
for iCell = 1:nCells
    fprintf(repmat('\b', 1, nChars));
    nChars = fprintf('ROI %d/%d\n', iCell, nCells);
    %   take traces of one ROI
    dataData = fDataMatrix(:,:,iCell);
    dataModel = fModelMatrix(:,:,iCell);
    %   filter the traces in time
    dataData = conv2(1, hGauss, dataData);
    dataModel = conv2(1, hGauss, dataModel);
    % normalize data to be between 0 and 1
    rangeF = minmax(dataData(:)');
    dataData = (dataData - rangeF(1))/diff(rangeF);
    rangeF = minmax(dataModel(:)');
    dataModel = (dataModel - rangeF(1))/diff(rangeF);

    if all(isnan(dataData(:))) || all(isnan(dataModel(:)))
        % this ROI traces were not made in the previous section
        continue;
    end
    
    for iGroup = 1:nGroups
        
        % calculate the trials' statistics
        tracesRData = dataData(trialsGoR{iGroup}, :);
        tracesLData = dataData(trialsGoL{iGroup}, :);
        tracesRModel = dataModel(trialsGoR{iGroup}, :);
        tracesLModel = dataModel(trialsGoL{iGroup}, :);
        
        % ROC analysis
        valuesRData = nanmean(tracesRData, 2);
        valuesLData = nanmean(tracesLData, 2);
        valuesRData = valuesRData(~isnan(valuesRData));
        valuesLData = valuesLData(~isnan(valuesLData));
        valuesRModel = nanmean(tracesRModel, 2);
        valuesLModel = nanmean(tracesLModel, 2);
        valuesRModel = valuesRModel(~isnan(valuesRModel));
        valuesLModel = valuesLModel(~isnan(valuesLModel));
        
        if ~any([isempty(valuesRData), isempty(valuesLData), isempty(valuesRModel), isempty(valuesLModel)])
            % roc statistics
            [rocRLData(iCell, iGroup), tprRLData{iCell, iGroup}, fprRLData{iCell, iGroup}] = rocArea(valuesRData, valuesLData);
            [rocRLModel(iCell, iGroup), tprRLModel{iCell, iGroup}, fprRLModel{iCell, iGroup}] = rocArea(valuesRModel, valuesLModel);
            
            % U-test
            [pValRLData(iCell, iGroup), hRLData(iCell, iGroup), ~] = ranksum(valuesRData, valuesLData);
            [pValRLModel(iCell, iGroup), hRLModel(iCell, iGroup), ~] = ranksum(valuesRModel, valuesLModel);
        end
        
    end % iGroup
    
end % iCell
fprintf('\n');

%% plotting is done here

nRows = 2;
nColumns = nGroups;

nBins = 15;
binCentres = linspace(1/nBins/2, 1-1/nBins/2, nBins);

% rocRLData(rocRLData==0) = NaN;
% rocRLModel(rocRLModel==0) = NaN;

for iGroup = 1:nGroups
    subplot(nRows, nColumns, iGroup);
    hold off;
    
    roc = rocRLData(:, iGroup);
    [rocAll] = hist(roc(~isnan(roc)), binCentres);
    [rocSig] = hist(roc(hRLData(:, iGroup)), binCentres);
    bar(binCentres, rocAll, 'w');
    hold on;
    bar(binCentres, rocSig, 'g');
    xlim([0 1]);
    ylim([0 50]); % this is a hack for a specific example
    plot([0.5 0.5], ylim, 'k:');
    title(groupLabels{iGroup});
    if iGroup == 1
        ylabel('R-L decision');
    end
    xlabel(sprintf('%1.0f/%1.0f significant', sum(rocSig), sum(rocAll)));
    if iGroup == 1
        legend('All', 'p<0.05');
    end
    axis square
    
    
    iRow = 2;
    subplot(nRows, nColumns, (iRow-1)*nColumns + iGroup);
    hold off;


    plot(rocRLData(:,iGroup), rocRLModel(:, iGroup), 'k.');
    hold on;
    
%     sigIdx = hRLData(:, iGroup);
%     nSigIdx = ~sigIdx;
%     plot(rocRLData(sigIdx,iGroup), rocRLModel(sigIdx, iGroup), 'g.');
%     hold on;
%     plot(rocRLData(nSigIdx,iGroup), rocRLModel(nSigIdx, iGroup), 'k.');
    
    validIdx = ~isnan(rocRLData(:,iGroup));
    rho = corrcoef(rocRLData(validIdx,iGroup), rocRLModel(validIdx, iGroup));
    rho = rho(2);
    plot([0 1], [0 1], 'k:');
    xlim([0 1]);
    ylim([0 1]);
    xlabel('d''_{Data}');
    if iGroup == 1
        ylabel('d''_{Model}');
        set(gca, 'YTick', [0 0.5 1]);
    end
    set(gca, 'YTick', [0 0.5 1], 'YTickLabel', []);
    title(sprintf('\\rho = %4.2f', rho));
    axis square
    
end

end
% ========================================================================
function [trialsOut, groupLabels] = groupTrials(trialsIn, cc)

groupLabels = {'High % L', 'Low % L', '0 %', 'Low % R', 'High % R', 'All %'};
bins = {[-50, -25];...
    [-12, -6];...
    [0];...
    [6, 12];...
    [25, 50];...
    cc};


% groupLabels = {'-50%', '-25%', '-12%', '-6%', '0%', '6%', '12%', '25%', '50%', 'All %'};
% bins = { -50; -25; -12; -6; 0; 6; 12; 25; 50; cc};


% groupLabels = {'0%', '\pm6%', '\pm12%', '\pm25%', '\pm50%', 'All %'};
% bins = {[0]; [-6, 6]; [-12, 12]; [-25, 25]; [-50, 50]; cc};

trialsIn = trialsIn(:)'; % making sure it is a row vector of cell (for cell2mat later)
nGroups = length(groupLabels);

trialsOut = cell(nGroups, 1);

for iGroup = 1:nGroups
    
    [~, ccInd, ~] = intersect(cc, bins{iGroup});
    trialsOut{iGroup} = unique(cell2mat(trialsIn(ccInd)));
    
end

end

% ========================================================================
function traces = estTraces(obj, res, trialIdx, iPlane, iROI, nPoints)

if nargin<6
    nPoints = 100;
end

tData = obj.times2p{iPlane};
fData = obj.data2p{iPlane}.F(:, iROI);

theMap = res{iPlane}(iROI).zThetaMap;
binCentres = res{iPlane}(iROI).zThetaBinCentres;

X = binCentres{2};
X(1) = X(1) - 0.001;
X(end) = X(end) + 0.001;
Y = binCentres{1};
Y(1) = Y(1) - 0.001;
Y(end) = Y(end) + 0.001;

[thMatrix, zMatrix, ~] = buildTraces(obj, 1:max(trialIdx), tData, fData, nPoints);
thMatrix = thMatrix(trialIdx, :);
zMatrix = zMatrix(trialIdx, :);

traces = nan(size(thMatrix));
for iTrial = 1:size(traces, 1)
    thValues = thMatrix(iTrial, :);
    zValues = zMatrix(iTrial, :);
    snippet = interp2(X, Y, theMap, thValues, zValues, 'linear');
    traces(iTrial,:) = snippet;
end

end % estTraces()
