function expVar = estExplainedVariance(obj, res)

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

%% Estimating the explained variance (mean signal per trial)

nCells = size(fDataMatrix, 3);

nChars = 0;

meanData = squeeze(nanmean(fDataMatrix, 2));
meanModel = squeeze(nanmean(fModelMatrix, 2));
rho = nan(nCells, 1);
expVar = nan(nCells, 1);

for iCell = 1:nCells
    fprintf(repmat('\b', 1, nChars));
    nChars = fprintf('ROI %d/%d\n', iCell, nCells);
    %   take traces of one ROI
    dataData = meanData(:,iCell);
    dataModel = meanModel(:,iCell);
    % normalize data to be between 0 and 1
    
    validIdx = ~(isnan(dataData) | isnan(dataModel));
    if sum(validIdx)>2
        tmp = corrcoef(dataData(validIdx), dataModel(validIdx));
        rho(iCell) = tmp(2);
        
        c = polyfit(dataModel(validIdx), dataData(validIdx), 1);
        dataModel = c(2)+c(1)*dataModel;
        expVar(iCell) = 1 - sum((dataModel(validIdx)-dataData(validIdx)).^2)/...
            sum((nanmean(dataModel(validIdx))-dataData(validIdx)).^2);
    end
    
end % iCell
fprintf('\n');

%% plotting is done here

figure;
subplot(1, 2, 1)
hist(expVar, 20);
xlabel('Explained Variance');
subplot(1, 2, 2)
hist(rho, 20);
xlabel('rho');

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
