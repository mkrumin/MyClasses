function rocAnalysis_SfN(obj)

%% Getting the trial indices

% keyboard
trialsAll = 1:length(obj.dataTMaze.contrastSequence);

cc = unique(obj.dataTMaze.contrastSequence);
nContrasts = length(cc);
for iCC = 1:nContrasts
    trialsGoR{iCC} = find(obj.dataTMaze.contrastSequence' == cc(iCC) & obj.dataTMaze.report == 'R');
    trialsGoL{iCC} = find(obj.dataTMaze.contrastSequence' == cc(iCC) & obj.dataTMaze.report == 'L');
    trialsC{iCC} = find(obj.dataTMaze.contrastSequence' == cc(iCC) & obj.dataTMaze.outcome == 'C');
    trialsW{iCC} = find(obj.dataTMaze.contrastSequence' == cc(iCC) & obj.dataTMaze.outcome == 'W');
    trialsStimR{iCC} = find(obj.dataTMaze.contrastSequence' == cc(iCC) & obj.dataTMaze.contrastSequence'>0);
    trialsStimL{iCC} = find(obj.dataTMaze.contrastSequence' == cc(iCC) & obj.dataTMaze.contrastSequence'<0);

end
% trialsGoR{length(cc)+1} = find(obj.dataTMaze.report == 'R');
% trialsGoL{length(cc)+1} = find(obj.dataTMaze.report == 'L');

[trialsGoR, groupLabels] = groupTrials(trialsGoR, cc);
trialsGoL = groupTrials(trialsGoL, cc);
trialsC = groupTrials(trialsC, cc);
trialsW = groupTrials(trialsW, cc);
trialsStimR = groupTrials(trialsStimR, cc);
trialsStimL = groupTrials(trialsStimL, cc);

nGroups = length(trialsGoR);
nPoints = 100;
iRow = 0;
cellIdx = [];
planeIdx = [];

%% Cutting F traces (trial-wise) for all the ROIs
fprintf('Analyzing experiment %s\n', obj.expRef);
for iPlane = obj.Planes
    fprintf('Plane %d: ', iPlane);
    tData = obj.times2p{iPlane}';
    nCells = obj.nROIs(iPlane);
    fMatrix{iPlane} = nan(length(trialsAll), nPoints, nCells);
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
        
        iRow = iRow + 1;
        
        [~, ~, fMatrix{iPlane}(:,:,iCell)] = ...
            buildTraces(obj, trialsAll, tData, fData, nPoints);
        
    end
    fprintf('\n');
end

fMatrix = cell2mat(reshape(fMatrix, 1, 1, []));
% keyboard;

%% Performing the decision-style analysis

nCells = size(fMatrix, 3);
fStd = [0.01 3];

nChars = 0;
for iCell = 1:nCells
    fprintf(repmat('\b', 1, nChars));
    nChars = fprintf('ROI %d/%d\n', iCell, nCells);
    %   take traces of one ROI
    data = fMatrix(:,:,iCell);
    %   filter the traces in time
    data = imgaussfilt(data, fStd);
    % normalize data to be between 0 and 1
    minF = min(data(:));
    maxF = max(data(:));
    data = (data - minF)/(maxF-minF);
    if all(all(isnan(data)))
        % this ROI traces were not made in the previous section
        continue;
    end
    
    for iGroup = 1:nGroups
        
        % calculate the trials' statistics
        tracesR = data(trialsGoR{iGroup}, :);
        tracesL = data(trialsGoL{iGroup}, :);
        
        tracesC = data(trialsC{iGroup}, :);
        tracesW = data(trialsW{iGroup}, :);
        
        tracesStimR = data(trialsStimR{iGroup}, :);
        tracesStimL = data(trialsStimL{iGroup}, :);

        % ROC analysis
        valuesR = nanmean(tracesR, 2);
        valuesL = nanmean(tracesL, 2);
        valuesR = valuesR(~isnan(valuesR));
        valuesL = valuesL(~isnan(valuesL));
        [rocRL(iCell, iGroup), tprRL{iCell, iGroup}, fprRL{iCell, iGroup}] = rocArea(valuesR, valuesL);

        valuesC = nanmean(tracesC, 2);
        valuesW = nanmean(tracesW, 2);
        valuesC = valuesC(~isnan(valuesC));
        valuesW = valuesW(~isnan(valuesW));
        [rocCW(iCell, iGroup), tprCW{iCell, iGroup}, fprCW{iCell, iGroup}] = rocArea(valuesC, valuesW);

        % U-test
        
        try
        [pValRL(iCell, iGroup), hRL(iCell, iGroup), stat] = ranksum(valuesR, valuesL);
        [pValCW(iCell, iGroup), hCW(iCell, iGroup), stat] = ranksum(valuesC, valuesW);
        end
%         valuesR = nanmean(tracesStimR, 2);
%         valuesL = nanmean(tracesStimL, 2);
%         valuesR = valuesR(~isnan(valuesR));
%         valuesL = valuesL(~isnan(valuesL));
%         [rocStimRL(iCell, iGroup), tprStimRL{iCell, iGroup}, fprStimRL{iCell, iGroup}] = rocArea(valuesR, valuesL);

    end
    
end
fprintf('\n');

%% plotting is done here

nRows = 2;
nColumns = nGroups;

nBins = 15;
binCentres = linspace(1/nBins/2, 1-1/nBins/2, nBins);

for iGroup = 1:nGroups
subplot(nRows, nColumns, iGroup);
hold off;

rocRL(rocRL==0) = NaN;
roc = rocRL(:, iGroup);
[rocAll] = hist(roc(~isnan(roc)), binCentres);
[rocSig] = hist(roc(hRL(:, iGroup)), binCentres);
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

rocCW(rocCW==0) = NaN;
roc = rocCW(:, iGroup);
[rocAll] = hist(roc(~isnan(roc)), binCentres);
[rocSig] = hist(roc(hCW(:, iGroup)), binCentres);
bar(binCentres, rocAll, 'w');
hold on;
bar(binCentres, rocSig, 'g');


xlim([0 1]);
plot([0.5 0.5], ylim, 'k:');
xlabel(sprintf('%1.0f/%1.0f significant', sum(rocSig), sum(rocAll)));
if iGroup == 1
    ylabel('C-W discrimination');
end
axis square

end

%%
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

