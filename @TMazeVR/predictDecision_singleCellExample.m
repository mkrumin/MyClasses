function predictDecision_singleCellExample(obj, iPlane, iCell)

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
end
% trialsGoR{length(cc)+1} = find(obj.dataTMaze.report == 'R');
% trialsGoL{length(cc)+1} = find(obj.dataTMaze.report == 'L');

[trialsGoR, groupLabels] = groupTrials(trialsGoR, cc);
trialsGoL = groupTrials(trialsGoL, cc);
trialsC = groupTrials(trialsC, cc);
trialsW = groupTrials(trialsW, cc);

nGroups = length(trialsGoR);
nPoints = 100;
iRow = 0;
cellIdx = [];
planeIdx = [];

%% Cutting F traces (trial-wise) for all the ROIs
fprintf('Plane %d: ', iPlane);
tData = obj.times2p{iPlane}';
fMatrix = nan(length(trialsAll), nPoints);
thMatrix = nan(length(trialsAll), nPoints);
zMatrix = nan(length(trialsAll), nPoints);

roiType = obj.data2p{iPlane}.ROI.CellClasses{iCell};

fData = obj.data2p{iPlane}.F(:, iCell);

[thMatrix, zMatrix, fMatrix] = ...
    buildTraces(obj, trialsAll, tData, fData, nPoints);
fMatrix = estTracesFromMap(obj, iPlane, iCell, nPoints);
% keyboard;

%% Performing the decision-style analysis

fStd = [0.01 3];

%   take traces of one ROI
data = fMatrix;
%   filter the traces in time
data = imgaussfilt(data, fStd);
F0 = prctile(data(:), 10);
data = (data-F0)/F0;
% normalize data to be between 0 and 1
% minF = min(data(:));
% maxF = max(data(:));
% data = (data - minF)/(maxF-minF);

% theta and z should be the same for all cells (on the same plane), the
% code can be potentially optimized
th = thMatrix;
z = zMatrix;
maxTh = max(abs(th(:)));
minZ = min(z(:));
maxZ = max(z(:));

for iGroup = 1:nGroups
    
    % get the trajectories of the groups
    thR{iGroup} = th(trialsGoR{iGroup}, :);
    thL{iGroup} = th(trialsGoL{iGroup}, :);
    zR{iGroup} = z(trialsGoR{iGroup}, :);
    zL{iGroup} = z(trialsGoL{iGroup}, :);
    
    thC{iGroup} = th(trialsC{iGroup}, :);
    thW{iGroup} = th(trialsW{iGroup}, :);
    zC{iGroup} = z(trialsC{iGroup}, :);
    zW{iGroup} = z(trialsW{iGroup}, :);
    
    % calculate the trials' statistics
    tracesR = data(trialsGoR{iGroup}, :);
    tracesL = data(trialsGoL{iGroup}, :);
    meanR(iGroup, :) = nanmean(tracesR, 1)';
    nTrialsR(iGroup) = sum(~isnan(tracesR(:,1)));
    semR(iGroup, :) = nanstd(tracesR, [], 1)'/sqrt(nTrialsR(iGroup));
    meanL(iGroup, :) = nanmean(tracesL, 1)';
    nTrialsL(iGroup) = sum(~isnan(tracesL(:,1)));
    semL(iGroup, :) = nanstd(tracesL, [], 1)'/sqrt(nTrialsL(iGroup));
    
    tracesC = data(trialsC{iGroup}, :);
    tracesW = data(trialsW{iGroup}, :);
    meanC(iGroup, :) = nanmean(tracesC, 1)';
    nTrialsC(iGroup) = sum(~isnan(tracesC(:,1)));
    semC(iGroup, :) = nanstd(tracesC, [], 1)'/sqrt(nTrialsC(iGroup));
    meanW(iGroup, :) = nanmean(tracesW, 1)';
    nTrialsW(iGroup) = sum(~isnan(tracesW(:,1)));
    semW(iGroup, :) = nanstd(tracesW, [], 1)'/sqrt(nTrialsW(iGroup));
    
    
    % ROC analysis
    valuesR = nanmean(tracesR, 2);
    valuesL = nanmean(tracesL, 2);
    valuesR = valuesR(~isnan(valuesR));
    valuesL = valuesL(~isnan(valuesL));
    [rocRL(iGroup), tprRL{iGroup}, fprRL{iGroup}] = rocArea(valuesR, valuesL);
    
    valuesC = nanmean(tracesC, 2);
    valuesW = nanmean(tracesW, 2);
    valuesC = valuesC(~isnan(valuesC));
    valuesW = valuesW(~isnan(valuesW));
    [rocCW(iGroup), tprCW{iGroup}, fprCW{iGroup}] = rocArea(valuesC, valuesW);

        % U-test
    [pValRL(iGroup), hRL(iGroup), stat] = ranksum(valuesR, valuesL);
    [pValCW(iGroup), hCW(iGroup), stat] = ranksum(valuesC, valuesW);

end

%% plotting is done here

pcAxes = subplot(1, 3, 1);
showPC(obj.dataTMaze, pcAxes);
expRefStr = strrep(obj.expRef, '_', '\_');
title(sprintf('%s\niROI #%d', expRefStr, iCell));

for iGroup = 1:nGroups

    subplot(1, 3, 2)
    hold off;
    binCentres = linspace(min([valuesL(:); valuesR(:)]), max([valuesL(:); valuesR(:)]), 30);
    [hL, bL] = hist(valuesL, binCentres, 'b');
    [hR, bR] = hist(valuesR, binCentres, 'r');
    bar(bL, -hL, 'b');
    hold on;
    bar(bR, hR, 'r');
    title(sprintf('p = %d', pValRL(iGroup)));
    axis square tight
    
    medL = median(valuesL);
    medR = median(valuesR);
    plot([medR, medR], [0, max(ylim)], 'k:', 'LineWidth', 3)
    plot([medL, medL], [0, min(ylim)], 'k:', 'LineWidth', 3)
    
    subplot(1, 3, 3);
    cla;
    hold off;
    xx = [fprRL{iGroup}, 1, 0];
    yy = [tprRL{iGroup}, 0, 0];
    patch(xx, yy, 'b', 'FaceAlpha', 0.05, 'LineStyle', 'none')
    hold on;
    plot(fprRL{iGroup}, tprRL{iGroup}, 'b', 'LineWidth', 2);
    plot([0 1], [0 1], 'k:')
    
    fErr = (fprRL{iGroup}-0.5).^2;
    medFInd = find(fErr==min(fErr), 1, 'first');
    tErr = (tprRL{iGroup}-0.5).^2;
    medTInd = find(tErr==min(tErr), 1, 'first');
    
    plot([0, fprRL{iGroup}(medTInd)], [0.5 0.5], 'k:', 'LineWidth', 2)
    plot(repmat(fprRL{iGroup}(medTInd), 1, 2), [0, 0.5], 'k:', 'LineWidth', 2)
    plot([0.5 0.5], [0, tprRL{iGroup}(medFInd)], 'k:', 'LineWidth', 2)
    plot([0, 0.5], repmat(tprRL{iGroup}(medFInd), 1, 2), 'k:', 'LineWidth', 2)
    hold off;
    title(sprintf('ROC = %4.2f', rocRL(iGroup)))
    axis equal tight
    
end




%%
function [trialsOut, groupLabels] = groupTrials(trialsIn, cc)

% groupLabels = {'High % L', 'Low % L', '0 %', 'Low % R', 'High % R', 'All %'};
% bins = [-100, -25;...
%     -12, -6;...
%     0, 0;...
%     6, 12;...
%     25, 100;...
%     -Inf, Inf];
%

% groupLabels = {'-50%', '-25%', '-12%', '-6%', '0%', '6%', '12%', '25%', '50%', 'All %'};
% bins = [-50, -50;...
%     -25, -25;...
%     -12, -12;...
%     -6, -6;...
%     0, 0;...
%     6, 6;...
%     12, 12;...
%     25, 25;...
%     50, 50;...
%     -Inf, Inf];

groupLabels = {'All %'};
bins = [-Inf, Inf];

trialsIn = trialsIn(:)'; % making sure it is a row vector of cell (for cell2mat later)
nGroups = length(groupLabels);

trialsOut = cell(nGroups, 1);

for iGroup = 1:nGroups
    
    ccInd = cc>=bins(iGroup, 1) & cc<=bins(iGroup, 2);
    trialsOut{iGroup} = unique(cell2mat(trialsIn(ccInd)));
    
end

