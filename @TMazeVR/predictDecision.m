function predictDecision(obj)

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
for iPlane = obj.Planes
    fprintf('Plane %d: ', iPlane);
    tData = obj.times2p{iPlane}';
    nCells = obj.nROIs(iPlane);
    fMatrix{iPlane} = nan(length(trialsAll), nPoints, nCells);
    thMatrix{iPlane} = nan(length(trialsAll), nPoints, nCells);
    zMatrix{iPlane} = nan(length(trialsAll), nPoints, nCells);
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
        
        [thMatrix{iPlane}(:,:,iCell), zMatrix{iPlane}(:,:,iCell), fMatrix{iPlane}(:,:,iCell)] = ...
            buildTraces(obj, trialsAll, tData, fData, nPoints);
        
    end
    fprintf('\n');
end

fMatrix = cell2mat(reshape(fMatrix, 1, 1, []));
thMatrix = cell2mat(reshape(thMatrix, 1, 1, []));
zMatrix = cell2mat(reshape(zMatrix, 1, 1, []));
% keyboard;

%% Performing the decision-style analysis

nCells = size(fMatrix, 3);
fStd = [0.01 3];

for iCell = 1:nCells
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
    
    % theta and z should be the same for all cells (on the same plane), the
    % code can be potentially optimized
    th = thMatrix(:,:,iCell);
    z = zMatrix(:,:,iCell);
    maxTh = max(abs(thMatrix(:)));
    minZ = min(zMatrix(:));
    maxZ = max(zMatrix(:));
    
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
        
        % U-test
        [pValRL(iGroup), hRL(iGroup), stat] = ranksum(nanmean(tracesR, 2), nanmean(tracesL, 2));
        [pValCW(iGroup), hCW(iGroup), stat] = ranksum(nanmean(tracesC, 2), nanmean(tracesW, 2));
        
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
        
    end
    
    %% plotting is done here
    
    nRows = 5;
    nColumns = nGroups+2;
    
    maxYRL = max(max(meanR(:)+semR(:)), max(meanL(:)+semL(:)));
    minYRL = min(min(meanR(:)-semR(:)), min(meanL(:)-semL(:)));

    maxYCW = max(max(meanC(:)+semC(:)), max(meanW(:)+semW(:)));
    minYCW = min(min(meanC(:)-semC(:)), min(meanW(:)-semW(:)));
    
    maxY = max(maxYRL, maxYCW);
    minY = min(minYRL, minYCW);

    iRow = 2;
    pcAxes = subplot(nRows, nColumns, nColumns*(iRow-1)+nGroups + [1,2]);
    showPC(obj.dataTMaze, pcAxes);
    expRefStr = strrep(obj.expRef, '_', '\_');
    title(sprintf('%s\niROI #%d', expRefStr, iCell));
    
    for iGroup = 1:nGroups
        
        iRow = 1;
        subplot(nRows, nColumns, nColumns*(iRow-1)+iGroup);
        
        plot(thR{iGroup}', zR{iGroup}', 'r:', thL{iGroup}', zL{iGroup}', 'b:');
        xlim([-maxTh maxTh]);
        ylim([minZ, maxZ]);
        title(groupLabels{iGroup})
        
        iRow = 2;
        xx = [1:nPoints, nPoints:-1:1]';
        yyR = [meanR(iGroup, :) + semR(iGroup, :), fliplr(meanR(iGroup, :) - semR(iGroup, :))];
        yyL = [meanL(iGroup, :) + semL(iGroup, :), fliplr(meanL(iGroup, :) - semL(iGroup, :))];
        
        subplot(nRows, nColumns, nColumns*(iRow-1)+iGroup);
        cla;
        patch(xx, yyR, 'r', 'FaceAlpha', 0.5, 'LineStyle', 'none');
        hold on;
        patch(xx, yyL, 'b', 'FaceAlpha', 0.5, 'LineStyle', 'none');
        plot(1:nPoints, meanR(iGroup, :), 'r', 'LineWidth', 2);
        plot(1:nPoints, meanL(iGroup, :), 'b', 'LineWidth', 2);
        hold off;
        
        xlim([1, 100]);
        ylim([minY, maxY]);
        warning off;
        set(gca, 'XTickLabels', '', 'YTickLabels', '');
        warning on;
        text(min(xlim), max(ylim), [' ', num2str(nTrialsL(iGroup))], ...
            'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', 'FontSize', 14, 'Color', 'b');
        text(max(xlim), max(ylim), [num2str(nTrialsR(iGroup)), ' '], ...
            'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', 'FontSize', 14, 'Color', 'r');
        
        xlabel(sprintf('%d (%3.1d)', uint8(hRL(iGroup)), pValRL(iGroup)));
        
        iRow = 3;
        subplot(nRows, nColumns, iGroup+(iRow-1)*nColumns);
        plot(fprRL{iGroup}, tprRL{iGroup}, 'b', 'LineWidth', 2);
        hold on;
        plot([0 1], [0 1], 'k:')
        hold off;
        title(sprintf('ROC = %4.2f', rocRL(iGroup)))
        
        iRow = 4;
        xx = [1:nPoints, nPoints:-1:1]';
        yyC = [meanC(iGroup, :) + semC(iGroup, :), fliplr(meanC(iGroup, :) - semC(iGroup, :))];
        yyW = [meanW(iGroup, :) + semW(iGroup, :), fliplr(meanW(iGroup, :) - semW(iGroup, :))];
        
        subplot(nRows, nColumns, nColumns*(iRow-1)+iGroup);
        cla;
        patch(xx, yyC, 'r', 'FaceAlpha', 0.5, 'LineStyle', 'none');
        hold on;
        patch(xx, yyW, 'b', 'FaceAlpha', 0.5, 'LineStyle', 'none');
        plot(1:nPoints, meanC(iGroup, :), 'r', 'LineWidth', 2);
        plot(1:nPoints, meanW(iGroup, :), 'b', 'LineWidth', 2);
        hold off;
        
        xlim([1, 100]);
        ylim([minY, maxY]);
        warning off;
        set(gca, 'XTickLabels', '', 'YTickLabels', '');
        warning on;
        text(min(xlim), max(ylim), [' ', num2str(nTrialsW(iGroup))], ...
            'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', 'FontSize', 14, 'Color', 'b');
        text(max(xlim), max(ylim), [num2str(nTrialsC(iGroup)), ' '], ...
            'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', 'FontSize', 14, 'Color', 'r');
        
        xlabel(sprintf('%d (%3.1d)', uint8(hCW(iGroup)), pValCW(iGroup)));
        
        iRow = 5;
        subplot(nRows, nColumns, iGroup+(iRow-1)*nColumns);
        plot(fprCW{iGroup}, tprCW{iGroup}, 'b', 'LineWidth', 2);
        hold on;
        plot([0 1], [0 1], 'k:')
        hold off;
        title(sprintf('ROC = %4.2f', rocCW(iGroup)))
        
    end
    
    
    drawnow;
    pause;
    
end

return;

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
groupLabels = {'-50%', '-25%', '-12%', '-6%', '0%', '6%', '12%', '25%', '50%', 'All %'};
bins = [-50, -50;...
    -25, -25;...
    -12, -12;...
    -6, -6;...
    0, 0;...
    6, 6;...
    12, 12;...
    25, 25;...
    50, 50;...
    -Inf, Inf];

trialsIn = trialsIn(:)'; % making sure it is a row vector of cell (for cell2mat later)
nGroups = length(groupLabels);

trialsOut = cell(nGroups, 1);

for iGroup = 1:nGroups
    
    ccInd = cc>=bins(iGroup, 1) & cc<=bins(iGroup, 2);
    trialsOut{iGroup} = unique(cell2mat(trialsIn(ccInd)));
    
end

