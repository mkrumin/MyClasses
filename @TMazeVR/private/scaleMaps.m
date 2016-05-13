function [fScaled] = scaleMaps(obj, fData, fModel)

cSequence = obj.dataTMaze.contrastSequence;
cValues = unique(cSequence);
[~, cIdx] = ismember(cSequence, cValues);
cIdx(~ismember(obj.dataTMaze.report, 'LR')) = 0;

[nBins, nCells, nTrials] = size(fData);

groups = unique(cIdx);
groups = groups(groups~=0);
nGroups = length(groups);
fScaled = nan(size(fModel));
%%
for iCell = 1:nCells
    nRows = floor(sqrt(nGroups+2));
    nColumns = ceil((nGroups+2)/nRows);
    subplot(nRows, nColumns, 1)
    for iGroup = 1:nGroups
        trialIdx = find(cIdx == groups(iGroup));
%         xx{iGroup} =  reshape(nanmean(fData(:, iCell, trialIdx)), [], 1);
%         yy{iGroup} =  reshape(nanmean(fModel(:, iCell, trialIdx)), [], 1);
        xx{iGroup} =  reshape(fData(:, iCell, trialIdx), [], 1);
        yy{iGroup} =  reshape(fModel(:, iCell, trialIdx), [], 1);
    end
    
    [coeffs, x0y0] = bilinfit(xx, yy);
    
    for iGroup = 1:nGroups
        trialIdx = find(cIdx == groups(iGroup));
        c = coeffs(iGroup, :);
        fScaled(:, iCell, trialIdx) = (fModel(:, iCell, trialIdx) -c(2))/c(1);
        
    end
    
end

return;

%% plotting and testing and debugging
for iGroup = 1:nGroups
    plot(xx{iGroup}, coeffs(iGroup, 1)*xx{iGroup}+coeffs(iGroup, 2), 'Color', h(iGroup).Color);
end
axis tight
xlimMain = xlim;
plot(xlimMain, xlimMain, ':k')
hold off;

for iGroup = 1:nGroups
    subplot(nRows, nColumns, iGroup+1)
    cla
    plot(xx{iGroup}, yy{iGroup}, 'o', 'Color', h(iGroup).Color)
    hold on;
    axis tight
    plot(xlimMain, xlimMain, ':k')
    xRange = [nanmin(xx{iGroup}), nanmax(xx{iGroup})];
    plot(xRange, coeffs(iGroup, 1)*xRange+coeffs(iGroup, 2), 'Color', h(iGroup).Color);
    pp = polyfit(xx{iGroup}, yy{iGroup}, 1);
    yyNew{iGroup} = (yy{iGroup}-coeffs(iGroup, 2))/coeffs(iGroup, 1);
    yyNewLocal{iGroup} = (yy{iGroup}-pp(2))/pp(1);
    plot(xRange, pp(1)*xRange+pp(2), ':', 'Color', h(iGroup).Color, 'LineWidth', 2);
    plot(xx{iGroup}, yyNew{iGroup}, '.k');
    plot(xx{iGroup}, yyNewLocal{iGroup}, '.', 'Color', h(iGroup).Color);
    
    %     hold off
    
    %     plot(xRange, pp(1)*xRange+pp(2), 'k');%, 'Color', h(iGroup).Color, 'LineWidth', 2);
    %     hold on;
    %     plot(xRange, coeffs(iGroup, 1)*xRange+coeffs(iGroup, 2), 'r');%'Color', h(iGroup).Color);
    
end
subplot(nRows, nColumns, nGroups+2);
figure
plot(cell2mat(xx'), cell2mat(yy'), '.k')
hold on;
plot(cell2mat(xx'), cell2mat(yyNew'), '.r')
plot(cell2mat(xx'), cell2mat(yyNewLocal'), '.c')

xlim(xlimMain)
ylim(xlimMain)
plot(xlimMain, xlimMain, ':k')
corrcoef(cell2mat(xx'), cell2mat(yy'))
corrcoef(cell2mat(xx'), cell2mat(yyNew'))
corrcoef(cell2mat(xx'), cell2mat(yyNewLocal'))


hold off

