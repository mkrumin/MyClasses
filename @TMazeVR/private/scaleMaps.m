function [fScaled, alpha] = scaleMaps(fData, fModel, cIdx)

[nBins, nCells, nTrials] = size(fData);

groups = unique(cIdx);
nGroups = length(groups);
%%
iCell = 7;
for iGroup = 1:nGroups
trialIdx = find(cIdx == groups(iGroup));
plot(reshape(nanmean(fData(:, iCell, trialIdx)), [], 1), ...
    reshape(nanmean(fModel(:, iCell, trialIdx)), [], 1), '.')
hold on;
end
axis equal
hold off;
%%
for iCell = 1:nCells
end
