function accMap = buildAccumMap(coords, signal, binEdges)

nDims = length(binEdges);
for iDim=1:nDims
    dBin = mean(diff(binEdges{iDim}));
    nBins(iDim) = length(binEdges{iDim})-1;
    cc(:,iDim) = round((coords(:,iDim)-binEdges{iDim}(1))/dBin+0.5);
end
% this is a ridiculously fast and not very intutive way to build the maps
% look into the help of sparse to understand how it works
accMap = full(sparse(cc(:,1), cc(:,2), signal, nBins(1), nBins(2)));

return;

% this code is about 100 times slower... (even after major optimization effort)

[nSamples, nDims] = size(coords);
indices = cell(nDims, 1);
nBins = nan(1, nDims);
% first, lets build indices matrices
for iDim = 1:nDims
    nBins(iDim) = length(binEdges{iDim})-1;
    indices{iDim} = false(nSamples, nBins(iDim));
    for iBin = 1:nBins(iDim)
        indices{iDim}(:, iBin) = ...
            coords(:, iDim)>=binEdges{iDim}(iBin) & ...
            coords(:, iDim)<binEdges{iDim}(iBin+1);
    end
end
accMap = nan(nBins);

% trying to be independent of the nDims (not to hard-code the
% dimensionality)
str = '[sub(1)';
for iDim = 2:nDims
    str = sprintf('%s, sub(%d)', str, iDim);
end
str = [str, '] = ind2sub(nBins, iElement);'];

nElements = prod(nBins);
for iElement = 1:nElements
    if nDims == 1
        accMap(iElement) = sum(signal(indices{1}(:, iElement)));
    elseif nDims ==2
        [sub(1), sub(2)] = ind2sub(nBins, iElement);
        accMap(iElement) = sum(signal(indices{1}(:, sub(1)) & indices{2}(:, sub(2))));
    else
        eval(str);
        ind = indices{1}(:, sub(1));
        for iDim=2:length(sub);
            ind = ind & indices{iDim}(:, sub(iDim));
        end
        accMap(iElement) = sum(signal(ind));
    end
end

end % buildAccumMap()