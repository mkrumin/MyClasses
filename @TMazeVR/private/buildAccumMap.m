function accMap = buildAccumMap(coords, signal, binEdges)

nDims = length(binEdges);
for iDim=1:nDims
    dBin = mean(diff(binEdges{iDim}));
    nBins(iDim) = length(binEdges{iDim})-1;
    cc(:,iDim) = round((coords(:,iDim)-binEdges{iDim}(1))/dBin+0.5);
end
% this is a ridiculously fast and not very intutive way to build the maps
% look into the help of sparse to understand how it works
switch nDims
    case 2
        accMap = full(sparse(cc(:,1), cc(:,2), signal, nBins(1), nBins(2)));
    case 1
        accMap = full(sparse(cc, ones(size(cc)), signal, nBins(1), 1));
    otherwise
        error('Only 1D and 2D maps supported for now');
end

return;
