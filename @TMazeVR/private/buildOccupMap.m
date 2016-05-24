function [occMap, binCentres] = buildOccupMap(coords, binEdges)

nDims = length(binEdges);
switch nDims
    case 2
        [occMap, binCentres] = hist3(coords, 'Edges', binEdges);
        occMap = occMap(1:end-1, 1:end-1);
        binCentres{1} = binCentres{1}(1:end-1);
        binCentres{2} = binCentres{2}(1:end-1);
    case 1
        [occMap, ~] = histcounts(coords, binEdges{1});
        binCentres{1} = (binEdges{1}(1:end-1) + binEdges{1}(2:end))/2;
    otherwise
        error('buildOccupMap() only supports 1-D and 2-D maps');
        
end

end % buildOccupMap();

