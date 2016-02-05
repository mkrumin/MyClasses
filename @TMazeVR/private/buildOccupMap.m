function [occMap, binCentres] = buildOccupMap(coords, binEdges)

[occMap, binCentres] = hist3(coords, 'Edges', binEdges);
occMap = occMap(1:end-1, 1:end-1);
binCentres{1} = binCentres{1}(1:end-1);
binCentres{2} = binCentres{2}(1:end-1);

end % buildOccupMap();

