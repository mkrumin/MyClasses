function [theMap, binCentres] = estMap(coords, signal, binEdges, hFilter, epsilon)

[occMap, binCentres] = buildOccupMap(coords, binEdges);
% occMap = imfilter(occMap, hFilter, 0);
sigAccumMap = buildAccumMap(coords, signal, binEdges);
% sigAccumMap = imfilter(sigAccumMap, hFilter, 0);

if nargin<5
    theMap = filterAndDivideMaps(occMap, sigAccumMap, hFilter);
else
    theMap = filterAndDivideMaps(occMap, sigAccumMap, hFilter, epsilon);
end    

end % estMap()
