function [theMap, binCentres] = estMap(coords, signal, binEdges, hSigma, epsilon)

[occMap, binCentres] = buildOccupMap(coords, binEdges);
% occMap = imfilter(occMap, hFilter, 0);
sigAccumMap = buildAccumMap(coords, signal, binEdges);
% sigAccumMap = imfilter(sigAccumMap, hFilter, 0);

hGauss{1} = ndGaussian(hSigma(1));
hGauss{2} = ndGaussian(hSigma(2))';

if nargin<5
    theMap = filterAndDivideMaps(occMap, sigAccumMap, hGauss);
else
    theMap = filterAndDivideMaps(occMap, sigAccumMap, hGauss, epsilon);
end    

end % estMap()
