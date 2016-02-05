function mapOut = filterAndDivideMaps(occMap, sigAccumMap, hGauss, epsilon)

if nargin<4 || isempty(epsilon)
    epsilon = 0.01;
end

meanSignal = sum(sigAccumMap(:))/sum(occMap(:));

% make a separable filter - the filtering will be significantly faster
hGauss1 = sum(hGauss, 2);
hGauss2 = sum(hGauss, 1);
    
%     mapOut = (imfilter(imfilter(sigAccumMap, hGauss1, 0), hGauss2, 0)+epsilon*meanSignal)./...
%         (imfilter(imfilter(occMap, hGauss1, 0), hGauss2, 0)+epsilon);

% filtering using conv2() is much (10x ?) faster, at least when the filter is separable
    mapOut = (conv2(hGauss1, hGauss2, sigAccumMap, 'same')+epsilon*meanSignal)./...
        (conv2(hGauss1, hGauss2, occMap, 'same')+epsilon);

% this is the older real 2D filtering - it was even hacked to make it faster
% mapOut = (myimfilter(sigAccumMap, hGauss, 0)+epsilon*meanSignal)./...
%     (myimfilter(occMap, hGauss, 0)+epsilon);
