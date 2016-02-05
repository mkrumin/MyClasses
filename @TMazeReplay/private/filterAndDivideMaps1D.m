function mapOut = filterAndDivideMaps1D(occMap, sigAccumMap, meanSignal, hGauss1, epsilon)

if nargin<5 || isempty(epsilon)
    epsilon = 0.01;
end

% meanSignal = sum(sigAccumMap(:))/sum(occMap(:));

% hGauss1 = hGauss{1};
% hGauss2 = hGauss{2};
    
%     mapOut = (imfilter(imfilter(sigAccumMap, hGauss1, 0), hGauss2, 0)+epsilon*meanSignal)./...
%         (imfilter(imfilter(occMap, hGauss1, 0), hGauss2, 0)+epsilon);

% filtering using conv2() is much (10x ?) faster, at least when the filter is separable
    mapOut = (conv2(hGauss1, 1, sigAccumMap, 'same')+epsilon*meanSignal)./...
        (conv2(hGauss1, 1, occMap, 'same')+epsilon);

% this is the older real 2D filtering - it was even hacked to make it faster
% mapOut = (myimfilter(sigAccumMap, hGauss, 0)+epsilon*meanSignal)./...
%     (myimfilter(occMap, hGauss, 0)+epsilon);
