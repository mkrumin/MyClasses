function mapOut = filterAndDivideMaps(occMap, sigAccumMap, hGauss, epsilon)

if nargin<4 || isempty(epsilon)
    epsilon = 0.01;
end

meanSignal = sum(sigAccumMap(:))/sum(occMap(:));
flatMap = ones(size(sigAccumMap));

if numel(hGauss)>2
    % this means that it is an actual filter, and not only the std values
% make a separable filter - the filtering will be significantly faster
hGauss1 = sum(hGauss, 2);
hGauss2 = sum(hGauss, 1);
    
%     mapOut = (imfilter(imfilter(sigAccumMap, hGauss1, 0), hGauss2, 0)+epsilon*meanSignal)./...
%         (imfilter(imfilter(occMap, hGauss1, 0), hGauss2, 0)+epsilon);

normMap = conv2(hGauss1, hGauss2, flatMap, 'same');

% filtering using conv2() is much (10x ?) faster, at least when the filter is separable
    mapOut = (conv2(hGauss1, hGauss2, sigAccumMap, 'same')./normMap+epsilon*meanSignal)./...
        (conv2(hGauss1, hGauss2, occMap, 'same')./normMap+epsilon);
else
    if isinf(hGauss(1))
        sigAccumMap = repmat(mean(sigAccumMap, 1), size(sigAccumMap, 1), 1);
        occMap = repmat(mean(occMap, 1), size(occMap, 1), 1);
        normMap = repmat(mean(flatMap, 1), size(flatMap, 1), 1);
    else
        sigAccumMap = conv2(ndGaussian(hGauss(1)), 1, sigAccumMap, 'same');
        occMap = conv2(ndGaussian(hGauss(1)), 1, occMap, 'same');
        normMap = conv2(ndGaussian(hGauss(1)), 1, flatMap, 'same');
    end
    if isinf(hGauss(2))
        sigAccumMap = repmat(mean(sigAccumMap, 2), 1, size(sigAccumMap, 2));
        occMap = repmat(mean(occMap, 2), 1, size(occMap, 2));
        normMap = repmat(mean(normMap, 2), 1, size(normMap, 2));
    else
        sigAccumMap = conv2(1, ndGaussian(hGauss(2)), sigAccumMap, 'same');
        occMap = conv2(1, ndGaussian(hGauss(2)), occMap, 'same');
        normMap = conv2(1, ndGaussian(hGauss(2)), normMap, 'same');
    end
    
    mapOut = (sigAccumMap./normMap+epsilon*meanSignal)./...
        (occMap./normMap+epsilon);
    
end