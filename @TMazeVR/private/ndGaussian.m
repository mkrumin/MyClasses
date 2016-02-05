function h = ndGaussian(x)

% x - std of the gaussian in each dimension (in pixels)

nDims = length(x);
sz = nan(nDims, 1);
for iDim = 1:nDims
    sz(iDim) = ceil(3*(x(iDim))*2+1);
end

switch nDims
    case 1
%         h = fspecial('gaussian', [sz(1) 1], x(1));
        t = -((sz(1)-1)/2):((sz(1)-1)/2);
%         h = 1/(sqrt(2*pi)*x(1))*exp(-t.^2/(2*x(1).^2));
        % normalization is not perfect due to finite length, so let's save a few flops
        % as there is normalization at the end of the function
        h = exp(-(t.^2)/(2*x(1)^2));
    case 2
        h = fspecial('gaussian', [sz(1) 1], x(1)) * ...
            fspecial('gaussian', [1 sz(2)], x(2));
    otherwise
        fprintf('ndGaussian not implementd for %d dimensions yet\n', nDims);
end

% normalizing the filter, so that filtering does not affect the mean signal
h = h/sum(h(:));