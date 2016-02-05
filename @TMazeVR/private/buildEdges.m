function binEdges = buildEdges(xMinMax, dx)

if length(xMinMax) == 1
    % building symmetric bins around zero with range of +-xRange
    dx = xMinMax/ceil(xMinMax/dx);
    binCentres = -xMinMax:dx:xMinMax;
    binEdges = [binCentres-dx/2, binCentres(end)+dx/2];
elseif length(xMinMax)==2
    xMinMax = sort(xMinMax);
    dx = diff(xMinMax)/ceil(diff(xMinMax)/dx);
    binCentres = xMinMax(1):dx:xMinMax(2);
    binEdges = [binCentres-dx/2, binCentres(end)+dx/2];
else
    binEdges = [];
    fprintf('buildEdges(): the range for building edges should be either one or two-value long vector\n');
end