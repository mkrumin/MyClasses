function eta = getZScaledETA(t, signal, zTimes, z, eventWindows, nPoints)

if nargin < 6
    nPoints = 100;
end

zMin = min(z);
zMax = max(z);
zAxis = linspace(zMin, zMax, nPoints);

nEvents = size(eventWindows, 1);

allSnippets = nan(nEvents, nPoints);

for iEvent = 1:nEvents
    if isnan(eventWindows(iEvent, 1))
        continue;
    end
    
    idx = (zTimes >= eventWindows(iEvent, 1) & zTimes <= eventWindows(iEvent, 2));
    zEvent = z(idx);
    sEvent = interp1(t, signal, zTimes(idx));
    [~, ind] = histc(zEvent, zAxis);
    snippet = nan(1, nPoints);
    for iBin = 1:nPoints
        snippet(iBin) = mean(sEvent(ind==iBin));
    end
    notNanIdx = ~isnan(snippet);
    if sum(notNanIdx)>1 % otherwise leave all the NaNs as they are
        snippet = interp1(zAxis(notNanIdx), snippet(notNanIdx), zAxis, 'linear');
    end
    allSnippets(iEvent, :) = snippet;
end

eta.z = zAxis;
eta.mean = nanmean(allSnippets);
eta.sem = nanstd(allSnippets)/sqrt(nEvents);
eta.std = nanstd(allSnippets);
eta.all = allSnippets;



