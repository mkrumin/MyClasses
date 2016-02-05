function eta = getTimeScaledETA(t, signal, eventWindows, nPoints)

if nargin < 4
    nPoints = 100;
end

nEvents = size(eventWindows, 1);

allSnippets = nan(nEvents, nPoints);

for iEvent = 1:nEvents
    tt = linspace(eventWindows(iEvent, 1), eventWindows(iEvent, 2), nPoints);
    allSnippets(iEvent, :) = interp1(t, signal, tt, 'linear');
end

eta.t = linspace(0, 1, nPoints);
eta.mean = nanmean(allSnippets);
eta.sem = nanstd(allSnippets)/sqrt(nEvents);
eta.std = nanstd(allSnippets);
eta.all = allSnippets;



