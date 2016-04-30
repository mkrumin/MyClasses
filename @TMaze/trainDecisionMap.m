function [map, binCenters, res] = trainDecisionMap(obj, opt)

% training the decision (Prob(goR)) map as a function of z-theta space
% It is based on the trajectories the animal took during the behavior
% map = trainDecisionMap(obj)
%      obj: TMaze object
%      map: a 2D map of the Prob(R)
%      binCenters: the coordinates of the z-theta dimensions


nContrasts = length(opt.cGroups);
map = cell(nContrasts, 1);
res = struct('std', [nan nan]);
[zMat, thMat, dMat] = getTraces(obj);
% filterStd = opt.filterStd;
binEdges = {[4.5:105.5]', [-30.5:30.5]'};

nTrials = size(zMat, 2);
nZ = length(binEdges{1})-1;
nTh = length(binEdges{2})-1;

accumMap = nan(nZ, nTh, nTrials);
occMap = nan(nZ, nTh, nTrials);
for iTrial = 1:nTrials
    if ~ismember(obj.report(iTrial), 'LR')
        %only build maps for finished trials
        continue;
    end
    accumMap(:,:,iTrial) = ...
        buildAccumMap([zMat(:, iTrial), thMat(:, iTrial)], dMat(:, iTrial), binEdges);
    [occMap(:,:,iTrial), binCenters] = ...
        buildOccupMap([zMat(:, iTrial), thMat(:, iTrial)], binEdges);
end

for iC = 1:nContrasts
    % trials belonging to a contrast group
    trialIdx = ismember(obj.contrastSequence, opt.cGroups{iC});
    % and also were finished
    trialIdx = find(trialIdx & ismember(obj.report, 'LR')');
%     z = zMat(:, trialIdx);
%     th = thMat(:, trialIdx);
%     d = dMat(:, trialIdx);
    
    filterStd = getOptStdCV(occMap(:,:,trialIdx), accumMap(:,:,trialIdx), Inf);
    map{iC} = filterAndDivideMaps(sum(occMap(:,:,trialIdx), 3), sum(accumMap(:,:,trialIdx), 3), ndGaussian(filterStd), 0.01);
    res(iC).std = filterStd;
end

%==================================================================================
function [z, th, d] = getTraces(obj, trIdx)

if nargin<2 || isempty(trIdx)
    nTrials = obj.nTrials;
    trIdx = 1:nTrials;
end
nPoints = 200;

z =nan(nPoints, length(trIdx));
th =nan(nPoints, length(trIdx));
d =nan(nPoints, length(trIdx));

[~, zInd] = ismember( 'Z', obj.SESSION.allTrials(1).pospars);
[~, thInd] = ismember( 'theta', obj.SESSION.allTrials(1).pospars);

for trialNumber = 1:length(trIdx)
    iTrial = trIdx(trialNumber);
    if ~ismember(obj.report(iTrial), 'LR')
        continue;
    end
    idx = find(obj.SESSION.allTrials(iTrial).trialActive);
    tmpZ = obj.SESSION.allTrials(iTrial).posdata(idx, zInd);
    tmpTh = obj.SESSION.allTrials(iTrial).posdata(idx, thInd);
    coords = evenInterp([tmpZ, tmpTh], nPoints);
    tmpZ = coords(:, 1);
    tmpTh = coords(:, 2);
    tmpD = ones(size(tmpZ))*(obj.report(iTrial) == 'R');
    z(:, trialNumber) = tmpZ;
    th(:, trialNumber) = tmpTh;
    d(:, trialNumber) = tmpD;
end

z = -z;
th = th*180/pi;


function coordsOut = evenInterp(coordsIn, nPoints)

if nargin<2 || isempty(nPoints)
    nPoints = 100;
end

dCoords = diff(coordsIn);
travel = [0; cumsum(sqrt(sum(dCoords.^2, 2))+eps(1000))];
travelOut = linspace(0, travel(end), nPoints);
coordsOut = interp1(travel, coordsIn, travelOut);

function filterStd = getOptStdCV(occMap, accumMap, cvFactor)

nTrials = size(occMap, 3);
filterStd = [1 1];


