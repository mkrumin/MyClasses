function [map, binCenters, stats] = trajDensityMap(obj, opt)

% training the decision (Prob(goR)) map as a function of z-theta space
% It is based on the trajectories the animal took during the behavior
% map = trainDecisionMap(obj)
%      obj: TMaze object
%      map: a 2D map of the Prob(R)
%      binCenters: the coordinates of the z-theta dimensions


nContrasts = length(opt.cGroups);
map = cell(nContrasts, 2);
[zMat, thMat, dMat] = getTraces(obj);
filterStd = opt.filterStd;
stats = struct('xx', [], 'yy', []);
for iC = 1:nContrasts
    % trials belonging to a contrast group
    trialIdx = ismember(obj.contrastSequence, opt.cGroups{iC});
    % and also were finished
    trialIdxL = find(trialIdx & ismember(obj.report, 'L')');
    trialIdxR = find(trialIdx & ismember(obj.report, 'R')');
    zL = zMat(:, trialIdxL);
    thL = thMat(:, trialIdxL);
%     dL = dMat(:, trialIdxL);
    zR = zMat(:, trialIdxR);
    thR = thMat(:, trialIdxR);
%     dR = dMat(:, trialIdxR);
    
    binEdges = {[4.5:105.5]', [-30.5:30.5]'};
    [occMapL, binCenters] = buildOccupMap([zL(:), thL(:)], binEdges);
    map{iC, 1} = filterMap(occMapL, ndGaussian(filterStd));
    [occMapR, ~] = buildOccupMap([zR(:), thR(:)], binEdges);
    map{iC, 2} = filterMap(occMapR, ndGaussian(filterStd));
    tmp = trajStats(map{iC, 1}, binCenters, opt.alpha);
    stats(iC, 1).xx = tmp.xx;
    stats(iC, 1).yy = tmp.yy;
    tmp = trajStats(map{iC, 2}, binCenters, opt.alpha);
    stats(iC, 2).xx = tmp.xx;
    stats(iC, 2).yy = tmp.yy;
    d = stats(iC, 2).xx(:, 1)-stats(iC, 1).xx(:, 3);
    zz = stats(iC, 2).yy;
    pp = csaps(zz, d, 1);
    zCross = fnzeros(pp);
    if ~isempty(zCross)
    zCross = zCross(1);
    stats(iC, 1).thCross = interp1(zz, stats(iC, 1).xx(:, 3), zCross, 'linear');
    stats(iC, 2).thCross = interp1(zz, stats(iC, 2).xx(:, 1), zCross, 'linear');
    stats(iC, 1).zCross = zCross;
    stats(iC, 2).zCross = zCross;
    else
    stats(iC, 1).thCross = NaN;
    stats(iC, 2).thCross = NaN;
    stats(iC, 1).zCross = NaN;
    stats(iC, 2).zCross = NaN;
    end

end

if ~opt.doPlotting
    return;
end
%% plotting

if iscell(obj.ExpRef)
    nameStr = cell2mat(obj.ExpRef)';
else
    nameStr = obj.ExpRef;
end
figure('Name', nameStr)

trajMaps = map;
nMaps = length(trajMaps(:));
for iMap = 1:nMaps
%     trajMaps{iMap} = bsxfun(@rdivide, map, max(map, [], 2));
    trajMaps{iMap} = bsxfun(@rdivide, trajMaps{iMap}, sum(trajMaps{iMap}, 2));
end

nMaps = size(trajMaps, 1);
for iMap = 1:nMaps
    mm = max(trajMaps{iMap, 1}(:)+trajMaps{iMap, 2}(:));

    st = stats(iMap, :);
    
    subplot(1, nMaps, iMap);
    im = zeros([size(trajMaps{1}), 3]);
    im(:, :, 1) = 1-trajMaps{iMap, 2}/mm; % red channel
    im(:, :, 2) = 1-trajMaps{iMap, 2}/mm-trajMaps{iMap, 1}/mm; % green channel
    im(:, :, 3) = 1-trajMaps{iMap, 1}/mm; % blue channel
    image(binCenters{2}, binCenters{1}, im);
    axis xy equal tight
    
    hold on;
    plot(st(1).xx(:, 2), st(1).yy, 'r', 'LineWidth', 1);
    plot(st(1).xx(:, [1 3]), st(1).yy, 'r:', 'LineWidth', 1);
    plot(st(2).xx(:, 2), st(2).yy, 'b', 'LineWidth', 1);
    plot(st(2).xx(:, [1 3]), st(2).yy, 'b:', 'LineWidth', 1);
    % colorbar
    
    zCross = st(1).zCross;
    thCross = (st(1).thCross + st(2).thCross)/2;
    plot(xlim, zCross*[1 1], 'k:')
    plot(thCross, zCross, 'ko', 'MarkerSize', 5, 'LineWidth', 2);
    
    title(opt.cLabels{iMap});

end

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

function mapOut = filterMap(mapIn, hGauss)

% make a separable filter - the filtering will be significantly faster
hGauss1 = sum(hGauss, 2);
hGauss2 = sum(hGauss, 1);
    
flatMap = ones(size(mapIn));
normMap = conv2(hGauss1, hGauss2, flatMap, 'same');

% filtering using conv2() is much (10x ?) faster, at least when the filter is separable
mapOut = conv2(hGauss1, hGauss2, mapIn, 'same')./normMap;

