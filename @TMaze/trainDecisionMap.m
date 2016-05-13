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
    
    filterStd = getOptStd(occMap(:,:,trialIdx), accumMap(:,:,trialIdx), opt.cvFactor);
    map{iC} = filterAndDivideMaps(sum(occMap(:,:,trialIdx), 3), sum(accumMap(:,:,trialIdx), 3), filterStd, 0.01);
    res(iC).std = filterStd;
end
reshape(res, size(map));

%==============================================================================
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

%=================================================
function coordsOut = evenInterp(coordsIn, nPoints)

if nargin<2 || isempty(nPoints)
    nPoints = 100;
end

dCoords = diff(coordsIn);
travel = [0; cumsum(sqrt(sum(dCoords.^2, 2))+eps(1000))];
travelOut = linspace(0, travel(end), nPoints);
coordsOut = interp1(travel, coordsIn, travelOut);

%==============================================================================
function [filterStd, cvErr] = getOptStd(occMap, accumMap, cvFactor)

[nZ, nTh, nTrials] = size(occMap);
filterStd = [nan nan];
cvErr = nan;
epsilon = 0.01;
doSpeedUp = true;

dataClass = 'double';
if doSpeedUp
    occMap = single(occMap);
    accumMap = single(accumMap);
    dataClass = 'single';
end

% define the grid of z and th for the search
pow = 3;
zGrid = linspace(0.5^(1/pow), nZ^(1/pow), 15).^pow;
thGrid = linspace(0.5^(1/pow), nTh^(1/pow), 15).^pow;
zGrid(end) = Inf;
thGrid(end) = Inf;
% zGrid = linspace(2^(1/pow), (nZ/3)^(1/pow), 15).^pow;
% thGrid = linspace(2^(1/pow), (nTh/3)^(1/pow), 15).^pow;

% divide trials for cross-validation
if (cvFactor>nTrials)
    cvFactor = nTrials;
end

groupings = nan(cvFactor, ceil(nTrials/cvFactor));
groupings(1:nTrials) = randperm(nTrials);

% extensive grid search for the optimal filter parameters
err = nan(length(zGrid), length(thGrid), cvFactor);
for iGroup = 1:cvFactor
    idxTest = groupings(iGroup, :);
    idxTest = idxTest(~isnan(idxTest));
    occTest = sum(occMap(:,:,idxTest), 3);
    accumTest = sum(accumMap(:,:,idxTest), 3);
    idxTrain = setdiff(1:nTrials, idxTest);
    occTrain = sum(occMap(:,:,idxTrain), 3);
    accumTrain = sum(accumMap(:,:,idxTrain), 3);
    normMap = ones(size(accumTrain), dataClass);
    %     meanSignal = sum(accumTrain(:))/sum(occTrain(:));
    meanSignal = 0.5;
    for iZ = 1:length(zGrid)
        if isinf(zGrid(iZ))
            occTrainF1 = repmat(mean(occTrain, 1), nZ, 1);
            accumTrainF1 = repmat(mean(accumTrain, 1), nZ, 1);
            normMapF1 = repmat(mean(normMap, 1), nZ, 1);
        else
            hGauss1 = ndGaussian(zGrid(iZ));
            occTrainF1 = conv2(hGauss1, 1, occTrain, 'same');
            accumTrainF1 = conv2(hGauss1, 1, accumTrain, 'same');
            normMapF1 = conv2(hGauss1, 1, normMap, 'same');
        end
        for iTh = 1:length(thGrid)
            if isinf(thGrid(iTh))
                occTrainF1 = repmat(mean(occTrain, 2), 1, nTh);
                accumTrainF1 = repmat(mean(accumTrain, 2), 1, nTh);
                normMapF1 = repmat(mean(normMap, 2), 1, nTh);
            else
                hGauss2 = ndGaussian(thGrid(iTh));
                occTrainF = conv2(1, hGauss2, occTrainF1, 'same');
                accumTrainF = conv2(1, hGauss2, accumTrainF1, 'same');
                normMapF = conv2(1, hGauss2, normMapF1, 'same');
            end
            
            map = (accumTrainF./normMapF+epsilon*meanSignal)./...
                (occTrainF./normMapF+epsilon);
            
            %             map = filterAndDivideMaps(occTrain, accumTrain, hGauss, 0.01);
            err(iZ, iTh, iGroup) = mapError(map, occTest, accumTest);
        end
    end
end

overallErr = mean(err, 3);
[cvErr, optIdx] = min(overallErr(:));
[optZ, optTh] = ind2sub([length(zGrid), length(thGrid)], optIdx);
filterStd = [zGrid(optZ), thGrid(optTh)];

%==============================================================================
function err = mapError(trainedMap, occMap, accumMap)

err = 0;
idx = find(occMap);
weights = occMap(idx);
signal = accumMap(idx)./weights;
model = trainedMap(idx);

err = sum((signal-model).^2.*weights)/sum(weights);

