function out = decodeZT(obj, resStruct, opt)

out = [];
doPlot = false;
res = resStruct.res;
nCells = 0;
for iPlane = obj.Planes
    nROIs = obj.nROIs(iPlane);
    nCells = nCells + nROIs;
end

zAxis = res{1}(1).zThetaBinCentres{1};
thAxis = res{1}(1).zThetaBinCentres{2};

%%
nTrials = obj.dataTMaze.nTrials;
nPlanes = max(obj.Planes);
ztMaps = nan(length(zAxis), length(thAxis), nCells);
iCell = 0;
errVals = [];
thVector = cell(nTrials, 1);
zVector = cell(nTrials, 1);
fVector = cell(nTrials, nPlanes);
tVector = cell(nTrials, 1);
for iPlane = obj.Planes
    nROIs = obj.nROIs(iPlane);
    F = obj.data2p{iPlane}.F;
    F0 = prctile(obj.data2p{iPlane}.F, 20);
    F = bsxfun(@rdivide, bsxfun(@minus, F, F0), F0);
    for iTrial = 1:nTrials
        [thVector{iTrial}, zVector{iTrial}, fVector{iTrial, iPlane}, tVector{iTrial}] = ...
            buildVectors(obj, iTrial, obj.times2p{iPlane}, F);
    end
    errVals = cat(2, errVals, [res{iPlane}.errVals]);
    for iROI = 1:nROIs
        iCell = iCell + 1;
        ztMaps(:, :, iCell) = res{iPlane}(iROI).zThetaMap;
        optStd(:, iCell) = res{iPlane}(iROI).optStd';
    end
end

errVals = reshape(errVals, 1, 1, []);
% fVector = reshape(fVector, 1, []);
allF = cell(nTrials, 1);
validTrials = true(nTrials, 1);
for iTrial = 1:nTrials
    allF{iTrial} = cell2mat(fVector(iTrial, :));
    if isempty(allF{iTrial})
        validTrials(iTrial) = false;
    end
end


%%
cells2Use = find(errVals<0.8);
for iTrial = 1:nTrials
    if ~validTrials(iTrial)
        continue;
    end
    trialIdx = setdiff(find(validTrials), iTrial);
    trainF = cell2mat(allF(trialIdx));
    trainF = trainF(:, cells2Use);
    sigmaF = std(trainF);
    sigmaF = reshape(sigmaF, 1, 1, []);
    testF = permute(allF{iTrial}(:, cells2Use), [1 3 2]);
    
    trainZ = cell2mat(zVector(trialIdx));
    trainTh = cell2mat(thVector(trialIdx));
    testZ = zVector{iTrial};
    testTh = thVector{iTrial};
    
    ztMaps = nan(length(zAxis), length(thAxis), length(cells2Use));
    for iCell = 1:length(cells2Use)
        hFilter = ndGaussian(optStd(:, cells2Use(iCell)));
        ztMaps(:, :, iCell) = ...
            estMap([trainZ, trainTh], trainF(:, iCell), {res{1}(1).zEdges, res{1}(1).thEdges}, hFilter);
    end
    
    flatMaps = reshape(ztMaps, [], length(cells2Use));
    flatMaps = permute(flatMaps, [3, 1, 2]);
    
    tmpL = bsxfun(@minus, flatMaps, testF);
    tmpL = tmpL.^2;
    tmpL = bsxfun(@rdivide, tmpL, 2*sigmaF.^2);
    % tmpL = bsxfun(@plus, tmpL, 1-errVals(cells2Use));
    logL = -1*sum(tmpL, 3);
    clear tmpL;
    logL = logL';
    logL = reshape(logL, length(zAxis), length(thAxis), []);
    
    occMap = buildOccupMap([trainZ, trainTh], {res{1}(1).zEdges, res{1}(1).thEdges});
    occMap = occMap/sum(occMap(:));
    
    out(iTrial).logL = logL;
    out(iTrial).occMap = occMap;
    out(iTrial).traj = [testZ, testTh];
    out(iTrial).axes = {zAxis, thAxis};
    
    
    if doPlot
        post = imgaussfilt3(logL, [1 1 3]);
        logOccMap = log(imgaussfilt((occMap), 0.001));
        % logOccMap = (imgaussfilt((log(occMap)), 3));
        % logOccMap = log(occMap);
        
        post = bsxfun(@plus, post, logOccMap);
        for i = 1:3:size(post, 3)
            frameTime = tic;
            LMap = post(:,:,i);
            imagesc(thAxis, zAxis, LMap);
            hold on;
            plot(testTh(i), testZ(i), 'o');
            %     [lMax, loc] = max(LMap(:));
            %     lMin = min(LMap(isfinite(LMap)));
            %     [zSub, thSub] = ind2sub(size(LMap), loc);
            %     plot(thAxis(thSub), zAxis(zSub), 'ro', 'MarkerSize', 20);
            %     cLevel = lMax-0.05*(lMax-lMin);
            %     contour(thAxis, zAxis, LMap, [cLevel cLevel], 'c');
            hold off;
            axis xy equal tight;
            title(sprintf('Trial %d/%d, %d%%', iTrial, nTrials, obj.dataTMaze.contrastSequence(iTrial)))
            colormap hot
            
            while toc(frameTime)<1/60
                ;
            end
            drawnow;
        end
    end
end


end

%%

