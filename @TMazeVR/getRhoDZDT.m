function rhoDZDT = getRhoDZDT(obj, res)

rhoDZDT = [];
for iPlane = obj.Planes
    fData = obj.data2p{iPlane}.F;
    tData = obj.times2p{iPlane};
    f0 = prctile(fData, 20);
    fData = bsxfun(@rdivide, bsxfun(@minus, fData, f0), f0);
    
    nCells = size(fData, 2);
    nTrials = obj.dataTMaze.nTrials;
    meanValuesModel = nan(nTrials, nCells);
    meanValuesData = nan(nTrials, nCells);
    for iTrial = 1:nTrials
        [dThVector, dZVector, fVector, ~] = buildDVectors(obj, iTrial, tData, fData);
        if ~isempty(dThVector)
            for iCell = 1:nCells
                dZAxis = res{iPlane}(iCell).dZdThetaBinCentres{1};
                dThAxis = res{iPlane}(iCell).dZdThetaBinCentres{2};
                theMap = res{iPlane}(iCell).dZdThetaMap;
                fModel = interp2(dThAxis, dZAxis, theMap, dThVector, dZVector, 'spline');%, 0);
                meanValuesModel(iTrial, iCell) = nanmean(fModel);
            end
        end
        meanValuesData(iTrial, :) = nanmean(fVector);
    end
    rho = nan(nCells, 1);
    for iCell = 1:nCells
        tmp = corrcoef(meanValuesModel(:, iCell), meanValuesData(:, iCell), 'rows', 'complete');
        rho(iCell) = tmp(2);
    end
    rhoDZDT = cat(1, rhoDZDT, rho);
end