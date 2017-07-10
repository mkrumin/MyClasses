function rhoRV = getRhoRV(obj, res)

rhoRV = [];
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
        [rVector, vVector, fVector, ~] = buildMotorVectors(obj, iTrial, tData, fData);
        if ~isempty(rVector)
            for iCell = 1:nCells
                vAxis = res{iPlane}(iCell).rvBinCentres{1};
                rAxis = res{iPlane}(iCell).rvBinCentres{2};
                theMap = res{iPlane}(iCell).rvMap;
                fModel = interp2(rAxis, vAxis, theMap, rVector, vVector, 'spline');%, 0);
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
    rhoRV = cat(1, rhoRV, rho);
end