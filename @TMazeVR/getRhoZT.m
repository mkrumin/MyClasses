function rhoZT = getRhoZT(obj, res)

rhoZT = [];
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
        [thVector, zVector, fVector, tVector] = buildVectors(obj, iTrial, tData, fData);
        if ~isempty(thVector)
            for iCell = 1:nCells
                zAxis = res{iPlane}(iCell).zThetaBinCentres{1};
                thAxis = res{iPlane}(iCell).zThetaBinCentres{2};
%                 thAxis(1) = thAxis(1)-0.001;
%                 thAxis(end) = thAxis(end)+0.001;
                theMap = res{iPlane}(iCell).zThetaMap;
                fModel = interp2(thAxis, zAxis, theMap, thVector, zVector, 'spline');%, 0);
                meanValuesModel(iTrial, iCell) = mean(fModel);
            end
        end
        meanValuesData(iTrial, :) = mean(fVector);
    end
    rho = nan(nCells, 1);
    for iCell = 1:nCells
        tmp = corrcoef(meanValuesModel(:, iCell), meanValuesData(:, iCell), 'rows', 'complete');
        rho(iCell) = tmp(2);
    end
    rhoZT = cat(1, rhoZT, rho);
end