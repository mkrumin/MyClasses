function [mapOut, binsOut] = calcSingleMap(obj, iPlane, iROI, options)

% vSfN of TrainMaps - based on v6, data and modeled traces overlayed, only
% scaled model presented in the traces

if nargin<4
    options = struct;
end
options = fillOptions(obj, options);

% nChunks = options.cvFactor;

fData = obj.data2p{iPlane}.F(:, iROI);
tData = obj.times2p{iPlane};

trialIdx = 1:length(obj.dataTMaze.contrastSequence);
validTrials = find(~isnan(obj.timesTmazeOnset));
trialIdx = intersect(trialIdx, validTrials);
[thAll, zAll, fAll] = buildVectors(obj, trialIdx, tData, fData);

zMax = max(abs(zAll));
zMin = min(abs(zAll));
thMax = max(abs(thAll));
% clear zAll thAll;

% now, a complicated way to define bin edges
thEdges = buildEdges(thMax, options.dTheta);
zEdges = buildEdges([zMin, zMax], options.dZ);


hGauss = ndGaussian(options.optStd);
[theMap, binCentres] = estMap([zAll, thAll], fAll, {zEdges, thEdges}, hGauss);

mapOut = theMap;
binsOut = binCentres;

end
% return;
% 
% [thMatrix, zMatrix, fMatrix] = buildTraces(obj, 1:length(obj.dataTMaze.contrastSequence), tData, fData);

%% ==========================================================================

function optOut = fillOptions(obj, optIn)

optOut = optIn;
if ~isfield(optIn, 'contrastGroups')
    contrasts = unique(obj.dataTMaze.contrastSequence);
    optOut.contrastGroups = mat2cell(contrasts(:), ones(size(contrasts)), 1);
    optOut.contrastGroupLabels = cell(size(optOut.contrastGroups));
    for iLabel = 1:length(optOut.contrastGroups)
        optOut.contrastGroupLabels{iLabel} = [num2str(contrasts(iLabel)), '%'];
    end
    optOut.contrastGroups{end+1} = contrasts;
    optOut.contrastGroupLabels{end+1} = 'All';
end

if ~isfield(optIn, 'decisionGroups')
    optOut.decisionGroups = {'L'; 'R'; 'LR'};
    optOut.decisionGroupLabels = {'Go L'; 'Go R'; 'All'};
end

if ~isfield(optIn, 'dZ')
    optOut.dZ = 3; % [cm]
end

if ~isfield(optIn, 'dTheta')
    optOut.dTheta = 2; % [deg]
end

if ~isfield(optIn, 'cvFactor')
    optOut.cvFactor = 5; %
end

if ~isfield(optIn, 'cvMethod')
    optOut.cvMethod = 'random'; % 'random' or 'block'
end

end %fillOptions();

