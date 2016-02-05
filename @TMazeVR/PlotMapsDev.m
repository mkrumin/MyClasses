function PlotMapsDev(obj, iPlane, iROI)

fData = obj.data2p{iPlane}.F(:, iROI);
contrasts = unique(obj.dataTMaze.contrastSequence);
%     contrasts = [50 25 12 6 0 -50 -25 -12 -6];
contrastGroups = {contrasts(contrasts<0);...
    0;...
    contrasts(contrasts>0);
    contrasts};
nContrastGroups = length(contrastGroups);
contrastGroupLabels = {'StimL All';
    '0% contrast';
    'StimR All';
    'All Trials'};

figHandle = gcf;
clf(figHandle, 'reset');
set(figHandle, 'Name', sprintf('ROI # %d processing...', iROI));
drawnow;

trialIdx = 1:length(obj.dataTMaze.contrastSequence);
[thAll, zAll, ~] = buildVectors(obj, trialIdx, obj.times2p{iPlane}, fData);

zMax = max(abs(zAll));
zMin = min(abs(zAll));
thMax = max(abs(thAll));
clear zAll thAll;


dTheta = 3; % [deg]
dZ = 5; % [cm]
thEdges = buildEdges(thMax, dTheta);
zEdges = buildEdges([zMin, zMax], dZ);
hGauss = ndGaussian([1, 1]);

fMap = cell(nContrastGroups, 1);
thVector = cell(nContrastGroups, 1);
zVector = cell(nContrastGroups, 1);
fVector = cell(nContrastGroups, 1);
for iContrastGroup = 1:nContrastGroups
    % subdividing the trials
    trialIdx = find(ismember(obj.dataTMaze.contrastSequence, contrastGroups{iContrastGroup}));
    if isempty(trialIdx)
        continue;
    end
    % extracting part of the data
    [thVector{iContrastGroup}, zVector{iContrastGroup}, fVector{iContrastGroup}] = buildVectors(obj, trialIdx, obj.times2p{iPlane}, fData);
    % estimating the map
    [fMap{iContrastGroup}, binCentres] = estMap([zVector{iContrastGroup}, thVector{iContrastGroup}], fVector{iContrastGroup}, {zEdges, thEdges}, hGauss);
end

% calculating the range of the colormap
minF = Inf;
maxF = -Inf;
for iContrastGroup = 1:nContrastGroups
    minF = min(minF, min(fMap{iContrastGroup}(:)));
    maxF = max(maxF, max(fMap{iContrastGroup}(:)));
end

% plotting
for iContrastGroup = 1:nContrastGroups
    subplot(1, nContrastGroups, iContrastGroup);
    h = plot(thVector{iContrastGroup}, zVector{iContrastGroup}, '.k', 'MarkerSize', 4);
    set(h, 'Tag', 'Trajectories'); % this will be used later in the WindowButtonDownFcn()
    %     xlim([-thMax, thMax]+[-0.1 0.1]);
    %     ylim([zMin-0.1 zMax+0.1]);
    hold on;
    imagesc(binCentres{2}, binCentres{1}, fMap{iContrastGroup}, 'AlphaData', 0.8);
    axis xy tight
    set(gca, 'CLim', [minF maxF]);
    box off;
    title(contrastGroupLabels{iContrastGroup});
end

set(figHandle, 'Name', sprintf('Cell # %d drawing...', iROI));
drawnow;

set(figHandle, 'Name', sprintf('Cell # %d done.', iROI));

set(figHandle, 'WindowButtonDownFcn', @figClickFcn);
data.obj = handle(obj);
data.iTrial = [];
data.iPlane = iPlane;
data.iROI = iROI;
set(figHandle, 'UserData', data);
end % PlotMaps()

%
