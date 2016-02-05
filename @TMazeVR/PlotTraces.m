function PlotTraces(obj, iPlane, iROI)

fData = obj.data2p{iPlane}.F(:, iROI);
fData = fData - min(fData);
fData = fData/max(fData)*63;
contrasts = unique(obj.dataTMaze.contrastSequence);
%     contrasts = [50 25 12 6 0 -50 -25 -12 -6];
contrastGroups = {contrasts(contrasts<0);...
    contrasts(contrasts<=-25);...
    contrasts(contrasts<0 & contrasts>-25);...
    0;...
    contrasts(contrasts>0 & contrasts<25);...
    contrasts(contrasts>=25);...
    contrasts(contrasts>0);
    contrasts};
nContrastGroups = length(contrastGroups);
contrastGroupLabels = {'StimL All'; 
    'StimL HighC';
    'StimL LowC';
    '0% contrast';
    'StimR LowC';
    'StimR HighC';
    'StimR All';
    'All Trials'};

figHandle = gcf;
clf(figHandle, 'reset');
set(figHandle, 'Name', sprintf('ROI # %d processing...', iROI));
drawnow;

trialIdx = 1:length(obj.dataTMaze.contrastSequence);
[thVector, zVector, ~] = buildVectors(obj, trialIdx, obj.times2p{iPlane}, fData);

zMax = max(abs(zVector));
zMin = min(abs(zVector));
thMax = max(abs(thVector));

% subplot(2, ceil(nContrasts/2), find(isnan(contrasts)));
% scatter(thVector, zVector, sVector, fVector, '.');
% hold on;
% set(gca, 'CLim', [0 63]);
% xlim([-thMax, thMax]+[-0.1 0.1]);
% ylim([0 zMax+1]);
% axis off;
% title(sprintf('all trials'));

for iContrastGroup = 1:nContrastGroups
    % subdividing the trials
    trialIdx = find(ismember(obj.dataTMaze.contrastSequence, contrastGroups{iContrastGroup}));
    if isempty(trialIdx)
        continue;
    end
    [thVector, zVector, fVector] = buildVectors(obj, trialIdx, obj.times2p{iPlane}, fData);
    sVector = ones(size(thVector));
    
    subplot(1, nContrastGroups, iContrastGroup);
    scatter(thVector, zVector, sVector, fVector, '.');
    hold on;
    set(gca, 'CLim', [0 63]);
    xlim([-thMax, thMax]+[-0.1 0.1]);
    ylim([zMin-0.1 zMax+0.1]);
    axis off;
    title(contrastGroupLabels{iContrastGroup});
end


set(figHandle, 'Name', sprintf('Cell # %d drawing...', iROI));
drawnow;

set(figHandle, 'Name', sprintf('Cell # %d done.', iROI));
end % PlotTraces()
