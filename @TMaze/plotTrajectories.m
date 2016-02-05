function figHandle = plotTrajectories(obj, options)

if nargin<2 || ~isfield(options, 'contrasts')
    options.contrasts = unique(abs(obj.contrastSequence));
end

nContrasts = length(options.contrasts);

% [~, zInd] = intersect(obj.SESSION.allTrials(1).pospars, 'Z');
% [~, thInd] = intersect(obj.SESSION.allTrials(1).pospars, 'theta');
% zMax = obj.EXP.roomLength;
% thMax = obj.EXP.restrictionAngle;
zMax = max(obj.posUniform.z(:));
thMax = max(abs(obj.posUniform.th(:)));

plotOptions(1).FaceColor = [0 0 1];
plotOptions(1).FaceAlpha = 0.4;
plotOptions(1).LineColor = [0 0 1];
plotOptions(1).DotColor = [0 0 1];

plotOptions(2).FaceColor = [1 0 0];
plotOptions(2).FaceAlpha = 0.4;
plotOptions(2).LineColor = [1 0 0];
plotOptions(2).DotColor = [1 0 0];

figHandle = figure('Name', 'Trajectories', 'Position', [50 50 1750 1000]);

for iContrast = 1:nContrasts
    % subdividing the trials
    stimLeftIdx = (obj.contrastSequence == -options.contrasts(iContrast));
    stimRightIdx = (obj.contrastSequence == options.contrasts(iContrast));
    wentLeftIdx = (stimRightIdx | stimLeftIdx) & (obj.report == 'L')';
    wentRightIdx = (stimRightIdx | stimLeftIdx) & (obj.report == 'R')';
    
    % plotting 'went right' vs 'went left'
    subplot(2, nContrasts, iContrast);
    plotOptions(1).idx = wentLeftIdx;
    plotOptions(2).idx = wentRightIdx;
    
    plotTrajFrame(obj, plotOptions);
    xlim([-thMax, thMax]+[-0.1 0.1]);
    ylim([0 zMax+1]);
    axis off;
    
    % plotting 'stimulus right' vs 'stimulus left'
    % Unnecessary to plot twice for 0% stimulus
    if (options.contrasts(iContrast) ~= 0)
        subplot(2, nContrasts, nContrasts+iContrast);
        plotOptions(1).idx = stimLeftIdx;
        plotOptions(2).idx = stimRightIdx;
    
        plotTrajFrame(obj, plotOptions);
        xlim([-thMax, thMax]+[-0.1 0.1]);
        ylim([0 zMax+1]);
    end
    axis off;
    
end

drawnow;

end % plotTrajectories()

