function plotTrajFrame(obj, options)

    nGroups = length(options);
    if isfield(options(1), 'nSamples')
        nSamples = options(1).nSamples;
    else
        nSamples = size(obj.posUniform.z, 1);
    end
    
    for iGroup = 1:nGroups
        stats(iGroup) = getStats(obj.posUniform, options(iGroup).idx);
    end
    % first plot the patches (the std areas)
    cla;
    for iGroup = 1:nGroups
%         xPatch = [stats(iGroup).thQuartiles(1:nSamples, 1); stats(iGroup).thMedian(nSamples); flipud(stats(iGroup).thQuartiles(1:nSamples, 2))];
%         yPatch = [stats(iGroup).zQuartiles(1:nSamples, 1); stats(iGroup).zMedian(nSamples); flipud(stats(iGroup).zQuartiles(1:nSamples, 2))];
        xPatch = [stats(iGroup).thQuartiles(1:nSamples, 1); flipud(stats(iGroup).thQuartiles(1:nSamples, 2))];
        yPatch = [stats(iGroup).zMedian(1:nSamples); flipud(stats(iGroup).zMedian(1:nSamples))];
        cPatch = zeros(size(xPatch));
        patch(xPatch, yPatch, cPatch, 'EdgeColor', 'none', 'FaceColor', options(iGroup).FaceColor, 'FaceAlpha', options(iGroup).FaceAlpha);
        hold on;
    end
    
    % then plot the trajectories
    for iGroup = 1:nGroups
        plot(obj.posUniform.th(1:nSamples, options(iGroup).idx), obj.posUniform.z(1:nSamples, options(iGroup).idx), ':', 'Color', options(iGroup).LineColor);
    end
    
    % then plot the medians
    for iGroup = 1:nGroups
        plot(stats(iGroup).thMedian(1:nSamples), stats(iGroup).zMedian(1:nSamples), 'Color', options(iGroup).LineColor, 'LineWidth', 2);
    end
    
    % then plot the ends of the trajectories
    for iGroup = 1:nGroups
        plot(obj.posUniform.th(nSamples, options(iGroup).idx), obj.posUniform.z(nSamples, options(iGroup).idx), '.', 'Color', options(iGroup).DotColor);
    end

end % plotTrajFrame()
