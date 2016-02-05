function stats = getStats(pos, idx)
    stats.zMedian = nanmedian(pos.z(:, idx), 2);
    stats.thMedian = nanmedian(pos.th(:, idx), 2);
    nSamples = size(pos.z, 1);
    stats.zQuartiles = nan(nSamples, 2);
    stats.thQuartiles = nan(nSamples, 2);
    for iSample = 1:nSamples
        notNans = idx & (~isnan(pos.z(iSample,:)))';
        stats.zQuartiles(iSample, :) = prctile(pos.z(iSample, notNans), [25, 75], 2);
        stats.thQuartiles(iSample, :) = prctile(pos.th(iSample, notNans), [25, 75], 2);
    end
end % getStats()