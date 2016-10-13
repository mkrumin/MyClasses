function ballData = getRunningSpeed(info)

filenames = dat.expFilePath(info.expRef, 'Timeline');

try 
    load(filenames{1});
catch
    load(filenames{2});
end

nEvents = Timeline.ballUDPCount;
ballData.t = Timeline.ballUDPTimes(1:nEvents);
ballData.timestamps = nan(nEvents, 1);
ballData.total = nan(nEvents, 1);
ballData.forward = nan(nEvents, 1);
ballData.rotation = nan(nEvents, 2);
ballData.sideways = nan(nEvents, 1);
% tic
for iEvent = 1:nEvents
    str = Timeline.ballUDPEvents{iEvent};
    vars = str2num(str);
    % vars == [timestamp, dax, day, dbx, dby]
    ballData.timestamps(iEvent) = vars(1);
    ballData.forward(iEvent) = -vars(4);
    ballData.rotation(iEvent,:) = [vars(3), vars(5)];
    ballData.sideways(iEvent) = vars(2);
%     ballData.total(iEvent) = norm([ballData.forward(iEvent), ballData.rotation(iEvent), ballData.sideways(iEvent)]);
end

nonZeroIdx = ballData.forward ~= 0 | [ballData.forward(2:end); 0] == 0 | ...
    [0; ballData.forward(1:end-1)] == 0;
ballData.forward = interp1(ballData.t(nonZeroIdx), ballData.forward(nonZeroIdx), ...
    ballData.t, 'linear', 'extrap');
nonZeroIdx = ballData.rotation(:,1) ~= 0 | ...
    [ballData.rotation(2:end,1); 0] == 0 | ...
    [0; ballData.rotation(1:end-1,1)] == 0;
ballData.rotation(:,1) = interp1(ballData.t(nonZeroIdx), ...
    ballData.rotation(nonZeroIdx,1), ballData.t, 'linear', 'extrap');
nonZeroIdx = ballData.rotation(:,2) ~= 0 | ...
    [ballData.rotation(2:end,2); 0] == 0 | ...
    [0; ballData.rotation(1:end-1,2)] == 0;
ballData.rotation(:,2) = interp1(ballData.t(nonZeroIdx), ...
    ballData.rotation(nonZeroIdx,2), ballData.t, 'linear', 'extrap');
ballData.rotation = mean(ballData.rotation,2);
nonZeroIdx = ballData.sideways ~= 0 | [ballData.sideways(2:end); 0] == 0 | ...
    [0; ballData.sideways(1:end-1)] == 0;
ballData.sideways = interp1(ballData.t(nonZeroIdx), ...
    ballData.sideways(nonZeroIdx), ballData.t, 'linear', 'extrap');

ballData.total = sqrt(ballData.forward.^2 + ballData.rotation.^2 + ballData.sideways.^2);
% toc

% figure
% plot(ballData.t, filtfilt([1 1 1 1 1 1]/6, 1, total))