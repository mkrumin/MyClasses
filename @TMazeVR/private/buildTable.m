function tbl = buildTable(obj)

%%
tt = [];
for iPlane = obj.Planes
    tt = cat(1, tt, obj.times2p{iPlane}');
end
tt = sort(tt, 'ascend');

F = [];
for iPlane = obj.Planes
    F = cat(2, F, interp1(obj.times2p{iPlane}', obj.data2p{iPlane}.F, tt));
end


pupil_x = obj.dataEye.results.x;
pupil_y = obj.dataEye.results.y;
pupil_area = obj.dataEye.results.area;

ttEye = 1:length(pupil_x);

% interpolating samples of bad fits
nanIdx = isnan(pupil_x);
pupil_x(nanIdx) = interp1(ttEye(~nanIdx), pupil_x(~nanIdx), ttEye(nanIdx));
nanIdx = isnan(pupil_y);
pupil_y(nanIdx) = interp1(ttEye(~nanIdx), pupil_y(~nanIdx), ttEye(nanIdx));
nanIdx = isnan(pupil_area);
pupil_area(nanIdx) = interp1(ttEye(~nanIdx), pupil_area(~nanIdx), ttEye(nanIdx));

pupil_x = smooth(pupil_x);
pupil_y = smooth(pupil_y);
pupil_area = smooth(pupil_area);

blink = obj.dataEye.results.blink;

% substituting during-blink data with NaNs
pupil_x(blink) = nan;
pupil_y(blink) = nan;
pupil_area(blink) = nan;

pupil_x = interp1(obj.timesEye, pupil_x, tt);
pupil_y = interp1(obj.timesEye, pupil_y, tt);
pupil_area = interp1(obj.timesEye, pupil_area, tt);

blink = interp1(obj.timesEye, double(blink), tt, 'nearest', 0);
blink = logical(blink);

p = polyfit([1001:length(obj.dataBall.t)]', obj.dataBall.t(1001:end), 1);
ballTimes = [1:length(obj.dataBall.t)]'*p(1) + p(2);

ball_rotation = interp1(ballTimes, obj.dataBall.rotation, tt, 'linear', 'extrap');
ball_sideways = interp1(ballTimes, obj.dataBall.sideways, tt, 'linear', 'extrap');
ball_forward = interp1(ballTimes, obj.dataBall.forward, tt, 'linear', 'extrap');
ball_rotation = smooth(ball_rotation, 31);
ball_sideways = smooth(ball_sideways, 31);
ball_forward = smooth(ball_forward, 31);

nTrials = obj.dataTMaze.nTrials;

x = nan(size(tt));
z = nan(size(tt));
theta = nan(size(tt));
contrast = nan(size(tt));
decision = nan(size(tt));
rewardedTrial = nan(size(tt));
t_fromTrialStart = nan(size(tt));
t_toTrialEnd = nan(size(tt));

for iTrial = 1:nTrials
    if isempty(obj.timesVRframes(iTrial).t)
        continue;
    end
    tTrial = obj.timesVRframes(iTrial).t(2:end-1);
    idx = obj.timesVRframes(iTrial).idx(2:end);
    tStart = tTrial(1);
    tEnd = tTrial(end);
    
    samples = find(tt>=tStart & tt<=tEnd);
    contrast(samples) = obj.dataTMaze.contrastSequence(iTrial);
    decision(samples) = (obj.dataTMaze.report(iTrial) == 'R') - (obj.dataTMaze.report(iTrial) == 'L');
    rewardedTrial(samples) = (obj.dataTMaze.outcome(iTrial) == 'C');
    
end

%%


%%
figure;

data = [contrast, decision, rewardedTrial, pupil_x, pupil_y, pupil_area, ball_rotation, ball_forward, ball_sideways];
data = bsxfun(@minus, data, nanmean(data));
data = bsxfun(@rdivide, data, nanstd(data));
plot(tt, data)
%%
tbl = table(F, 'VariableNames', {'F'});

