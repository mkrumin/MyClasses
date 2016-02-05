function obj = getPosUniformT(obj, nSamples)

if nargin<2
    nSamples = 100;
end

trials = obj.SESSION.allTrials(1:obj.nTrials);
nTrials = length(trials);
[~, zInd] = intersect(trials(1).pospars, 'Z');
[~, xInd] = intersect(trials(1).pospars, 'X');
[~, thInd] = intersect(trials(1).pospars, 'theta');
[~, spInd] = intersect(trials(1).pospars, 'speed');

%% ---- time-based and interpolated single trajectories approach ----
% time is uniformly sampled at nSamples points, 
% and we get z(t), x(t), theta(t), and speed(t) vectors of length nSamples

framesUniform = linspace(0, 1, nSamples);
zUniform = nan(length(framesUniform), nTrials);
xUniform = nan(length(framesUniform), nTrials);
thUniform = nan(length(framesUniform), nTrials);
spUniform = nan(length(framesUniform), nTrials);
for iTrial = 1:nTrials
    zStart = trials(iTrial).posdata(1, zInd);
    indStart = find(trials(iTrial).posdata(:, zInd) ~= zStart, 1, 'first') - 1;
    zVector = -trials(iTrial).posdata(indStart:end, zInd);
    xVector = trials(iTrial).posdata(indStart:end, xInd);
    thVector = trials(iTrial).posdata(indStart:end, thInd);
    spVector = trials(iTrial).posdata(indStart:end, spInd);
    
    frameAxis = linspace(0, 1, length(zVector));
    
    zUniform(:, iTrial) = interp1(frameAxis, zVector, framesUniform);
    xUniform(:, iTrial) = interp1(frameAxis, xVector, framesUniform);
    thUniform(:, iTrial) = interp1(frameAxis, thVector, framesUniform);
    spUniform(:, iTrial) = interp1(frameAxis, spVector, framesUniform);
end

obj.posUniform.z = zUniform;
obj.posUniform.x = xUniform;
obj.posUniform.th = thUniform;
obj.posUniform.sp = spUniform;

end % getPosUniformT()
