function obj = discrClassAnalysis(obj, options)

% this is a method of TMazeVR class
% it performs discriminant classification analysis
Planes = obj.Planes;
nTrials = obj.dataTMaze.nTrials;


zLocations = 10:20:90;
tol = 1;
fDur = 0.5; % sec

thVector = cell(nTrials, max(Planes));
zVector = cell(nTrials, max(Planes));
fVector = cell(nTrials, max(Planes));
tVector = cell(nTrials, max(Planes));

validTrials = [1:nTrials]';
for iPlane = Planes
    for iTrial = 1:nTrials
        tData = obj.times2p{iPlane};
        fData = obj.data2p{iPlane}.F;
        
        [thVector{iTrial, iPlane}, zVector{iTrial, iPlane}, fVector{iTrial, iPlane}, tVector{iTrial, iPlane}] =...
            buildVectors(obj, iTrial, tData, fData);
        if isempty(thVector{iTrial, iPlane})
            validTrials = setdiff(validTrials, iTrial);
        end
    end
end

% exclude timed-out trials
validTrials = setdiff(validTrials, find(obj.dataTMaze.report == 'T'));
nValidTrials = length(validTrials);

theta = nan(nValidTrials, max(Planes), length(zLocations));
F = cell(nValidTrials, max(Planes), length(zLocations));
for iZ = 1:length(zLocations)
    zz = zLocations(iZ);
    for iPlane = Planes
        for trialInd = 1:nValidTrials
            iTrial = validTrials(trialInd);
            dt = mean(diff(tVector{iTrial, iPlane}));
            nSamples = round(fDur/dt);
            idx = find(zVector{iTrial, iPlane}>zz-tol & zVector{iTrial, iPlane}<zz+tol);
            fIdx = [(max(1, idx(1)-nSamples):idx(1)-1)'; idx];
            theta(trialInd, iPlane, iZ) = mean(thVector{iTrial, iPlane}(idx));
            F{trialInd, iPlane, iZ} = mean(fVector{iTrial, iPlane}(fIdx, :));
        end
    end
end

theta = theta(:, min(Planes), :);
F = cell2mat(F);
R = obj.dataTMaze.report(validTrials)';
C = obj.dataTMaze.contrastSequence(validTrials);



% get aligned vectors of z, theta, F, FR for all the finished trials
% sample theta and f at specific z locations
% maybe use integral of fover the last dt period of time
%
% construct classifiers (d'(theta) and d'(f)) for each of the selected z values