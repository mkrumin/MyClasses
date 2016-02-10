function resOut = discrClassAnalysis(obj, options)

% this is a method of TMazeVR class
% it performs discriminant classification analysis

if nargin<2
    options = struct;
end

if ~isfield(options, 'zz');
    options.zz = 5:15:95; % cm
end
if ~isfield(options, 'fWin')
    options.fWin = 1; %sec
end
if ~isfield(options, 'thCGroups');
    options.thCGroups = {0; [-6 6 -12 12]; [-25 25 -50 50]};
end
if ~isfield(options, 'fCGroups');
    options.fCGroups = {[-50 -25]; [-12 -6]; 0; [6 12]; [25 50]};
end
if ~isfield(options, 'useFR');
    options.useFR = false;
end

Planes = obj.Planes;
nTrials = obj.dataTMaze.nTrials;


zLocations = options.zz;
tol = 0.5; % cm
fDur = options.fWin; % sec

thVector = cell(nTrials, max(Planes));
zVector = cell(nTrials, max(Planes));
fVector = cell(nTrials, max(Planes));
tVector = cell(nTrials, max(Planes));

validTrials = [1:nTrials]';
for iPlane = Planes
    for iTrial = 1:nTrials
        tData = obj.times2p{iPlane};
        if options.useFR
            fData = obj.data2p{iPlane}.FR;
        else
            fData = obj.data2p{iPlane}.F;
        end
        
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
intF = cell(nValidTrials, max(Planes), length(zLocations));
for iZ = 1:length(zLocations)
    zz = zLocations(iZ);
    for iPlane = Planes
        for trialInd = 1:nValidTrials
            iTrial = validTrials(trialInd);
            dt = mean(diff(tVector{iTrial, iPlane}));
            nSamples = round(fDur/2/dt);
            N = 1;
            idx = find(zVector{iTrial, iPlane}>zz-N*tol & zVector{iTrial, iPlane}<zz+N*tol);
            while isempty(idx)
                N = N*2;
                idx = find(zVector{iTrial, iPlane}>zz-N*tol & zVector{iTrial, iPlane}<zz+N*tol);
            end
            fIdx = [(max(1, idx(1)-nSamples)):min(size(fVector{iTrial, iPlane}, 1), idx(end)+nSamples)]';
            theta(trialInd, iPlane, iZ) = mean(thVector{iTrial, iPlane}(idx));
            F{trialInd, iPlane, iZ} = mean(fVector{iTrial, iPlane}(fIdx, :));
            intFIdx = [(1:idx(1)-1)'; idx];
            intF{trialInd, iPlane, iZ} = mean(fVector{iTrial, iPlane}(intFIdx, :));
        end
    end
end

theta = squeeze(theta(:, min(Planes), :));
F = cell2mat(F);
accF = cell2mat(intF);
R = obj.dataTMaze.report(validTrials)';
C = obj.dataTMaze.contrastSequence(validTrials);

%%
% keyboard

%%
% get results for theta
nZ = length(zLocations);
contrastGroups = options.thCGroups;
nC = length(contrastGroups);
thRes.rArea = nan(nZ, nC);
thRes.pVal = nan(nZ, nC);
for iZ = 1:nZ
    for iC = 1:nC
        idxL = R=='L' & ismember(C, contrastGroups{iC});
        idxR = R=='R' & ismember(C, contrastGroups{iC});
        if sum(idxL) && sum(idxR) % if we have samples from both decisions
            thRes.rArea(iZ, iC) = rocArea(theta(idxL, iZ), theta(idxR, iZ));
            thRes.pVal(iZ, iC) = ranksum(theta(idxL, iZ), theta(idxR, iZ));
        end
    end
end

% get results for F and accF

contrastGroups = options.fCGroups;
nZ = length(zLocations);
nC = length(contrastGroups);
nCells = size(F, 2);
fRes.rArea = nan(nZ, nC, nCells);
fRes.pVal = nan(nZ, nC, nCells);
accFRes.rArea = nan(nZ, nC, nCells);
accFRes.pVal = nan(nZ, nC, nCells);
for iCell = 1:nCells
    if any(isnan(F(:, iCell, :)))
        continue;
    end
    for iZ = 1:nZ
        for iC = 1:nC
            idxL = R=='L' & ismember(C, contrastGroups{iC});
            idxR = R=='R' & ismember(C, contrastGroups{iC});
            if sum(idxL) && sum(idxR) % there are samples from both distributions
                fRes.rArea(iZ, iC, iCell) = rocArea(F(idxL, iCell, iZ), F(idxR, iCell, iZ));
                fRes.pVal(iZ, iC, iCell) = ranksum(F(idxL, iCell, iZ), F(idxR, iCell, iZ));
                accFRes.rArea(iZ, iC, iCell) = rocArea(accF(idxL, iCell, iZ), accF(idxR, iCell, iZ));
                accFRes.pVal(iZ, iC, iCell) = ranksum(accF(idxL, iCell, iZ), accF(idxR, iCell, iZ));
            end
        end
    end
end % iCell

resOut = struct;
resOut.options = options;
resOut.th = thRes;
resOut.F = fRes;
resOut.accF = accFRes;

%%
return

%%
figure
nRows = length(zLocations);
contrastGroups = options.thCGroups;

nContrasts = length(contrastGroups);
nColumns = nContrasts;
edges = [-30:3:30];
for iZ = 1:size(theta, 2)
    for iContrast = 1:nContrasts
        subplot(nRows, nColumns, sub2ind([nColumns, nRows], iContrast, nRows-iZ+1))
        idxL = R=='L' & ismember(C, contrastGroups{iContrast});
        hL = histogram(theta(idxL, iZ), edges);
        % hL.DisplayStyle = 'stairs';
        
        hold on;
        idxR = R=='R' & ismember(C, contrastGroups{iContrast});
        hR = histogram(theta(idxR, iZ), edges);
        % hR.DisplayStyle = 'stairs';
        xlim([-30, 30])
                    if iZ == nRows % if top Row in subplots
                title(sprintf('C = %s %%', num2str(contrastGroups{iContrast})))
            end
            rArea = NaN;
            p = NaN;
            if sum(idxL) && sum(idxR)
            rArea = rocArea(theta(idxL, iZ), theta(idxR, iZ));
            p = ranksum(theta(idxL, iZ), theta(idxR, iZ));
            end
            hx = xlabel(sprintf('ROC = %4.2f, pVal = %d', rArea, p));
            if p<0.001
                hx.FontWeight = 'bold';
                hx.Color = 'r';
            end
            if (iContrast==1) % meaninf we are in column 1
                ylabel(sprintf('z = %1.0f[cm]', zLocations(iZ)));
            end

        drawnow;
        %         pause(1)
    end
end

%%
figure
nRows = length(zLocations);
contrastGroups = options.fCGroups;
nContrasts = length(contrastGroups);
nColumns = nContrasts;
nCells = size(F, 2);
for iCell = 1:nCells
    set(gcf, 'Name', sprintf('iCell = %d', iCell));
    mm = minmax(reshape(F(:, iCell, :), 1, []));
    edges = linspace(mm(1), mm(2), 21);
    for iZ = 1:size(theta, 2)
        for iContrast = 1:nContrasts
            subplot(nRows, nColumns, sub2ind([nColumns, nRows], iContrast, nRows-iZ+1))
            idxL = R=='L' & ismember(C, contrastGroups{iContrast});
            %         hL = histogram(F(idxL, iCell, iZ), edges);
            nL = histcounts(F(idxL, iCell, iZ), edges);
            % hL.DisplayStyle = 'stairs';
            
            %         hold on;
            idxR = R=='R' & ismember(C, contrastGroups{iContrast});
            %         hR = histogram(F(idxR, iCell, iZ), edges);
            nR = histcounts(F(idxR, iCell, iZ), edges);
            
            xx = (edges(2:end)+edges(1:end-1))/2;
            bar(xx, nL, 'b')
            hold on;
            bar(xx, -nR, 'r')
            hold off;
            % hR.DisplayStyle = 'stairs';
            xlim([mm])
            %         drawnow;
            if iZ == nRows % if top Row in subplots
                title(sprintf('C = %s %%', num2str(contrastGroups{iContrast})))
            end
            rArea = NaN;
            p = NaN;
            if sum(idxL) && sum(idxR)
            rArea = rocArea(F(idxL, iCell, iZ), F(idxR, iCell, iZ));
            p = ranksum(F(idxL, iCell, iZ), F(idxR, iCell, iZ));
            end
            hx = xlabel(sprintf('ROC = %4.2f, pVal = %5.3f', rArea, p));
            if p<0.01
                hx.FontWeight = 'bold';
                hx.Color = 'r';
            end
            if (iContrast==1) % meaninf we are in column 1
                ylabel(sprintf('z = %1.0f[cm]', zLocations(iZ)));
            end
            %         pause(1)
        end
%     fitcdiscr(F(:, iCell, iZ), R)
    end
    pause
end % iCell
%%
%
% construct classifiers (d'(theta) and d'(f)) for each of the selected z values