function PlotMultipanelZ(obj, iPlane, iROI)

stimRightTrials = obj.dataTMaze.contrastSequence > 0;
stimLeftTrials = obj.dataTMaze.contrastSequence < 0;
stimZeroTrials = obj.dataTMaze.contrastSequence == 0;
highContrastTrials = abs(obj.dataTMaze.contrastSequence) >= 25;
lowContrastTrials = abs(obj.dataTMaze.contrastSequence) < 25;
lowContrastTrials = lowContrastTrials & ~stimZeroTrials;
wentRightTrials = (obj.dataTMaze.report == 'R')';
wentLeftTrials =  (obj.dataTMaze.report == 'L')';

nTrials = length(wentRightTrials);
validTrials = ~isnan(obj.timesTmazeOnset(1:nTrials));
rows = struct;
columns = struct;
rows(1).idx = validTrials & wentLeftTrials;
rows(1).label = 'Go L';
rows(2).idx = validTrials & wentRightTrials;
rows(2).label = 'Go R';
rows(3).idx = validTrials;
rows(3).label = 'All';

columns(1).idx = validTrials & stimLeftTrials;
columns(1).label = 'Stim L All';
columns(2).idx = validTrials & stimLeftTrials & highContrastTrials;
columns(2).label = 'Stim L high C';
columns(3).idx = validTrials & stimLeftTrials & lowContrastTrials;
columns(3).label = 'Stim L low C';
columns(4).idx = validTrials & stimZeroTrials;
columns(4).label = 'Stim L 0%';
columns(5).idx = validTrials & stimRightTrials & lowContrastTrials;
columns(5).label = 'Stim R low C';
columns(6).idx = validTrials & stimRightTrials & highContrastTrials;
columns(6).label = 'Stim R high C';
columns(7).idx = validTrials & stimRightTrials;
columns(7).label = 'Stim R All';
columns(8).idx = validTrials;
columns(8).label = 'All';

%             figure
nSamples = 100;
nr = length(rows);
nc = length(columns);

zTimes = [];
z = [];
for iTrial = 1:nTrials
    if ~isempty(obj.timesVRframes(iTrial).t)
        zTimes = cat(1, zTimes, obj.timesVRframes(iTrial).t(2:end-1));
        frameIdx = obj.timesVRframes(iTrial).idx(2:end);
        z = cat(1, z, -obj.dataTMaze.SESSION.allTrials(iTrial).posdata(frameIdx, 3));
    end
end

eta = getZScaledETA(obj.times2p{iPlane}, obj.data2p{iPlane}.F(:, iROI), ...
    zTimes, z, [obj.timesTmazeOnset, obj.timesTmazeOffset], nSamples);

for iRow = 1:nr
    for iColumn = 1:nc
        
        subplot(nr, nc, (iRow-1)*nc+iColumn)
        idx = rows(iRow).idx & columns(iColumn).idx;
        %                     eta = getTimeScaledETA(obj.times2p{iPlane}, obj.data2p{iPlane}.F(:, iROI), [obj.timesTmazeOnset(idx), obj.timesTmazeOffset(idx)], nSamples);
        %                     imagesc(eta.t, 1:size(eta.all, 1), eta.all);
        imagesc(eta.z, 1:sum(idx), eta.all(find(idx), :));
        if iRow == 1
            title(columns(iColumn).label);
        end
        if iColumn == 1
            ylabel(rows(iRow).label);
        end
        if iRow == nr && iColumn == nc
            title(sprintf('Cell %d/%d', iROI, size(obj.data2p{iPlane}.F, 2)));
        end
        set(gca, 'CLim', [min(obj.data2p{iPlane}.F(:, iROI)), max(obj.data2p{iPlane}.F(:, iROI))]);
        
    end
end

end % PlotMulitpanel()
