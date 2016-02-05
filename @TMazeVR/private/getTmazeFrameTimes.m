function frameTimes = getTmazeFrameTimes(info)


% first try to load already analized data, and only if it fails do the
% full analysis

try
    % first trying to load a local copy of Timeline
    load(fullfile(info.folderTLLocal, [info.expRef, '_vrTimes']));
    return;
catch
    try
    load(fullfile(info.folderTL, [info.expRef, '_vrTimes']));
    return;
    catch
        % will load the whole movie
    end
end

try
    % first trying to load a local copy of Timeline
    load(fullfile(info.folderTLLocal, info.basenameTL));
catch
    load(fullfile(info.folderTL, info.basenameTL));
end

nInputs=length(Timeline.hw.inputs);
for iInput=1:nInputs
    if isequal(Timeline.hw.inputs(iInput).name, 'photoDiode')
        ind=iInput;
        break;
    end
end

phd=Timeline.rawDAQData(:, ind);
t = Timeline.rawDAQTimestamps;
% slight filtering to avoid double-crossings of the threshold
% phd = filtfilt([1 1 1]/3, 1, phd);
phd=(phd-min(phd))/(max(phd)-min(phd));
thrUp=0.1; % 
thrDown=0.7; %

above=phd>thrUp;
deltas=[0; diff(above)];
goUpTimes=Timeline.rawDAQTimestamps(deltas==1);

below=phd<thrDown;
deltas=[0; diff(below)];
goDownTimes=Timeline.rawDAQTimestamps(deltas==1);

stimStartTimes = [];
stimEndTimes = [];
expStartTime = [];
expEndTime = [];
instr = '';
for iUDP = 1:Timeline.mpepUDPCount
    [instr, animal, iseries, iexp, irepeat, istim, dur] = ...
        strread(Timeline.mpepUDPEvents{iUDP}, '%s %s %d %d %d %d %d' );

    switch instr{1}
        case 'ExpStart'
            expStartTime = Timeline.mpepUDPTimes(iUDP);
        case 'StimStart'
            stimStartTimes = [stimStartTimes; Timeline.mpepUDPTimes(iUDP)];
        case 'StimEnd'
            stimEndTimes = [stimEndTimes; Timeline.mpepUDPTimes(iUDP)];
        case 'ExpEnd'
            expEndTime = Timeline.mpepUDPTimes(iUDP);
    end
end

nTrials = length(stimStartTimes);
frameTimes = struct('idx', [], 't', []);
frameTimes(nTrials, 1).idx = [];
frameTimes(nTrials, 1).t= [];

% finding the folder with the T-Maze videos, might be local or on the server
vrFolder = fullfile(info.folderTLLocal, [info.expRef, '_VR']);
if ~exist(vrFolder, 'dir')
    vrFolder = fullfile(info.folderTL, [info.expRef, '_VR']);
end

% not analysing the first trial (the recording starts late), and the last
% trial (which is aborted and we don't have the video of it)

nChars = 0;
for iTrial = 2:nTrials-1
    fprintf(repmat('\b', 1, nChars));
    nChars = fprintf('%d/%d', iTrial, nTrials-1);
    vrFile = sprintf('%s_%03d.mj2', info.expRef, iTrial);
    vr = VideoReader(fullfile(vrFolder, vrFile));
    vrFrames = read(vr);
    
    [nRows, nColumns, ~, nVideoFrames] = size(vrFrames);
    syncSquare = round(double(squeeze(vrFrames(nRows, nColumns, 1, :)))/255);
    % find wich frames are not sequences of identical frames
    flipIdx = [1; find(diff(syncSquare))+1];
    
    % the logic is the following here:
    % the first UP after the StimStart UDP is the appearance of the T-Maze
    % the first DOWN after the StimStart UDP is the appearance of the Gray
    % screen (and the offset of the T-Maze from the previous trial)
    % the first DOWN after the StimEnd UDP is the offset of the last frame of that Trial 
    % the last UP before the StimEnd UDP is the onset of the last frame of that Trial 
    
    upIdx = find(goUpTimes>stimStartTimes(iTrial) & goUpTimes<stimEndTimes(iTrial));
    downIdx = find(goDownTimes>stimStartTimes(iTrial) & goDownTimes<stimEndTimes(iTrial));
    
    tNextDown = goDownTimes(find(goDownTimes>stimEndTimes(iTrial), 1, 'first'));
    
    if (length(flipIdx)-length(upIdx)-length(downIdx))
        warning('Trial %d buggy, nFrames in the video does not match number of flips of the photodiode\n Returning empty data for this trial\n', iTrial);
        continue;
    else
        % these are the frame indices in the video
        frameTimes(iTrial).idx = flipIdx;
        % these are their relative onset timestamps
        % the last timestamp is the offset of the last T-Maze frame (swithc
        % to the gray screen)
        frameTimes(iTrial).t = [reshape([goDownTimes(downIdx); goUpTimes(upIdx)], [], 1); tNextDown];
    end
%     syncSquare(flipIdx)
end
fprintf('\n');
    
save(fullfile(info.folderTLLocal, [info.expRef, '_vrTimes']), 'frameTimes');
% save(fullfile(info.folderTL, [info.expRef, '_vrTimes']));

return;

%% plotting for debugging purposes

% plot single trial with photodiode flip detections
figure
plot(t, phd);
hold on;
stem(goDownTimes(downIdx), ones(size(downIdx)),'r')
stem(goUpTimes(upIdx), ones(size(upIdx)),'g')
xlim([stimStartTimes(iTrial) stimEndTimes(iTrial)]);


% plot the whole signal with UDP timings
figure
plot(t, phd);
hold on;
stem(expStartTime, 1, 'g');
stem(stimStartTimes, ones(size(stimStartTimes)), 'g:');
stem(stimEndTimes, ones(size(stimEndTimes)), 'r:');
plot(expEndTime, ylim, 'r');



