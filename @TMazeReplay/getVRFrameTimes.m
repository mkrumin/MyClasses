function times = getVRFrameTimes(obj)

% This function will return the actual monitor frame times
% assumes the sync square was in the 'Flicker' mode

% times is a cell array of the size nStims x nRepeats
% to get monitor frame times for a specific repeat of a specific stimulus
% call times{iStim, iRepeat}

Timeline = obj.dataTL;
nInputs=length(Timeline.hw.inputs);
for iInput=1:nInputs
    if isequal(Timeline.hw.inputs(iInput).name, 'photoDiode')
        ind=iInput;
        break;
    end
end

phd=Timeline.rawDAQData(:, ind);
phd=(phd-min(phd))/(max(phd)-min(phd));
thr=0.5; % using one threshold here

above=phd>thr;
deltas=[0; diff(above)];
goingUpTimes=Timeline.rawDAQTimestamps(deltas==1);
goingDownTimes=Timeline.rawDAQTimestamps(deltas==-1);


startIdx=[];
endIdx=[];
for iUDP=1:Timeline.mpepUDPCount
    if isequal(Timeline.mpepUDPEvents{iUDP}(1:9), 'StimStart')
        startIdx=[startIdx; iUDP];
    end
    if isequal(Timeline.mpepUDPEvents{iUDP}(1:7), 'StimEnd')
        endIdx=[endIdx; iUDP];
    end
end

% the first UP after the StimStart UDP

stimOnsets=[];
for iUDP=1:length(startIdx)
    tmp=min(goingUpTimes(goingUpTimes>Timeline.mpepUDPTimes(startIdx(iUDP))));
    stimOnsets=[stimOnsets; tmp];
end
% the last DOWN before the StimEnd UDP
stimOffsets=[];
delayConst = 0.05; % a delay constant, for the case when the 'StimEnd' udp arrived before the stimulus finished to play
for iUDP=1:length(endIdx)
    tmp=max(goingDownTimes(goingDownTimes<(Timeline.mpepUDPTimes(endIdx(iUDP)) + delayConst)));
    stimOffsets=[stimOffsets; tmp];
end

for iStim = 1:length(stimOnsets)
    ups = goingUpTimes(goingUpTimes>=stimOnsets(iStim) & goingUpTimes<=stimOffsets(iStim));
    downs = goingDownTimes(goingDownTimes>=stimOnsets(iStim) & goingDownTimes<=stimOffsets(iStim));
    times{iStim} = sort([ups(:); downs(:)]);
end

% reorganizing to stimuli and repeats
% pp = ProtocolLoad(info.subject, str2num(info.expDate(info.expDate~='-')), info.exp);
% times = times(pp.seqnums);
% figure
% plot(Timeline.rawDAQTimestamps, phd);
% hold on;
% plot(times{2}, 0.5, 'r.');
% plot(times{3}, 0.5, 'g.');
