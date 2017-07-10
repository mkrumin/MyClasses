function [preF, postF, preT, postT] = buildVectorsPrePost(obj, trialIdx, tData, fData)

preDur = 1; % [sec] - duration of pre-trial-onset response

preF = cell(length(trialIdx), 1);
preT = cell(length(trialIdx), 1);
postF = cell(length(trialIdx), 1);
postT = cell(length(trialIdx), 1);
for trialNum = 1:length(trialIdx)
    iTrial = trialIdx(trialNum);
    tt = obj.timesVRframes(iTrial).t;
    if isempty(tt)
        continue;
    end
    tt = tt(2:end-1);
    idxPost = find(tData>=tt(1) & tData <=tt(end));
    postF{trialNum} = fData(idxPost, :);
    postT{trialNum} = tData(idxPost) - tt(1);
    idxPre = find(tData<tt(1) & tData >= (tt(1)-preDur));
    preF{trialNum} = fData(idxPre, :);
    preT{trialNum} = tData(idxPre) - tt(1);
end


end % buildVectors()
