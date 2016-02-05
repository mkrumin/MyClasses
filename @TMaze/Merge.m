function obj = Merge(obj, addObj)
nTrialsOriginal = obj.nTrials;
nTrialsAdditional = addObj.nTrials;
obj.SESSION.stimSequence = cat(2, obj.SESSION.stimSequence(1:nTrialsOriginal), addObj.SESSION.stimSequence(1:nTrialsAdditional));
obj.SESSION.contrastSequence = cat(1, obj.SESSION.contrastSequence(1:nTrialsOriginal), addObj.SESSION.contrastSequence(1:nTrialsAdditional));
obj.SESSION.optiStimAmplitude = cat(1, obj.SESSION.optiStimAmplitude(1:nTrialsOriginal), addObj.SESSION.optiStimAmplitude(1:nTrialsAdditional));
obj.SESSION.optiStimAmplitude2 = cat(1, obj.SESSION.optiStimAmplitude2(1:nTrialsOriginal), addObj.SESSION.optiStimAmplitude2(1:nTrialsAdditional));
obj.SESSION.probR = cat(2, obj.SESSION.probR(1:nTrialsOriginal), addObj.SESSION.probR(1:nTrialsAdditional));
obj.SESSION.Log = cat(2, obj.SESSION.Log(1:nTrialsOriginal), addObj.SESSION.Log(1:nTrialsAdditional));
obj.SESSION.allTrials = cat(2, obj.SESSION.allTrials(1:nTrialsOriginal), addObj.SESSION.allTrials(1:nTrialsAdditional));

if ischar(obj.ExpRef)
    obj.ExpRef = {obj.ExpRef};
end
if ischar(addObj.ExpRef)
    addObj.ExpRef = {addObj.ExpRef};
end
obj.ExpRef = cat(1, obj.ExpRef, addObj.ExpRef);

obj.stimSide = cat(2, obj.stimSide, addObj.stimSide);
obj.outcome = cat(2, obj.outcome, addObj.outcome);
obj.report = cat(2, obj.report, addObj.report);
obj.contrastSequence = cat(1, obj.contrastSequence, addObj.contrastSequence); % signed contrast (with negative being left-side stimuli)

obj.pcData = struct('cc', [], 'nn', [], 'pp', [], 'sem', []);
obj.pcFit = struct('modelType', '', 'modelStr', '', 'pars', [], 'Likelihood', [], 'nFits', [], 'parsStart', [], 'parsMin', [], 'parsMax', []);


end % Merge
