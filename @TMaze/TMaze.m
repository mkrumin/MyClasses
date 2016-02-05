classdef TMaze
    properties
        SESSION;    % this structure includes all the session data
        EXP;        % this structure includes the session parameters
        ExpRef;     % a string defining the ExpRef of the session (cell array for merged sessions)
        stimSide = ''; % e.g. 'RRLLRLRLLRLLLLRRRLLRLL'
        outcome = '';   % e.g. 'WWCCFCCCCWCCCWWCCCTCCCCU'
        report = '';    %  e.g. 'RRLLLRLRLLLRLTTRLLRLRLTFU'
        contrastSequence; % signed contrast (with negative being left-side stimuli)
        pcData = struct('cc', [], 'nn', [], 'pp', [], 'sem', []); % summary of behavioral data for fitting PC
        pcFit = struct('modelType', '', 'modelStr', '', 'pars', [], 'Likelihood', [], 'nFits', [], 'parsStart', [], 'parsMin', [], 'parsMax', []);
        posUniform; % posdata on a uniform grid (vectors of the same length, time-rescaled according to trial beginning/end)
    end
    
    methods
        function obj = TMaze(ExpRef)
            obj.ExpRef = ExpRef;
            filenames = dat.expFilePath(ExpRef, 'TMaze');
            try
                data = load(filenames{1});
            catch
                try
                    data = load(filenames{2});
                catch
                    fprintf('No TMaze data found for ExpRef %s\n', ExpRef);
                    return;
                end
            end
            
            if ~isfield(data.SESSION, 'probR')
                data.SESSION.probR = data.EXP.probRight*ones(size(data.SESSION.stimSequence));
            end
            
            obj.SESSION = data.SESSION;
            obj.EXP = data.EXP;
            
            nTrials = length(obj.SESSION.allTrials);
            % excluding the last trial (without checking - bad code),
            % which is very likely to be 'USER ABORT' trial
            if isequal(obj.SESSION.allTrials(nTrials).info.outcome, 'USER ABORT')
                nTrials = nTrials -1;
            end
            obj.stimSide = '';
            obj.outcome = '';
            obj.report = '';
            for iTrial = 1:nTrials
                obj.stimSide(iTrial) = obj.SESSION.stimSequence{iTrial}(1);
                obj.outcome(iTrial) = obj.SESSION.allTrials(iTrial).info.outcome(1);
            end
            obj.report = obj.outcome;
            obj.report(obj.outcome == 'C') = obj.stimSide(obj.outcome == 'C');
            obj.report(obj.outcome == 'W') = 'R'+'L' - obj.stimSide(obj.outcome == 'W');
            obj.contrastSequence = obj.SESSION.contrastSequence(1:nTrials);
            obj.contrastSequence(obj.stimSide=='L') = -obj.contrastSequence(obj.stimSide=='L');
        end % TMaze constructor
        
    end % methods
    

end
