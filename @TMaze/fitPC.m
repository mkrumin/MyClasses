function obj = fitPC(obj, modelType)

if nargin<2
    modelType = 'erf_psycho_2gammas';
end
addpath('\\zserver\Code\Psychofit\');

parsStart = [];
parsMin = [];
parsMax = [];
nFits = 10;

switch modelType
    case 'erf_psycho_2gammas'
        modelStr = 'erf\_asym';
        parsStart = [mean(obj.pcData.cc), 10, 0.1, 0.1];
        parsMin = [-100, 0, 0, 0];
        parsMax = [100, 100, 1, 1];
        [pars L]= mle_fit_psycho([obj.pcData.cc; obj.pcData.nn; obj.pcData.pp], modelType, ...
            parsStart, parsMin, parsMax, nFits);
    case 'erf_psycho'
        modelStr = 'erf';
        parsStart = [mean(obj.pcData.cc), 10, 0.05 ];
        parsMin = [min(obj.pcData.cc) 0 0];
        parsMax = [max(obj.pcData.cc) 100 0.50];
        [pars L] = mle_fit_psycho([obj.pcData.cc; obj.pcData.nn; obj.pcData.pp], modelType, ...
            parsStart, parsMin, parsMax, nFits);
    otherwise
        fprintf('''%s'' model is unsupported for fitting psychometric curve for this data\n', modelType);
        fprintf('Please use either ''erf_psycho'' or ''erf_psycho_2gammas''\n');
        fprintf('No psychometric curve was fit\n');
        return;
end

if length(obj.pcFit)==1 && isequal(obj.pcFit.modelType, '')
    iModel = 1;
else
    iModel = length(obj.pcFit)+1;
end

obj.pcFit(iModel).modelType = modelType;
obj.pcFit(iModel).modelStr = modelStr;
obj.pcFit(iModel).pars = pars;
obj.pcFit(iModel).Likelihood = L;
obj.pcFit(iModel).nFits = nFits;
obj.pcFit(iModel).parsStart = parsStart;
obj.pcFit(iModel).parsMin = parsMin;
obj.pcFit(iModel).parsMax = parsMax;

end % fitPC()
