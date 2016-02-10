function obj = getPCData(obj)

cc = unique(obj.contrastSequence);
nContrasts = length(cc);
cc = reshape(cc, 1, nContrasts);
nn = nan(1, nContrasts);
pp = nan(1, nContrasts);
sem = nan(1, nContrasts);
confInt = nan(2, nContrasts);
for iContrast = 1:nContrasts
    idx = obj.contrastSequence == cc(iContrast);
    % excluding time-outs and fails from the calculation
    idx = idx & (obj.report =='R' | obj.report == 'L')'; 
    nn(iContrast) = sum(idx);
%     pp(iContrast) = sum(obj.report(idx) == 'R')/nn(iContrast);
    [pp(iContrast), confInt(:, iContrast)] = binofit(sum(obj.report(idx) == 'R'), nn(iContrast), 0.05);
    % calculation of the SEM is based on the Binomial
    % distribution var(x) = n*p*q formula;
    sem(iContrast) = sqrt(pp(iContrast)*(1-pp(iContrast))/nn(iContrast));
end

obj.pcData.cc = cc;
obj.pcData.nn = nn;
obj.pcData.pp = pp;
obj.pcData.conf = confInt;
obj.pcData.sem = sem;


end % gatPCData
