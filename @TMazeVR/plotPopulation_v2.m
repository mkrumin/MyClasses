function plotPopulation_v2(obj)
%%
% keyboard
trialsAll = 1:length(obj.dataTMaze.contrastSequence);
trialsRR = find(obj.dataTMaze.contrastSequence'>0 & obj.dataTMaze.report == 'R');
trialsRL = find(obj.dataTMaze.contrastSequence'>0 & obj.dataTMaze.report == 'L');
trialsLR = find(obj.dataTMaze.contrastSequence'<0 & obj.dataTMaze.report == 'R');
trialsLL = find(obj.dataTMaze.contrastSequence'<0 & obj.dataTMaze.report == 'L');

cc = unique(obj.dataTMaze.contrastSequence);
for iCC = 1:length(cc)
    trialsGoR{iCC} = find(obj.dataTMaze.contrastSequence' == cc(iCC) & obj.dataTMaze.report == 'R');
    trialsGoL{iCC} = find(obj.dataTMaze.contrastSequence' == cc(iCC) & obj.dataTMaze.report == 'L');
end
trialsGoR{length(cc)+1} = find(obj.dataTMaze.report == 'R');
trialsGoL{length(cc)+1} = find(obj.dataTMaze.report == 'L');

nGroups = length(trialsGoR);

nPoints = 100;
iRow = 0;
cellIdx = [];
planeIdx = [];

for iPlane = obj.Planes
    tData = obj.times2p{iPlane}';
    nCells = obj.nROIs(iPlane);
    for iCell = 1:nCells
        roiType = obj.data2p{iPlane}.ROI.CellClasses{iCell};
%         if roiType ~= 's'
%             % only process 'soma' ROIs
%             continue;
%         end
        
        fData = obj.data2p{iPlane}.F(:, iCell);
        if any(isnan(fData))
            continue;
        end
        
        cellIdx = cat(1, cellIdx, iCell);
        planeIdx = cat(1, planeIdx, iPlane);
        
        iRow = iRow + 1;

        [~, ~, fMatrix] = buildTraces(obj, trialsAll, tData, fData, nPoints);
%         fMatrix = fMatrix-nanmin(fMatrix(:));
%         fMatrix = fMatrix/nanmax(fMatrix(:));

        trace = nanmean(fMatrix(trialsRR, :), 1);
%         trace = trace - min(trace);
%         trace = trace/max(trace);
        tracesRR(iRow, :) = trace;
%         t2pRR(iRow) = find(trace == max(trace));

        trace = nanmean(fMatrix(trialsRL, :), 1);
        tracesRL(iRow, :) = trace;
%         t2pRL(iRow) = find(trace == max(trace));

        trace = nanmean(fMatrix(trialsLR, :), 1);
        tracesLR(iRow, :) = trace;
%         t2pLR(iRow) = find(trace == max(trace));

        trace = nanmean(fMatrix(trialsLL, :), 1);
        tracesLL(iRow, :) = trace;
%         t2pLL(iRow) = find(trace == max(trace));

        for iGroup = 1:nGroups
            tracesGoR{iGroup}(iRow, :) = nanmean(fMatrix(trialsGoR{iGroup}, :), 1);
            tracesGoL{iGroup}(iRow, :) = nanmean(fMatrix(trialsGoL{iGroup}, :), 1);
        end
    end
end

%%
fStd = [0.1 3];
fTracesRR = imgaussfilt(tracesRR, fStd);
fTracesRL = imgaussfilt(tracesRL, fStd);
fTracesLR = imgaussfilt(tracesLR, fStd);
fTracesLL = imgaussfilt(tracesLL, fStd);
for iGroup = 1:nGroups
    fTracesGoR{iGroup} = imgaussfilt(tracesGoR{iGroup}, fStd);
    fTracesGoL{iGroup} = imgaussfilt(tracesGoL{iGroup}, fStd);
end

minValues = min([fTracesRR, fTracesRL, fTracesLR, fTracesLL], [], 2);
% maxValues = max([fTracesRR, fTracesRL, fTracesLR, fTracesLL], [], 2);
maxValues = max([fTracesRR, fTracesRL, fTracesLR, fTracesLL], [], 2);

fTracesRR = bsxfun(@rdivide, bsxfun(@minus, fTracesRR, minValues), maxValues-minValues);
fTracesRL = bsxfun(@rdivide, bsxfun(@minus, fTracesRL, minValues), maxValues-minValues);
fTracesLR = bsxfun(@rdivide, bsxfun(@minus, fTracesLR, minValues), maxValues-minValues);
fTracesLL = bsxfun(@rdivide, bsxfun(@minus, fTracesLL, minValues), maxValues-minValues);
for iGroup = 1:nGroups
    fTracesGoR{iGroup} = bsxfun(@rdivide, bsxfun(@minus, fTracesGoR{iGroup}, minValues), maxValues-minValues);
    fTracesGoL{iGroup} = bsxfun(@rdivide, bsxfun(@minus, fTracesGoL{iGroup}, minValues), maxValues-minValues);
end

[maxRR, t2pRR] = max(fTracesRR, [], 2);
[maxRL, t2pRL] = max(fTracesRL, [], 2);
[maxLR, t2pLR] = max(fTracesLR, [], 2);
[maxLL, t2pLL] = max(fTracesLL, [], 2);
for iGroup = 1:nGroups
    [maxGoR(:, iGroup), t2pGoR(:, iGroup)] = max(fTracesGoR{iGroup}, [], 2);
    [maxGoL(:, iGroup), t2pGoL(:, iGroup)] = max(fTracesGoL{iGroup}, [], 2);
end

thr = 1;
idxRR = find(maxRR >= thr);
idxRL = find(maxRL >= thr);
idxLR = find(maxLR >= thr);
idxLL = find(maxLL >= thr);

idxGoR = find(maxRR >= thr | maxLR >= thr);
idxGoL = find(maxRL >= thr | maxLL >= thr);

%% plotting
cc(nGroups+1) = NaN;
figure('Name', obj.expRef);
colormap jet;
for iGroup = 1:nGroups
[~, I] = sort(t2pGoL(idxGoR, iGroup));
IR = idxGoR(I);
[~, I] = sort(t2pGoL(idxGoL, iGroup));
IL = idxGoL(I);

subplot(2, nGroups, iGroup);
imagesc(cat(1, fTracesGoL{iGroup}(IR, :), fTracesGoL{iGroup}(IL, :)));
axis xy;
title(num2str(cc(iGroup)))

[~, I] = sort(t2pGoR(idxGoR, iGroup));
IR = idxGoR(I);
[~, I] = sort(t2pGoR(idxGoL, iGroup));
IL = idxGoL(I);

subplot(2, nGroups, nGroups+iGroup);
imagesc(cat(1, fTracesGoR{iGroup}(IR, :), fTracesGoR{iGroup}(IL, :)));
axis xy;

end

%%
return;
%% plotting 

% figure('Name', 'unsorted');
% subplot(2,2,1);
% imagesc(fTracesRR);
% axis xy
% title('stim R, went R');
% subplot(2,2,2);
% imagesc(fTracesRL);
% axis xy
% title('stim R, went L');
% subplot(2,2,3);
% imagesc(fTracesLR);
% axis xy
% title('stim L, went R');
% subplot(2,2,4);
% imagesc(fTracesLL);
% axis xy
% title('stim L, went L');
% colormap jet

figure('Name', obj.expRef);
[~, I] = sort(t2pRR(idxRR));
I = idxRR(I);

subplot(4,4,1);
imagesc(fTracesRR(I, :));
title('stim R, went R');
subplot(4,4,2);
imagesc(fTracesRL(I, :));
title('stim R, went L');
subplot(4,4,5);
imagesc(fTracesLR(I, :));
title('stim L, went R');
subplot(4,4,6);
imagesc(fTracesLL(I, :));
title('stim L, went L');
% colormap jet

[~, I] = sort(t2pLL(idxLL));
I = idxLL(I);

subplot(4, 4, 11);
imagesc(fTracesRR(I, :));
title('stim R, went R');
subplot(4, 4, 12);
imagesc(fTracesRL(I, :));
title('stim R, went L');
subplot(4, 4, 15);
imagesc(fTracesLR(I, :));
title('stim L, went R');
subplot(4, 4, 16);
imagesc(fTracesLL(I, :));
title('stim L, went L');
% colormap jet

[~, I] = sort(t2pRL(idxRL));
I = idxRL(I);

subplot(4, 4, 3);
imagesc(fTracesRR(I, :));
title('stim R, went R');
subplot(4, 4, 4);
imagesc(fTracesRL(I, :));
title('stim R, went L');
subplot(4, 4, 7);
imagesc(fTracesLR(I, :));
title('stim L, went R');
subplot(4, 4, 8);
imagesc(fTracesLL(I, :));
title('stim L, went L');
% colormap jet

[~, I] = sort(t2pLR(idxLR));
I = idxLR(I);

subplot(4, 4, 9);
imagesc(fTracesRR(I, :));
title('stim R, went R');
subplot(4, 4, 10);
imagesc(fTracesRL(I, :));
title('stim R, went L');
subplot(4, 4, 13);
imagesc(fTracesLR(I, :));
title('stim L, went R');
subplot(4, 4, 14);
imagesc(fTracesLL(I, :));
title('stim L, went L');
colormap jet

annotation('line', [0.51 0.51], [0 1]);
annotation('line', [0 1], [0.51 0.51]);

%% 
ch = get(gcf, 'Children');
for iCh = 1:length(ch)
    if isequal(get(ch(iCh), 'Type'), 'axes')
        set(ch(iCh), 'XTickLabel', '', 'Box', 'off');
        set(ch(iCh), 'CLim', [0 1]);
        set(ch(iCh), 'TitleFontWeight', 'normal');
        set(ch(iCh), 'YDir', 'normal');
%         set(ch(iCh), 'CLimMode', 'auto');
        
    end
end

