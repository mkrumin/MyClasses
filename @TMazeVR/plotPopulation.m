function plotPopulation(obj)

% keyboard
trialsAll = 1:length(obj.dataTMaze.contrastSequence);
trialsRR = find(obj.dataTMaze.contrastSequence'>0 & obj.dataTMaze.report == 'R');
trialsRL = find(obj.dataTMaze.contrastSequence'>0 & obj.dataTMaze.report == 'L');
trialsLR = find(obj.dataTMaze.contrastSequence'<0 & obj.dataTMaze.report == 'R');
trialsLL = find(obj.dataTMaze.contrastSequence'<0 & obj.dataTMaze.report == 'L');
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

        trace = nanmean(fMatrix(trialsRR, :));
%         trace = trace - min(trace);
%         trace = trace/max(trace);
        tracesRR(iRow, :) = trace;
%         t2pRR(iRow) = find(trace == max(trace));

        trace = nanmean(fMatrix(trialsRL, :));
        tracesRL(iRow, :) = trace;
%         t2pRL(iRow) = find(trace == max(trace));

        trace = nanmean(fMatrix(trialsLR, :));
        tracesLR(iRow, :) = trace;
%         t2pLR(iRow) = find(trace == max(trace));

        trace = nanmean(fMatrix(trialsLL, :));
        tracesLL(iRow, :) = trace;
%         t2pLL(iRow) = find(trace == max(trace));
    end
end

%%
fStd = [0.1 3];
fTracesRR = imgaussfilt(tracesRR, fStd);
fTracesRL = imgaussfilt(tracesRL, fStd);
fTracesLR = imgaussfilt(tracesLR, fStd);
fTracesLL = imgaussfilt(tracesLL, fStd);

minValues = min([fTracesRR, fTracesRL, fTracesLR, fTracesLL], [], 2);
% maxValues = max([fTracesRR, fTracesRL, fTracesLR, fTracesLL], [], 2);
maxValues = max([fTracesRR, fTracesRL, fTracesLR, fTracesLL], [], 2);

fTracesRR = bsxfun(@rdivide, bsxfun(@minus, fTracesRR, minValues), maxValues-minValues);
fTracesRL = bsxfun(@rdivide, bsxfun(@minus, fTracesRL, minValues), maxValues-minValues);
fTracesLR = bsxfun(@rdivide, bsxfun(@minus, fTracesLR, minValues), maxValues-minValues);
fTracesLL = bsxfun(@rdivide, bsxfun(@minus, fTracesLL, minValues), maxValues-minValues);

[maxRR, t2pRR] = max(fTracesRR, [], 2);
[maxRL, t2pRL] = max(fTracesRL, [], 2);
[maxLR, t2pLR] = max(fTracesLR, [], 2);
[maxLL, t2pLL] = max(fTracesLL, [], 2);

thr = 1;
idxRR = find(maxRR >= thr);
idxRL = find(maxRL >= thr);
idxLR = find(maxLR >= thr);
idxLL = find(maxLL >= thr);

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

