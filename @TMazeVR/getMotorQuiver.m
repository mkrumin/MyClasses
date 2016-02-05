function [vMap, rMap] = getMotorQuiver(obj, binCentres)

%%
dZ = median(diff(binCentres{1}));
zEdges = [binCentres{1}, binCentres{1}(end)+dZ] - dZ/2;
dTheta = median(diff(binCentres{2}));
thEdges = [binCentres{2}, binCentres{2}(end)+dTheta] - dTheta/2;

dtBall = median(diff(obj.dataBall.t));
rotGain = 1/8.3 * obj.dataTMaze.EXP.aGain * 180/pi; 
rotGain = rotGain * dtBall; % weird hack to scale properly, for some reason is only needed for the rotation!
velGain = -1/53 * obj.dataTMaze.EXP.zGain;

rotGain = rotGain/dtBall;
velGain = velGain/dtBall;

r = rotGain * obj.dataBall.rotation;
v = velGain * obj.dataBall.forward;
nTrials = length(obj.timesVRframes);
[thVector, zVector, vVector, ~] = buildVectors(obj, 1:nTrials, obj.dataBall.t, v);
[~, ~, rVector, ~] = buildVectors(obj, 1:nTrials, obj.dataBall.t, r);

coords = [zVector, thVector];
rAccumMap = buildAccumMap(coords, rVector, {zEdges, thEdges});
vAccumMap = buildAccumMap(coords, vVector, {zEdges, thEdges});
occMap = buildOccupMap(coords, {zEdges, thEdges});

% rMap = rAccumMap./occMap;
% vMap = vAccumMap./occMap;

h = ndGaussian([0.7 0.7]);
rMap = filterAndDivideMaps(occMap, rAccumMap, h);
vMap = filterAndDivideMaps(occMap, vAccumMap, h);

return;

zAxis = binCentres{1};
thAxis = binCentres{2};

figure;
subplot(1, 2, 1)
imagesc(thAxis, zAxis, rMap);
title('Rotation [deg/sec]');
xlabel('\theta [deg]');
ylabel('z [cm]');
axis xy;
colorbar;

subplot(1, 2, 2)
imagesc(thAxis, zAxis, vMap);
title('Running Speed [cm/sec]');
xlabel('\theta [deg]');
ylabel('z [cm]');
axis xy;
colorbar;

