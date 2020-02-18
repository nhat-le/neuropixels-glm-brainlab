function visualizePSTH()
global trialData
% Make flattened arrays of spikes
spikes_all = {trialStruct.sptrain};
spikesFlat = spikes;

% Bin
[x, y, spikeBinnedCounts] = BinSpikeTimes(spikesFlat, trials_stim_times(:,1), [-1 3], 400);

% Plot
h = imagesc(spikeBinnedCounts);
xData = get(h, 'XData');
yData = get(h, 'YData');
xWorldLimit = [-1 3];
yWorldLimit = yData;
I = imref2d(size(spikeBinnedCounts), xWorldLimit, yWorldLimit);
imagesc(x, y, spikeBinnedCounts)
%imshow(spikeBinnedCounts, I);
colormap(flipud(gray));