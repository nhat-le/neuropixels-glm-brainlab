%% Files of interest here
files = dir('ws_wvar\*.mat');


%%

allVISCollection = [];
allACACollection = [];
allSCCollection = [];

allVISNumbers = [];
allACANumbers = [];
allSCNumbers = [];

%% Define relevant windows
stimWindowMs = [0 200]; %ms
rewWindowMs = [0 200]; %ms
choiceWindowMs = [-100 0]; %ms

stimWindowPts = ws_all(1).stimOn.tr;
stimWindowPts = find(stimWindowPts >= stimWindowMs(1) & ...
    stimWindowPts <= stimWindowMs(2));
rewWindowPts = ws_all(1).posFeedback.tr;
rewWindowPts = find(rewWindowPts >= stimWindowMs(1) & ...
    rewWindowPts <= stimWindowMs(2));
choiceWindowPts = ws_all(1).leftResponse.tr;
choiceWindowPts = find(choiceWindowPts >= choiceWindowMs(1) & ...
    choiceWindowPts <= choiceWindowMs(2));

%%

for i = 1:numel(files)
    
    load(fullfile(files(i).folder, files(i).name))
    str = files(i).name;
    startIndex = regexp(str, '[0-9]+[^0-9]*');
    area = str(startIndex+4:end-4);
    fprintf('Processing file %d of %d: %s...\n', i, numel(files), area);
    
    % To visualize results of GLM fitting for multiple cells


    %% Collect all stimOn durations
    tpts = 1:200; %ws_all(1).stimOn.tr;
    cellsStim = findSignificantCells(ws_all, wvars_all, 'stimOn', stimWindowPts);
    cellsRew = findSignificantCells(ws_all, wvars_all, 'posFeedback', rewWindowPts);
    cellsLeft = findSignificantCells(ws_all, wvars_all, 'leftResponse', choiceWindowPts);
    cellsRight = findSignificantCells(ws_all, wvars_all, 'rightResponse', choiceWindowPts);
    
    
    %% Put in collection
    collection = [cellsStim; cellsRew; cellsLeft; cellsRight];
    
    if contains(area, 'VISp')
        allVISCollection = [allVISCollection collection];
        allVISNumbers = [allVISNumbers size(collection, 2)];
    elseif contains(area, 'ACA')
        allACACollection = [allACACollection collection];
        allACANumbers = [allACANumbers size(collection, 2)];
    elseif contains(area, 'SC')
        allSCCollection = [allSCCollection collection];
        allSCNumbers = [allSCNumbers size(collection, 2)];
    end

    
    
    
end

%% Superposition of stimulus and reward
allVISSuper = allVISCollection(1,:) .* allVISCollection(2,:);
allACCSuper = allACACollection(1,:) .* allACACollection(2,:);
allSCSuper = allSCCollection(1,:) .* allSCCollection(2,:);


% allVISCollection = [allVISCollection; allVISSuper];
% allACACollection = [allACACollection; allACCSuper];
% allSCCollection = [allSCCollection; allSCSuper];



%% Visualize results
figure(1);
subplot(311);
imagesc(allVISCollection)
colormap(flipud(gray))
title('V1')




format_axes
% Plot session demarcations
cumallVIS = cumsum(allVISNumbers);
for i = 1:numel(cumallVIS)
    hold on
    vline(cumallVIS(i), 'r')
end

subplot(312)
imagesc(allACACollection)
title('ACC')
format_axes
cumallACA = cumsum(allACANumbers);
for i = 1:numel(cumallACA)
    hold on
    vline(cumallACA(i), 'r')
end



subplot(313)
imagesc(allSCCollection)
title('SC')
format_axes
xlabel('Cell #')
cumallSC = cumsum(allSCNumbers);
for i = 1:numel(cumallSC)
    hold on
    vline(cumallSC(i), 'r')
end






function result = findSignificantCells(ws, wvars, field, tpts)
    cells = [];

    for i = 1:numel(ws)
        % Stimulus kernels
        kernel = ws(i).(field).data;
        dev = wvars(i).(field).data;

        % Within 200 ms, is there a bin with significant increase?
        activity_window = kernel(tpts);
        dev_window = dev(tpts);
        upper_bound = activity_window + dev_window;
        lower_bound = activity_window - dev_window;

        if any(upper_bound > 0 & lower_bound > 0)
            cells = [cells i];
        end

    end
    
    % Make 0/1 result vector
    result = zeros(1, numel(ws));
    result(cells) = 1;


end

function format_axes
xlim([0, 700])
yticks([1,2,3,4])
yticklabels({'Stimulus', 'Reward', 'Choice Left', 'Choice Right'})
set(gca, 'FontSize', 16)
end

