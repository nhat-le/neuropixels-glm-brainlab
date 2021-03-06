%% Files of interest
area = 'SCm';
files = dir(['GLMResults/*' area '.mat']);
fprintf('Number of files = %d\n', numel(files));


%% Define relevant windows
stimWindowMs = [0 200]; %ms
rewWindowMs = [0 200]; %ms
choiceWindowMs = [-100 0]; %ms

load(fullfile('GLMResults', files(1).name));
stimWindowPts = ws(1).stimOn.tr;
stimWindowPts = find(stimWindowPts >= stimWindowMs(1) & ...
    stimWindowPts <= stimWindowMs(2));
rewWindowPts = ws(1).posFeedback.tr;
rewWindowPts = find(rewWindowPts >= stimWindowMs(1) & ...
    rewWindowPts <= stimWindowMs(2));
choiceWindowPts = ws(1).leftResponse.tr;
choiceWindowPts = find(choiceWindowPts >= choiceWindowMs(1) & ...
    choiceWindowPts <= choiceWindowMs(2));

%% Loop through the files
sigCellsStim = [];
sigCellsRew = [];
sigCellsLeft = [];
sigCellsRight = [];

for i = 1:numel(files)
    load(fullfile('GLMResults', files(i).name));

    sigStim = isSignificant(ws, wvars, 'stimOn', stimWindowPts);
    sigRew = isSignificant(ws, wvars, 'posFeedback', rewWindowPts);
    sigLeft = isSignificant(ws, wvars, 'leftResponse', choiceWindowPts);
    sigRight = isSignificant(ws, wvars, 'rightResponse', choiceWindowPts);
    
    sigCellsStim(i) = sigStim;
    sigCellsRew(i) = sigRew;
    sigCellsLeft(i) = sigLeft;
    sigCellsRight(i) = sigRight;
end

%% Visualize
sigCells = [sigCellsStim; sigCellsRew; sigCellsLeft; sigCellsRight];
figure;
imagesc(sigCells)
colormap(flipud(gray))
title(area)
format_axes


function result = isSignificant(ws, wvars, field, tpts)
% Returns true if the cell is significant in the defined field
% Stimulus kernels
kernel = ws.(field).data;
dev = wvars.(field).data;

% Within 200 ms, is there a bin with significant increase?
activity_window = kernel(tpts);
dev_window = dev(tpts);
upper_bound = activity_window + dev_window;
lower_bound = activity_window - dev_window;

if any(upper_bound > 0 & lower_bound > 0)
    result = 1;
else
    result = 0;
end

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
%xlim([0, 700])
yticks([1,2,3,4])
yticklabels({'Stimulus', 'Reward', 'Choice Left', 'Choice Right'})
set(gca, 'FontSize', 16)
end

