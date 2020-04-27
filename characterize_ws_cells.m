% To visualize results of GLM fitting for multiple cells

% Load data
load('Cori1214\ws_wvar_VISp_Cori1214VISp.mat');
fields = fieldnames(ws_all);
wsVIS = ws_all;
wvarVIS = wvars_all;

load('Cori1214\ws_wvar_VISp_Cori1214ACA.mat');
wsACA = ws_all;
wvarACA = wvars_all;

load('Cori1214\ws_wvar_VISp_Cori1214VISp.mat');
wsSC = ws_all;
wvarSC = wvars_all;



%% Collect all stimOn durations
tpts = 1:200; %ws_all(1).stimOn.tr;
cellsStimVIS = findSignificantCells(wsVIS, wvarVIS, 'stimOn', tpts);
cellsRewVIS = findSignificantCells(wsVIS, wvarVIS, 'posFeedback', tpts);
cellsLeftVIS = findSignificantCells(wsVIS, wvarVIS, 'leftResponse', tpts);
cellsRightVIS = findSignificantCells(wsVIS, wvarVIS, 'rightResponse', tpts);



cellsStimACA = findSignificantCells(wsACA, wvarACA, 'stimOn', tpts);
cellsRewACA = findSignificantCells(wsACA, wvarACA, 'posFeedback', tpts);
cellsLeftACA = findSignificantCells(wsACA, wvarACA, 'leftResponse', tpts);
cellsRightACA = findSignificantCells(wsACA, wvarACA, 'rightResponse', tpts);

cellsStimSC = findSignificantCells(wsSC, wvarSC, 'stimOn', tpts);
cellsRewSC = findSignificantCells(wsSC, wvarSC, 'posFeedback', tpts);
cellsLeftSC = findSignificantCells(wsSC, wvarSC, 'leftResponse', tpts);
cellsRightSC = findSignificantCells(wsSC, wvarSC, 'rightResponse', tpts);

function cells = findSignificantCells(ws, wvars, field, tpts)
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


end


