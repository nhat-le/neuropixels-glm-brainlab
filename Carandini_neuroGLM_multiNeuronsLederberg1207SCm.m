%% Load the raw data
%rawData = load('exampleData.mat'); % run tutorial_exampleData to generate this
% nTrials = rawData.nTrials; % number of trials
unitOfTime = 'ms';
binSize = 1; % TODO some continuous observations might need up/down-sampling if binSize is not 1!?

%% Specify the fields to load
expt = buildGLM.initExperiment(unitOfTime, binSize, 'Carandini', []);

% Events
expt = buildGLM.registerTiming(expt, 'trialstart', 'Trial start');
expt = buildGLM.registerTiming(expt, 'stimOn', 'Stimulus onset');
expt = buildGLM.registerTiming(expt, 'stimOnLeft', 'Left Stimulus onset');
expt = buildGLM.registerTiming(expt, 'stimOnRight', 'Right Stimulus onset');
expt = buildGLM.registerTiming(expt, 'goCue', 'Go Cue');
expt = buildGLM.registerTiming(expt, 'response', 'Response');
expt = buildGLM.registerTiming(expt, 'feedback', 'Feedback');
expt = buildGLM.registerTiming(expt, 'posFeedback', 'Positive Feedback');
expt = buildGLM.registerTiming(expt, 'negFeedback', 'Negative Feedback');
expt = buildGLM.registerTiming(expt, 'stimOnLeftHigh', 'stimOnLeftHigh');
expt = buildGLM.registerTiming(expt, 'stimOnLeftLow', 'stimOnLeftLow');
expt = buildGLM.registerTiming(expt, 'stimOnRightHigh', 'stimOnRightHigh');
expt = buildGLM.registerTiming(expt, 'stimOnRightLow', 'stimOnRightLow');
expt = buildGLM.registerTiming(expt, 'leftResponse', 'leftResponse');
expt = buildGLM.registerTiming(expt, 'rightResponse', 'rightResponse');
expt = buildGLM.registerSpikeTrain(expt, 'sptrain', 'Our Neuron'); % Spike train!!!

%% Load the Cori dataset
folder = 'Data\Lederberg_2017-12-07';
spikes_times = readNPY(fullfile(folder, 'spikes.times.npy'));
spikes_clusters = readNPY(fullfile(folder, 'spikes.clusters.npy'));
clusters_annotation = readNPY(fullfile(folder, 'clusters._phy_annotation.npy'));
good_clusters = find(clusters_annotation >= 2) - 1;

% List of trials to include
trialIndices = readNPY(fullfile(folder, 'trials.included.npy'));
trialIndices = find(trialIndices);

trialData.trials_feedback_times = readNPY(fullfile(folder, 'trials.feedback_times.npy'));
trialData.trials_feedback_types = readNPY(fullfile(folder, 'trials.feedbackType.npy'));
trialData.trials_gocue_times = readNPY(fullfile(folder, 'trials.goCue_times.npy'));
trialData.trials_start = readNPY(fullfile(folder, 'trials.intervals.npy'));
trialData.trials_response_times = readNPY(fullfile(folder, 'trials.response_times.npy'));
trialData.trials_stim_times = readNPY(fullfile(folder, 'trials.visualStim_times.npy'));

% Stimulus and choice information
trialData.trials_choice = readNPY(fullfile(folder, 'trials.response_choice.npy'));
trialData.trials_left_contrast = readNPY(fullfile(folder, 'trials.visualStim_contrastLeft.npy'));
trialData.trials_right_contrast = readNPY(fullfile(folder, 'trials.visualStim_contrastRight.npy'));

% Load brain areas
load cell_areas_Lederberg1207.mat

%% Build the trial structure
cluster_id = 103; %good_clusters(13);
cluster_id_other = SCm;
trialData.spikes = spikes_times(spikes_clusters == cluster_id);

trialData.spikes_other = cell(1, numel(cluster_id_other));
for i = 1:numel(cluster_id_other)
    trialData.spikes_other{i} = spikes_times(spikes_clusters == cluster_id_other(i));
end

trialStruct = makeTrialStruct(trialData);

expt.trial = trialStruct;

%% Build 'designSpec' which specifies how to generate the design matrix
% Each covariate to include in the model and analysis is specified.
dspec = buildGLM.initDesignSpec(expt);

binfun = expt.binfun;
bs = basisFactory.makeSmoothTemporalBasis('boxcar', 100, 10, binfun);
bsStim = basisFactory.makeSmoothTemporalBasis('boxcar', 300, 10, binfun);
bs2 = basisFactory.makeSmoothTemporalBasis('raised cosine', 400, 40, binfun);
bs3 = basisFactory.makeSmoothTemporalBasis('boxcar', 400, 40, binfun);
bshist = basisFactory.makeSmoothTemporalBasis('boxcar', 300, 40, binfun);

dspec = buildGLM.addCovariateTiming(dspec, 'stimOn', 'stimOn', ...
      'Stimulus onset', bsStim);


dspec = buildGLM.addCovariateTiming(dspec, 'goCue', 'goCue', 'Go Cue', bs3, -100);
dspec = buildGLM.addCovariateTiming(dspec, 'leftResponse', 'leftResponse', 'Left Response', bs3, -100);
dspec = buildGLM.addCovariateTiming(dspec, 'rightResponse', 'rightResponse', 'Right Response', bs3, -100);
dspec = buildGLM.addCovariateTiming(dspec, 'posFeedback', 'posFeedback', 'posFeedback', bs3, -100);
dspec = buildGLM.addCovariateTiming(dspec, 'negFeedback', 'negFeedback', 'negFeedback', bs3, -100);

%bs.B = 0.1 * bs.B;

%% Compile the data into 'DesignMatrix' structure
wml_all = zeros(251, numel(cluster_id_other));
wvar_all = zeros(251, numel(cluster_id_other));
dspec2 = buildGLM.addCovariateSpiketrain(dspec, 'hist', 'sptrain', 'History filter', bshist);
labels = {dspec2.covar.label};
labels{2,1} = {};
ws_all = struct(labels{:});
wvars_all = struct(labels{:});

%% Get the spike trains back to regress against
for i = 1:numel(cluster_id_other)
    name = sprintf('sptrain%d', i);
    dspecCell = buildGLM.addCovariateSpiketrain(dspec, 'hist', name, 'History filter', bshist);
    dmCell = buildGLM.compileSparseDesignMatrix(dspecCell, trialIndices);
    
    fprintf('Processing cell %d of %d\n...', i, numel(cluster_id_other));
    y = buildGLM.getBinnedSpikeTrain(expt, name, dmCell.trialIndices);


    %% Do some processing on the design matrix
    dmCell = buildGLM.addBiasColumn(dmCell); % DO NOT ADD THE BIAS TERM IF USING GLMFIT

    %% Least squares for initialization
    tic
    wInit = dmCell.X \ y;
    toc

    %% Use matRegress for Poisson regression
    % it requires `fminunc` from MATLAB's optimization toolbox
    fnlin = @nlfuns.exp; % inverse link function (a.k.a. nonlinearity)
    lfunc = @(w)(glms.neglog.poisson(w, dmCell.X, y, fnlin)); % cost/loss function

    opts = optimoptions(@fminunc, 'Algorithm', 'trust-region', ...
        'GradObj', 'on', 'Hessian','on');

    [wml, nlogli, exitflag, ostruct, grad, hessian] = fminunc(lfunc, wInit, opts);
    wvar = diag(inv(hessian));
    
    wml_all(:, i) = wml;
    wvar_all(:, i) = wvar;
    
    % Collect the weights
    ws = buildGLM.combineWeights(dmCell, wml_all(:, i));
    wvars = buildGLM.combineWeights(dmCell, wvar_all(:, i));
    
    ws_all(i) = ws;
    wvars_all(i) = wvars;

    

end
%% Visualize
for i = 1:numel(cluster_id_other)
    fig = figure(i); clf;
    nCovar = numel(dspecCell.covar);
    count = 1;
    lst = [1, 2, 3, 4, 5, 6, 7];
    for kCov = lst
        label = dspecCell.covar(kCov).label;
        subplot(numel(lst), 1, count);
        %plot(ws.(label).tr, (ws.(label).data));
        errorbar(ws_all(i).(label).tr, ws_all(i).(label).data, sqrt(wvars_all(i).(label).data));
        hline(0)
        %ylim([-10 10])
        title(label);
        count = count + 1;
        ylim([-10 10])
    end
end

save('ws_wvar_Lederberg1207SCm.mat', 'ws_all', 'wvars_all');

function spikeTrialTimes = splitSpikeTimesToTrials(spikes, tstarts, tends)
spikeTrialTimes = {};
for i = 1:numel(tstarts)
    spikeTrial = spikes(spikes > tstarts(i) & spikes < tends(i));
    spikeTrialTimes{i} = spikeTrial;
end

end

function [x, y, spikeBinnedCounts] = BinSpikeTimes(spikes, tstarts, window, nbins)
    spikeBinnedCounts = zeros(numel(tstarts), nbins - 1);
    for i = 1:numel(tstarts)
        spikeTrial = spikes(spikes > tstarts(i) + window(1) & spikes < tstarts(i) + window(2));
        edges = linspace(tstarts(i) + window(1), tstarts(i) + window(2), nbins);
        counts = histogram(spikeTrial, edges);
        spikeBinnedCounts(i,:) = counts.Values;
    end
    y = 1:numel(tstarts);
    x = linspace(window(1), window(2), nbins);
end










