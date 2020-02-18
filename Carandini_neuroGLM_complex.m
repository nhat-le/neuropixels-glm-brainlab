%% Load the raw data
%rawData = load('exampleData.mat'); % run tutorial_exampleData to generate this
% nTrials = rawData.nTrials; % number of trials
unitOfTime = 'ms';
binSize = 1; % TODO some continuous observations might need up/down-sampling if binSize is not 1!?

global trialData

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
%expt = buildGLM.registerSpikeTrain(expt, 'sptrain2', 'Neighbor Neuron');
% expt = buildGLM.registerValue(expt, 'coh', 'Coherence'); % information on the trial, but not associated with time
% expt = buildGLM.registerValue(expt, 'choice', 'Direction of Choice');

%% Load the Cori dataset
folder = 'Cori_2016-12-14';
spikes_times = readNPY(fullfile(folder, 'spikes.times.npy'));
spikes_clusters = readNPY(fullfile(folder, 'spikes.clusters.npy'));
clusters_annotation = readNPY(fullfile(folder, 'clusters._phy_annotation.npy'));
good_clusters = find(clusters_annotation >= 2) - 1;

trials_included = readNPY(fullfile(folder, 'trials.included.npy'));
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

% IMPORTANT: Only include indicated trials
trialData.trials_feedback_times = trialData.trials_feedback_times(trials_included);
trialData.trials_gocue_times = trialData.trials_gocue_times(trials_included);
trialData.trials_start = trialData.trials_start(trials_included,:);
trialData.trials_response_times = trialData.trials_response_times(trials_included);
trialData.trials_stim_times = trialData.trials_stim_times(trials_included);
trialData.trials_feedback_types = trialData.trials_feedback_types(trials_included);
trialData.trials_choice = trialData.trials_choice(trials_included);
trialData.trials_left_contrast = trialData.trials_left_contrast(trials_included);
trialData.trials_right_contrast = trialData.trials_right_contrast(trials_included);

%% Build the trial structure
cluster_id = 72; %good_clusters(13);
trialData.spikes = spikes_times(spikes_clusters == cluster_id);

nTrials = sum(trials_included);
trialStruct = makeTrialStruct(nTrials);

expt.trial = trialStruct;

%% Visualize PSTH




%% Build 'designSpec' which specifies how to generate the design matrix
% Each covariate to include in the model and analysis is specified.
dspec = buildGLM.initDesignSpec(expt);

binfun = expt.binfun;
bs = basisFactory.makeSmoothTemporalBasis('boxcar', 100, 10, binfun);
bsStim = basisFactory.makeSmoothTemporalBasis('boxcar', 300, 10, binfun);
bs2 = basisFactory.makeSmoothTemporalBasis('raised cosine', 400, 40, binfun);
bs3 = basisFactory.makeSmoothTemporalBasis('boxcar', 400, 40, binfun);
bshist = basisFactory.makeSmoothTemporalBasis('boxcar', 100, 40, binfun);


% Add covariates
%dspec = buildGLM.addCovariateTiming(dspec, 'trialstart', 'trialstart', 'Trial start', bs3, -100);
% dspec = buildGLM.addCovariateTiming(dspec, 'stimOnLeftHigh', 'stimOnLeftHigh', ...
%     'Stimulus on left (high contrast)', bs3, -100);
% dspec = buildGLM.addCovariateTiming(dspec, 'stimOnLeftLow', 'stimOnLeftLow', ...
%     'Stimulus on left (low contrast)', bs3, -100);
% dspec = buildGLM.addCovariateTiming(dspec, 'stimOnRightHigh', 'stimOnRightHigh', ...
%     'Stimulus on right (high contrast)', bs3, -100);
% dspec = buildGLM.addCovariateTiming(dspec, 'stimOnRightLow', 'stimOnRightLow', ...
%     'Stimulus on right (low contrast)', bs3, -100);
dspec = buildGLM.addCovariateTiming(dspec, 'stimOnLeft', 'stimOnLeft', ...
     'Stimulus onset', bsStim);
 dspec = buildGLM.addCovariateTiming(dspec, 'stimOnRight', 'stimOnRight', ...
     'Stimulus onset', bsStim);


dspec = buildGLM.addCovariateTiming(dspec, 'goCue', 'goCue', 'Go Cue', bs3, -100);
dspec = buildGLM.addCovariateTiming(dspec, 'leftResponse', 'leftResponse', 'Left Response', bs3, -100);
dspec = buildGLM.addCovariateTiming(dspec, 'rightResponse', 'rightResponse', 'Right Response', bs3, -100);
dspec = buildGLM.addCovariateTiming(dspec, 'posFeedback', 'posFeedback', 'posFeedback', bs3, -100);
dspec = buildGLM.addCovariateTiming(dspec, 'negFeedback', 'negFeedback', 'negFeedback', bs3, -100);

%bs.B = 0.1 * bs.B;

%% Instantaneous Raw Signal without basis
% dspec = buildGLM.addCovariateRaw(dspec, 'LFP', [], bs);

%% Spike history
dspecHist = buildGLM.addCovariateSpiketrain(dspec, 'hist', 'sptrain', 'History filter', bshist);

%% Coupling filter
%dspec = buildGLM.addCovariateSpiketrain(dspec, 'coupling', 'sptrain2', 'Coupling from neuron 2');

%% Duration boxcar
% dspec = buildGLM.addCovariateBoxcar(dspec, 'dots', 'dotson', 'dotsoff', 'Motion dots stim');

%% Acausal Timing Event
% bs = basisFactory.makeSmoothTemporalBasis('boxcar', 300, 8, binfun);
%                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           
% dspec = buildGLM.addCovariateTiming(dspec, 'saccade', [], [], bs, offset);

%% Coherence
% a box car that depends on the coh value
% bs = basisFactory.makeSmoothTemporalBasis('raised cosine', 200, 10, binfun);
% stimHandle = @(trial, expt) trial.coh * basisFactory.boxcarStim(binfun(trial.dotson), binfun(trial.dotsoff), binfun(trial.duration));

% dspec = buildGLM.addCovariate(dspec, 'cohKer', 'coh-dep dots stimulus', stimHandle, bs);

%% 2-D eye position
% bs = basisFactory.makeSmoothTemporalBasis('raised cosine', 40, 4, binfun);
% dspec = buildGLM.addCovariateRaw(dspec, 'eyepos', [], bs);

%buildGLM.summarizeDesignSpec(dspec); % print out the current configuration

%% Compile the data into 'DesignMatrix' structure
trialIndices = 1:nTrials-1; % use all trials except the last one
dm = buildGLM.compileSparseDesignMatrix(dspec, trialIndices);
dm_hist = buildGLM.compileSparseDesignMatrix(dspecHist, trialIndices);

%% Visualize the design matrix
endTrialIndices = cumsum(binfun([expt.trial(trialIndices).duration]));
X = dm.X(1:endTrialIndices(3),:);
mv = max(abs(X), [], 1); mv(isnan(mv)) = 1;
X = bsxfun(@times, X, 1 ./ mv);
figure(742); clf; imagesc(X');
%buildGLM.visualizeDesignMatrix(dm, 1); % optionally plot the first trial

%% Get the spike trains back to regress against
y = buildGLM.getBinnedSpikeTrain(expt, 'sptrain', dm.trialIndices);

%% Do some processing on the design matrix
%dm = buildGLM.removeConstantCols(dm);
% colIndices = buildGLM.getDesignMatrixColIndices(dspec, 'LFP');
% dm = buildGLM.zscoreDesignMatrix(dm, [colIndices{:}]);

dm = buildGLM.addBiasColumn(dm); % DO NOT ADD THE BIAS TERM IF USING GLMFIT
dm_hist = buildGLM.addBiasColumn(dm_hist);

%% Least squares for initialization
tic
wInit = dm.X \ y;
wInitHist = dm_hist.X \ y;
toc

%% Use matRegress for Poisson regression
% it requires `fminunc` from MATLAB's optimization toolbox
fnlin = @nlfuns.exp; % inverse link function (a.k.a. nonlinearity)
lfunc = @(w)(glms.neglog.poisson(w, dm.X, y, fnlin)); % cost/loss function
lfuncHist = @(w)(glms.neglog.poisson(w, dm_hist.X, y, fnlin)); % cost/loss function

opts = optimoptions(@fminunc, 'Algorithm', 'trust-region', ...
    'GradObj', 'on', 'Hessian','on');

[wml, nlogli, exitflag, ostruct, grad, hessian] = fminunc(lfunc, wInit, opts);
[wmlH, nlogliH, exitflagH, ostructH, gradH, hessianH] = fminunc(lfuncHist, wInitHist, opts);
wvar = diag(inv(hessian));
wvarH = diag(inv(hessianH));

%% Alternative maximum likelihood Poisson estimation using glmfit
%[wml, dev, stats] = glmfit(dm.X, y, 'poisson', 'link', 'log');
%wvar = stats.se.^2;

%% Determine log-likelihood and AIC
AIC_model0 = nlogli + size(dm.X, 2);
AIC_model1 = nlogliH + size(dm_hist.X, 2);

%% Visualize
ws = buildGLM.combineWeights(dm_hist, wmlH);
wvars = buildGLM.combineWeights(dm_hist, wvarH);

fig = figure(2913); clf;
nCovar = numel(dspecHist.covar);
for kCov = 1:nCovar
    label = dspecHist.covar(kCov).label;
    subplot(nCovar, 1, kCov);
    %plot(ws.(label).tr, (ws.(label).data));
    errorbar(ws.(label).tr, ws.(label).data, sqrt(wvars.(label).data));
    hline(0)
    ylim([-10 10])
    title(label);
end

%% Compare different stimulus types
% figure;
% subplot(1,2,1)
% title('Stimulus');
% l1 = plot(ws.stimOnLeftHigh.tr, exp(ws.stimOnLeftHigh.data));
% hold on
% l2 = plot(ws.stimOnLeftLow.tr, exp(ws.stimOnLeftLow.data));
% l3 = plot(ws.stimOnRightHigh.tr, exp(ws.stimOnRightHigh.data));
% l4 = plot(ws.stimOnRightLow.tr, exp(ws.stimOnRightLow.data));
% legend([l1, l2, l3, l4], {'LeftHigh', 'LeftLow', 'RightHigh', 'RightLow'});
% 
% 
% subplot(1,2,2)
% title('Response')
% plot(ws.leftResponse.tr, exp(ws.leftResponse.data));
% hold on
% plot(ws.rightResponse.tr, exp(ws.rightResponse.data));


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










