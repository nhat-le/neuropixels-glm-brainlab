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
expt = buildGLM.registerTiming(expt, 'goCue', 'Go Cue');
expt = buildGLM.registerTiming(expt, 'response', 'Response');
expt = buildGLM.registerTiming(expt, 'feedback', 'Feedback');

expt = buildGLM.registerSpikeTrain(expt, 'sptrain', 'Our Neuron'); % Spike train!!!
%expt = buildGLM.registerSpikeTrain(expt, 'sptrain2', 'Neighbor Neuron');
% expt = buildGLM.registerValue(expt, 'coh', 'Coherence'); % information on the trial, but not associated with time
% expt = buildGLM.registerValue(expt, 'choice', 'Direction of Choice');

%% Load the Cori dataset
spikes_times = readNPY('Cori_2016-12-14\spikes.times.npy');
spikes_clusters = readNPY('Cori_2016-12-14\spikes.clusters.npy');
cluster_id = 1047;

trials_included = readNPY('Cori_2016-12-14/trials.included.npy');
trials_feedback_times = readNPY('Cori_2016-12-14/trials.feedback_times.npy');
trials_gocue_times = readNPY('Cori_2016-12-14/trials.goCue_times.npy');
trials_start = readNPY('Cori_2016-12-14/trials.intervals.npy');
trials_response_times = readNPY('Cori_2016-12-14/trials.response_times.npy');
trials_stim_times = readNPY('Cori_2016-12-14/trials.visualStim_times.npy');

% Stimulus and choice information
trials_choice = readNPY('Cori_2016-12-14/trials.response_choice.npy');
trials_left_contrast = readNPY('Cori_2016-12-14/trials.visualStim_contrastLeft.npy');
trials_right_contrast = readNPY('Cori_2016-12-14/trials.visualStim_contrastRight.npy');

trials_feedback_times_included = trials_feedback_times(trials_included);
trials_gocue_times_included = trials_gocue_times(trials_included);
trials_start_included = trials_start(trials_included,:);
trials_response_times_included = trials_response_times(trials_included);
trials_stim_times_included = trials_stim_times(trials_included);

spikes = spikes_times(spikes_clusters == cluster_id);

%% Build the trial structure
nTrials = sum(trials_included);
trialStruct = struct();
trialStruct(nTrials).duration = 0; % preallocate

badTrials = [];

for idTrial = 1:nTrials
    %Duration
    trialStruct(idTrial).duration = ...
        1000 * (trials_start_included(idTrial, 2) - trials_start_included(idTrial, 1));
    trialStruct(idTrial).trialstart = 0;
    
    % Stimulus
    if trials_left_contrast > 0.5
        trialStruct(idTrial).stimOnLeftHigh = ...
        1000 * (trials_stim_times_included(idTrial) - trials_start_included(idTrial, 1));
    elseif trials_left_contrast > 0
        trialStruct(idTrial).stimOnLeftLow = ...
        1000 * (trials_stim_times_included(idTrial) - trials_start_included(idTrial, 1));
    end
    
    if trials_right_contrast > 0.5
        trialStruct(idTrial).stimOnRightHigh = ...
        1000 * (trials_stim_times_included(idTrial) - trials_start_included(idTrial, 1));
    elseif trials_right_contrast > 0
        trialStruct(idTrial).stimOnRightLow = ...
        1000 * (trials_stim_times_included(idTrial) - trials_start_included(idTrial, 1));
    end
    
    
    trialStruct(idTrial).stimOn = ...
        1000 * (trials_stim_times_included(idTrial) - trials_start_included(idTrial, 1));
    trialStruct(idTrial).goCue = ...
        1000 * (trials_gocue_times_included(idTrial) - trials_start_included(idTrial, 1));
    trialStruct(idTrial).response = ...
        1000 * (trials_response_times_included(idTrial) - trials_start_included(idTrial, 1));
    trialStruct(idTrial).feedback = ...
        1000 * (trials_feedback_times_included(idTrial) - trials_start_included(idTrial, 1));
    assert(trialStruct(idTrial).stimOn > 0);
    assert(trialStruct(idTrial).goCue > trialStruct(idTrial).stimOn);
    try
        assert(trialStruct(idTrial).response > trialStruct(idTrial).goCue);
        assert(trialStruct(idTrial).feedback > trialStruct(idTrial).response);
    catch
        badTrials = [badTrials idTrial];
    end
    assert(trialStruct(idTrial).feedback < trialStruct(idTrial).duration);
    
    
    spike_trial = spikes(spikes > trials_start_included(idTrial, 1) & ...
        spikes < trials_start_included(idTrial, 2));
    
    trialStruct(idTrial).sptrain = 1000 * (spike_trial - trials_start_included(idTrial, 1));
    if ~isempty(trialStruct(idTrial).sptrain)
        assert(min(trialStruct(idTrial).sptrain) > 0);
        assert(max(trialStruct(idTrial).sptrain) <  trialStruct(idTrial).duration);
    end
end

expt.trial = trialStruct;

% TODO: Visualize PSTH



%% Build 'designSpec' which specifies how to generate the design matrix
% Each covariate to include in the model and analysis is specified.
dspec = buildGLM.initDesignSpec(expt);

binfun = expt.binfun;
bs = basisFactory.makeSmoothTemporalBasis('boxcar', 100, 10, binfun);
bs2 = basisFactory.makeSmoothTemporalBasis('raised cosine', 400, 40, binfun);
bs3 = basisFactory.makeSmoothTemporalBasis('boxcar', 500, 40, binfun);
bshist = basisFactory.makeSmoothTemporalBasis('boxcar', 100, 40, binfun);

% Add covariates
dspec = buildGLM.addCovariateTiming(dspec, 'trialstart', 'trialstart', 'Trial start', bs3, -100);
dspec = buildGLM.addCovariateTiming(dspec, 'stimOn', 'stimOn', 'Stimulus on', bs3, -100);
dspec = buildGLM.addCovariateTiming(dspec, 'goCue', 'goCue', 'Go Cue', bs3, -100);
dspec = buildGLM.addCovariateTiming(dspec, 'response', 'response', 'Response', bs3, -100);
dspec = buildGLM.addCovariateTiming(dspec, 'feedback', 'feedback', 'Feedback', bs3, -100);


%bs.B = 0.1 * bs.B;

%% Instantaneous Raw Signal without basis
% dspec = buildGLM.addCovariateRaw(dspec, 'LFP', [], bs);

%% Spike history
dspec = buildGLM.addCovariateSpiketrain(dspec, 'hist', 'sptrain', 'History filter', bshist);

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
dm = buildGLM.addBiasColumn(dm);

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

%dm = buildGLM.addBiasColumn(dm); % DO NOT ADD THE BIAS TERM IF USING GLMFIT

%% Least squares for initialization
tic
wInit = dm.X \ y;
toc

%% Use matRegress for Poisson regression
% it requires `fminunc` from MATLAB's optimization toolbox
addpath('matRegress')

fnlin = @nlfuns.exp; % inverse link function (a.k.a. nonlinearity)
lfunc = @(w)(glms.neglog.poisson(w, dm.X, y, fnlin)); % cost/loss function

opts = optimoptions(@fminunc, 'Algorithm', 'trust-region', ...
    'GradObj', 'on', 'Hessian','on');

[wml, nlogli, exitflag, ostruct, grad, hessian] = fminunc(lfunc, wInit, opts);
wvar = diag(inv(hessian));
ws = buildGLM.combineWeights(dm, wml);
wvar = buildGLM.combineWeights(dm, wvar);

%% Alternative maximum likelihood Poisson estimation using glmfit
%[wml, dev, stats] = glmfit(dm.X, y, 'poisson', 'link', 'log');
%wvar = stats.se.^2;


%% Visualize
fig = figure(2913); clf;
nCovar = numel(dspec.covar);
for kCov = 1:nCovar
    label = dspec.covar(kCov).label;
    subplot(nCovar, 1, kCov);
    %plot(ws.(label).tr, exp(ws.(label).data));
    %plot(ws.(label).tr, exp(ws.(label).data));
    %hline(1);
    errorbar(ws.(label).tr, ws.(label).data, sqrt(wvar.(label).data));
    title(label);
end

%% Evaluate log-likelihood and AIC
