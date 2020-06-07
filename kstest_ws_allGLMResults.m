%% Load the raw data
unitOfTime = 'ms';
binSize = 1;

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

%% Load the master table
[master,txt,~] = xlsread('labeled_clusters_all_session.csv');
txt = txt(2:end,:);

%%
area = 'ACA';
dirs = txt(:,end);
% Directories containing the area
subdirs = unique(dirs(master(:,2) >= 2 & strcmp(txt(:,7), area), :));
Nunits = numel(dirs(master(:,2) >= 2 & strcmp(txt(:,7), area), :));
tic;

pvallst = [];
totalSpikes = [];
folderNames = {};
clustIDs = [];

pd = makedist('Uniform','Lower',0,'Upper',1);
ctr = 1;


for ifolder = 7% 1:numel(subdirs)
    folderName = subdirs{ifolder};
    subtbl = master(master(:,2) >= 2 & strcmp(dirs, folderName) & ...
        strcmp(txt(:,7), area), :);
    good_clusters = subtbl(:, 4);

    %% Load the Cori dataset
    folder = fullfile('Data', folderName);
    spikes_times = readNPY(fullfile(folder, 'spikes.times.npy'));
    spikes_clusters = readNPY(fullfile(folder, 'spikes.clusters.npy'));

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
    trialData.trials_decision_times = readNPY(fullfile(folder, 'trials.decision_times.npy'));
    
    %% Build the trial structure
    trialData.spikes = nan;

    trialData.spikes_other = cell(1, numel(good_clusters));
    for i = 1:numel(good_clusters)
        trialData.spikes_other{i} = spikes_times(spikes_clusters == good_clusters(i));
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

    %% Get the spike trains back to regress against
    for i = 17%1:numel(good_clusters)
        name = sprintf('sptrain%d', i);
        dspecCell = buildGLM.addCovariateSpiketrain(dspec, 'hist', name, 'History filter', bshist);
        dmCell = buildGLM.compileSparseDesignMatrix(dspecCell, trialIndices);

        fprintf('Processing cell %d of %d...\n', i, numel(good_clusters));
        y = buildGLM.getBinnedSpikeTrain(expt, name, dmCell.trialIndices);

        dmCell = buildGLM.addBiasColumn(dmCell);
        wfilename = sprintf('%s_unit%d_%s.mat', folderName, good_clusters(i), area)
        load(fullfile('GLMResults', wfilename));
        
        % Transform
        %pd = makedist('Uniform','Lower',0,'Upper',1);

        if size(dmCell.X, 2) == numel(wml)
            ypred = exp(dmCell.X * wml);
            events = find(y);

    %         ycumsum = cumsum(ypred);
    %         tvals2 = ycumsum(events(2:end)) - ycumsum(events(1:end-1));

            tvals = [];
            for eventid = 1:numel(events) - 1
                tvals(eventid) = sum(ypred((events(eventid)) : (events(eventid + 1) - 1)));
            end

            zvals = 1 - exp(-tvals);
            if numel(zvals) == 0
                p = -1;
            else

                % ks statistic
                [h,p] = kstest(zvals, 'CDF', pd);
                fprintf('   p-val = %.4f\n', p);
            end
        else
            p = -2;
            disp('Dimension mismatch');
        end
        pvallst(ctr) = p;
        totalSpikes(ctr) = sum(y);
        foldernames{ctr} = folderName;
        clustIDs(ctr) = good_clusters(i);
        ctr = ctr + 1;
        disp(numel(pvallst))
        toc;
        
        %% KS plot
        plot(sort(zvals), (1:numel(zvals)) / numel(zvals), 'LineWidth', 2)
        hold on
        err = 1.36 / (numel(zvals))^0.5;
        plot([0, 1], [0, 1], 'k--', 'LineWidth', 2)
        plot([0, 1], [err, 1+err], 'k--', 'LineWidth', 2)
        plot([0, 1], [-err, 1-err], 'k--', 'LineWidth', 2)
        mymakeaxis(gca, 'yticks',[0, 0.5, 1], 'xticks', [0, 0.5, 1],...
    'y_label', 'Cumulative Distribution Function',...
    'x_label', 'Quantiles', 'font_size', 20)
        %set(gca, 'FontSize', 20);
%         
%         %% Histogram
%         nbins = 50;
%         hist(zvals, nbins)
%         expheight = numel(zvals)/nbins;
%         std = sqrt(numel(zvals) * (1/nbins) * (1-1/nbins)); 
%         hold on;
%         plotHorizontal(expheight);
%         plotHorizontal(expheight + std * 1.96);
%         plotHorizontal(expheight - std * 1.96);

    end
    clear trialData;
toc;
end


%% Save
savefilename = 'kstestresults_VISp_05032020.mat';
if ~exist(savefilename, 'file')
    save(savefilename, 'pvallst', 'totalSpikes', 'foldernames', 'clustIDs')
end
%%
% load('kstestresults_ACA.mat');
% for ifolder = 1:numel(dirs)
%     folderName = subdirs{ifolder};
%     subtbl = master(master(:,2) >= 2 & strcmp(dirs, folderName) & ...
%         strcmp(txt(:,7), area), :);
%     good_clusters = subtbl(:, 4);
% toc;
% end



