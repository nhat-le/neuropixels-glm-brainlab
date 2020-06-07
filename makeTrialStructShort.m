function trialStruct = makeTrialStructShort(trialData)
nTrials = numel(trialData.trials_gocue_times);
trialStruct = struct();
%trialStruct(nTrials).duration = 0; % preallocate

badTrials = [];

for idTrial = 1:nTrials
    %Duration
    trialStruct(idTrial).duration = ...
        1000 * (trialData.trials_start(idTrial, 2) - trialData.trials_stim_times(idTrial) - 0.5);
%     trialStruct(idTrial).trialstart = 0;
    
    % Stimulus
    trialStruct(idTrial).stimOn = 0;
    
    if trialData.trials_left_contrast(idTrial) > 0
        trialStruct(idTrial).stimOnLeft = 0;
    end
    
    if trialData.trials_right_contrast(idTrial) > 0
        trialStruct(idTrial).stimOnRight = 0;
    end 
    
    if trialData.trials_left_contrast(idTrial) > 0.5
        trialStruct(idTrial).stimOnLeftHigh = 0;
    elseif trialData.trials_left_contrast(idTrial) > 0
        trialStruct(idTrial).stimOnLeftLow = 0;
    end
    
    if trialData.trials_right_contrast(idTrial) > 0.5
        trialStruct(idTrial).stimOnRightHigh = 0;
    elseif trialData.trials_right_contrast(idTrial) > 0
        trialStruct(idTrial).stimOnRightLow = 0;
    end
    
    trialStruct(idTrial).goCue = 1000 * (trialData.trials_gocue_times(idTrial) - ...
        trialData.trials_stim_times(idTrial));
    
    % Choice
%     if trialData.trials_choice(idTrial) > 0
%         trialStruct(idTrial).leftResponse = ...
%         1000 * (trialData.trials_response_times(idTrial) - trialData.trials_start(idTrial, 1));
%         %assert(trialStruct(idTrial).leftResponse > trialStruct(idTrial).goCue);
%     elseif trialData.trials_choice(idTrial) < 0
%         trialStruct(idTrial).rightResponse = ...
%         1000 * (trialData.trials_response_times(idTrial) - trialData.trials_start(idTrial, 1));
%     end
    if trialData.trials_choice(idTrial) > 0 && ...
            trialData.trials_decision_times(idTrial) - trialData.trials_stim_times(idTrial) < 3
        trialStruct(idTrial).leftResponse = ...
            1000 * max(trialData.trials_decision_times(idTrial) - trialData.trials_stim_times(idTrial), 0);
        
    elseif trialData.trials_choice(idTrial) < 0 && ...
            trialData.trials_decision_times(idTrial) - trialData.trials_stim_times(idTrial) < 3
        trialStruct(idTrial).rightResponse = ...
            1000 * max(trialData.trials_decision_times(idTrial) - trialData.trials_stim_times(idTrial), 0);
    end
    % Feedback
    trialStruct(idTrial).feedback = ...
        1000 * (trialData.trials_feedback_times(idTrial) - trialData.trials_stim_times(idTrial));
    if trialData.trials_feedback_types(idTrial) > 0
        trialStruct(idTrial).posFeedback = ...
            1000 * (trialData.trials_feedback_times(idTrial) - trialData.trials_stim_times(idTrial));
    else
        trialStruct(idTrial).negFeedback = ...
            1000 * (trialData.trials_feedback_times(idTrial) - trialData.trials_stim_times(idTrial));
    end

%     assert(trialStruct(idTrial).feedback < trialStruct(idTrial).duration);
    
    
    spike_trial = trialData.spikes(trialData.spikes > trialData.trials_stim_times(idTrial) & ...
        trialData.spikes < trialData.trials_start(idTrial, 2) - 0.5);
    
%     spike_trial_other = trialData.spikes_other(trialData.spikes_other > trialData.trials_start(idTrial, 1) & ...
%         trialData.spikes_other < trialData.trials_start(idTrial, 2));
    
    trialStruct(idTrial).sptrain = 1000 * (spike_trial - ...
        trialData.trials_stim_times(idTrial));
%     trialStruct(idTrial).sptrain2 = 1000 * (spike_trial_other - ...
%         trialData.trials_start(idTrial, 1));
    
    % Spike trains of others
    for j = 1:numel(trialData.spikes_other)
        spOther = trialData.spikes_other{j};
        spike_trial_other = spOther(spOther > trialData.trials_stim_times(idTrial) & ...
            spOther < trialData.trials_start(idTrial, 2) -0.5);
        othername = sprintf('sptrain%d', j);
        %disp(othername);
        trialStruct(idTrial).(othername) = 1000 * (spike_trial_other - ...
        trialData.trials_stim_times(idTrial));
    end
    
    
    if ~isempty(trialStruct(idTrial).sptrain)
        assert(min(trialStruct(idTrial).sptrain) > 0);
%         assert(max(trialStruct(idTrial).sptrain) <  trialStruct(idTrial).duration);
    end
end