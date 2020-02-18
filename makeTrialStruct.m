function trialStruct = makeTrialStruct(nTrials)
global trialData
trialStruct = struct();
trialStruct(nTrials).duration = 0; % preallocate

badTrials = [];

for idTrial = 1:nTrials
    %Duration
    trialStruct(idTrial).duration = ...
        1000 * (trialData.trials_start(idTrial, 2) - trialData.trials_start(idTrial, 1));
    trialStruct(idTrial).trialstart = 0;
    
    % Stimulus
    trialStruct(idTrial).stimOn = ...
        1000 * (trialData.trials_stim_times(idTrial) - trialData.trials_start(idTrial, 1));
    
    if trialData.trials_left_contrast(idTrial) > 0
        trialStruct(idTrial).stimOnLeft = ...
            1000 * (trialData.trials_stim_times(idTrial) - trialData.trials_start(idTrial, 1));
    end
    
    if trialData.trials_right_contrast(idTrial) > 0
        trialStruct(idTrial).stimOnRight = ...
            1000 * (trialData.trials_stim_times(idTrial) - trialData.trials_start(idTrial, 1));
    end 
    
    if trialData.trials_left_contrast(idTrial) > 0.5
        trialStruct(idTrial).stimOnLeftHigh = ...
        1000 * (trialData.trials_stim_times(idTrial) - trialData.trials_start(idTrial, 1));
    elseif trialData.trials_left_contrast(idTrial) > 0
        trialStruct(idTrial).stimOnLeftLow = ...
        1000 * (trialData.trials_stim_times(idTrial) - trialData.trials_start(idTrial, 1));
    end
    
    if trialData.trials_right_contrast(idTrial) > 0.5
        trialStruct(idTrial).stimOnRightHigh = ...
        1000 * (trialData.trials_stim_times(idTrial) - trialData.trials_start(idTrial, 1));
    elseif trialData.trials_right_contrast(idTrial) > 0
        trialStruct(idTrial).stimOnRightLow = ...
        1000 * (trialData.trials_stim_times(idTrial) - trialData.trials_start(idTrial, 1));
    end
    
    trialStruct(idTrial).goCue = ...
        1000 * (trialData.trials_gocue_times(idTrial) - trialData.trials_start(idTrial, 1));
    
    % Choice
    if trialData.trials_choice(idTrial) > 0
        trialStruct(idTrial).leftResponse = ...
        1000 * (trialData.trials_response_times(idTrial) - trialData.trials_start(idTrial, 1));
        %assert(trialStruct(idTrial).leftResponse > trialStruct(idTrial).goCue);
    elseif trialData.trials_choice(idTrial) < 0
        trialStruct(idTrial).rightResponse = ...
        1000 * (trialData.trials_response_times(idTrial) - trialData.trials_start(idTrial, 1));
    end
    
    % Feedback
    trialStruct(idTrial).feedback = ...
        1000 * (trialData.trials_feedback_times(idTrial) - trialData.trials_start(idTrial, 1));
    if trialData.trials_feedback_types(idTrial) > 0
        trialStruct(idTrial).posFeedback = ...
            1000 * (trialData.trials_feedback_times(idTrial) - trialData.trials_start(idTrial, 1));
    else
        trialStruct(idTrial).negFeedback = ...
            1000 * (trialData.trials_feedback_times(idTrial) - trialData.trials_start(idTrial, 1));
    end

    assert(trialStruct(idTrial).feedback < trialStruct(idTrial).duration);
    
    
    spike_trial = trialData.spikes(trialData.spikes > trialData.trials_start(idTrial, 1) & ...
        trialData.spikes < trialData.trials_start(idTrial, 2));
    
    trialStruct(idTrial).sptrain = 1000 * (spike_trial - trialData.trials_start(idTrial, 1));
    if ~isempty(trialStruct(idTrial).sptrain)
        assert(min(trialStruct(idTrial).sptrain) > 0);
        assert(max(trialStruct(idTrial).sptrain) <  trialStruct(idTrial).duration);
    end
end