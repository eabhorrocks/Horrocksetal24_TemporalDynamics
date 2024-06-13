function trials = genTrialStructFromAllenData(stim_table, forceStrings)

trialParams = fieldnames(stim_table);

for ifield = 1:numel(trialParams)
    tempCellArray = stim_table.(trialParams{ifield});
    if ~iscell(tempCellArray)
        tempCellArray = num2cell(tempCellArray);
    end
    
    [trials(1:length(tempCellArray)).(trialParams{ifield})] = tempCellArray{:};
end

for itrial = 1:numel(trials)
    trials(itrial).stimulus_name = {trials(itrial).stimulus_name};
end

if forceStrings
    
    for itrial = 1:numel(trials)
        for ifield = 1:numel(trialParams)
            if ~isnumeric(trials(itrial).(trialParams{ifield}))
                [newval, test] = str2num(trials(itrial).(trialParams{ifield}));
                if test
                    trials(itrial).(trialParams{ifield}) = newval;
                end
            end
        end
        
    end
end
