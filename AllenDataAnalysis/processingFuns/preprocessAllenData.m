function [units, trials, runInfo, pupilInfo, sessionInfo] = preprocessAllenData(sessionFile, requiredStims, options)
%% TO DO:
% 'ITI' trials?

if ~iscell(requiredStims)
    error('requiredStims must be a cell or cell array')
end
% load session file
load(sessionFile); % this is just the spike times dict atm

% generate units struct with spike times and generic info
units = genUnitStructFromAllenData(unitInformation,unitIndex,spikeTimes);
for iunit=1:numel(units)
    units(iunit).genotype = {sessionInfo.genotype};
end
% process stimulus information into trials struct
trials = genTrialStructFromAllenData(stim_table, false);

if ~ismember(requiredStims, 'all')
% extract only required trials
trials = trials(ismember([trials.stimulus_name], requiredStims));
end

if isstruct(pupilInfo)
% calculate pupil dilation as mean of width and height
pupilInfo.pupilDilation = (pupilInfo.pupilWidth + pupilInfo.pupilHeight)./2;
pupilInfo.pupilArea = pi.*(pupilInfo.pupilWidth).*(pupilInfo.pupilHeight);
end
% get mean run speed and pupil dilation
trials = getRunAndPupilInfo(trials, runInfo, pupilInfo, options);

% save(outputFile, 'units', 'trials', 'runInfo', 'pupilInfo');


end % end of function
