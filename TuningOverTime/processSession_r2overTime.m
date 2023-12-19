function processSession_r2overTime(inputFileName,outputFileName,dataDir)

%% load session

load(fullfile(dataDir,inputFileName))

%% process trials
% fitler trials to the ones we want to analyse and sort according to
% behavioural state

trialsSpeed2D = trials.Speed2D;
trialsSpeed2D(1) = [];
tsd = trialsSpeed2D;

% split trials by state according to wheel data
run_idx = find(cellfun(@(x) prop(x>0.5)>=0.75 & mean(x)>3, {tsd.WheelSpeed}));
stat_idx = find(cellfun(@(x) prop(x<3)>=0.75 & mean(x)<0.5, {tsd.WheelSpeed}));

[tsd.runFlag] = deal([nan]);

[tsd(stat_idx).runFlag] = deal([0]);
[tsd(run_idx).runFlag] = deal([1]);

temp_tsd = tsd([tsd.numDots1]==573); % remove blnk trials.
temp_tsd = temp_tsd(~isnan([temp_tsd.runFlag]));
allStimConds = [vertcat(temp_tsd.VelX1), vertcat(temp_tsd.Contrast1), vertcat(temp_tsd.runFlag)];
[uniqueConds, ~, ic] = unique(allStimConds, 'rows');
tally = accumarray(ic, 1);
Result = [uniqueConds tally];

temp_tsd = temp_tsd([temp_tsd.Contrast1]==1); % only full contrast trials


%% get binned spike counts (10ms bins) for each trial

for itrial =1 :numel(temp_tsd)
    temp_tsd(itrial).start_time = temp_tsd(itrial).PDstart;
    temp_tsd(itrial).absVel = abs(temp_tsd(itrial).VelX1);
end

for iunit = 1:numel(units)
    units(iunit).spiketimes = units(iunit).spike_times;
end


% statTsd = temp_tsd([temp_tsd.runFlag]==0);
% runTsd = temp_tsd([temp_tsd.runFlag]==1);
options.intervalStart = -0.2;
options.intervalEnd = 1.8;
options.binSpacing=0.01;

[anUnits, cond] = getBinnedSpikeCounts(temp_tsd, units, {'absVel', 'runFlag'}, options);

for iunit = 1:numel(units)
    units(iunit).allSpikes = anUnits(iunit).allSpikes;
end

%% downsample to min # trials available
% set so equal # trials per condition

minTrial = min(min(cellfun(@(x) size(x,2), units(1).allSpikes)));

tempUnits = units;

for iunit = 1:numel(tempUnits)
    tempUnits(iunit).allSpikes = cellfun(@(x) x(:,1:minTrial), tempUnits(iunit).allSpikes,'UniformOutput', false);
end


%% get cross-validated R2 (tuning strength) over time
% options
r2opts.nPerms = 10;
r2opts.randFlag = 1;
r2opts.validMeans = 1;
r2opts.kval = 3;
r2opts.nShuffle = 100;

binSize = 0.01;
windowLength = 0.2;
binVector = -0.2:0.01:1.8;
binVector = round(binVector,2);
binLengths = windowLength/binSize;
binStarts = 1:1:181; % first bin of each window

% get the tuning strength and significance over time
for iunit = 1:numel(tempUnits)
    % stationary trials
    sca = tempUnits(iunit).allSpikes(:,1)';
    tempUnits(iunit).stat_ints = calcR2overtime_paper(sca,binLengths,binStarts,r2opts);
    % locomotion trials
    sca = tempUnits(iunit).allSpikes(:,2)';
    tempUnits(iunit).run_ints = calcR2overtime_paper(sca,binLengths,binStarts,r2opts);
end

%% use tuning strength vectors to get tuning timing info

for iunit = 1:numel(tempUnits)

    [tempUnits(iunit).stat_tunedStart, tempUnits(iunit).stat_tunedFinish,...
        tempUnits(iunit).stat_nBouts, tempUnits(iunit).stat_tunedDuration, ...
        tempUnits(iunit).stat_tuneFlagWeak, tempUnits(iunit).stat_tuneFlagMed,...
        tempUnits(iunit).stat_tuneFlagStrong] = ...
        getTuningTimeInfo_paper([tempUnits(iunit).stat_ints.R2], [tempUnits(iunit).stat_ints.R2_p]);


    [tempUnits(iunit).run_tunedStart, tempUnits(iunit).run_tunedFinish,...
        tempUnits(iunit).run_nBouts, tempUnits(iunit).run_tunedDuration, ...
        tempUnits(iunit).run_tuneFlagWeak, tempUnits(iunit).run_tuneFlagMed,...
        tempUnits(iunit).run_tuneFlagStrong] = ...
        getTuningTimeInfo_paper([tempUnits(iunit).run_ints.R2], [tempUnits(iunit).run_ints.R2_p]);
end


  %% save processed data
    
    save(fullfile(dataDir,outputFileName),'tempUnits','minTrial')
    
end
