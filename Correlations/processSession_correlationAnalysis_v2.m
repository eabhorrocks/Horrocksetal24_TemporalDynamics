function processSession_correlationAnalysis_v2(inputFileName,outputFileName,dataDir)

%% load session
addpath('X:\ibn-vision\DATA\PROJECTS\AllenData\gitcopy\AnalyseData\PopulationAnalyses')

load(fullfile(dataDir,inputFileName))

%% filter trials to those for analysis and classify according to behavioural state

trialDur = 1;
tolerance = 0.05; % within 50ms
tsd = trials.Speed2D;
trialdurs = [tsd.PDend]-[tsd.PDstart];
invalidDurs_idx = abs(trialdurs-1)>tolerance;
tsd(invalidDurs_idx)=[];

tsd = tsd([tsd.Contrast1]==1); % only using full contrast trials

stateTrialType = 'changepoints'; % 'normal', 'strict', 'changepoints'

switch stateTrialType
    case 'normal'

        % normal paper criteria
        stat_idx = find(cellfun(@(x) prop(x<3)>=0.75 & mean(x)<0.5, {tsd.WheelSpeed}));
        run_idx = find(cellfun(@(x) prop(x>0.5)>=0.75 & mean(x)>3, {tsd.WheelSpeed}));


    case 'strict'
        % stricter criteria
        stat_idx = find(cellfun(@(x) prop(x<0.5)>=1 & mean(x)<0.5, {tsd.WheelSpeed}));
        run_idx = find(cellfun(@(x) prop(x>0.5)>=0.9 & mean(x)>3, {tsd.WheelSpeed}));

    case 'changepoints'


wheelOn=[]; wheelOff=[];
for iwheel = 1:numel(wheel)
    wheel(iwheel).rawSpeedInterp = cat(1,0, wheel(iwheel).rawSpeedInterp);
    wheel(iwheel).eTimeInterp = cat(1,wheel(iwheel).eTimeInterp(1)-0.01, wheel(iwheel).eTimeInterp);
    wheel(iwheel).rawSpeedInterp(end+1) = 0;
    wheel(iwheel).eTimeInterp(end+1) = wheel(iwheel).eTimeInterp(end)+0.01;


    wheelOn(iwheel) = wheel(iwheel).eTimeInterp(1);
    wheelOff(iwheel) = wheel(iwheel).eTimeInterp(end);
end

wheelSpeed = cat(1,wheel.rawSpeedInterp);
wheelZSpd = zscore(wheelSpeed);
wheelTime = cat(1,wheel.eTimeInterp);
data = wheelZSpd;
zThres = 0.005; % moving standard deviations exceeded/fell below an empirical threshold of 0.005
timestamps = wheelTime;
inSampleRate = 100;
smoothWin = 2; % moving standard deviation of speed (2s in Lohani)
changeDur=5; % minimum duration of the state change in seconds (5s in Lohani)
timeBetween=0.5; % minimum time between off and the next on in seconds

% lohani use 3s buffer from start/end points for sustained


[~,OnTStamp ,OffTStamp ] =changepoints(data, zThres,timestamps,inSampleRate,...
    smoothWin,changeDur,timeBetween);

% check if intervals occur inbetween distinct recordings - if so, set off
% to end time on wheel/on to 1st wheel time
nEpochs = numel(OnTStamp);

for iepoch = 1:nEpochs
    for iwheel = 1:numel(wheel)
        % if epoch goes passed where wheel recording ends (for this stim set)
        % set its end time to wheel end + create a new epoch that starts
        % with next wheel
        if OnTStamp(iepoch)<wheelOff(iwheel) && OffTStamp(iepoch)>wheelOff(iwheel)
            tempOff = OffTStamp(iepoch);
            % end this epoch at end of wheel
            OffTStamp(iepoch)=wheelOff(iwheel);
            % create new epoch at start of next wheel and with original end
            OnTStamp(numel(OnTStamp)+1) = wheelOn(iwheel+1);  
            OffTStamp(numel(OffTStamp)+1) = tempOff;
        end
    end

end

locoInterval_meanSpeed=[];
locoInterval_duration=[];
locoInterval_timeSinceLast=[];

for i=1:numel(OnTStamp)
    locoInterval_meanSpeed(i) = mean(wheelSpeed(find(wheelTime==OnTStamp(i)):find(wheelTime==OffTStamp(i))));
    locoInterval_duration(i) = OffTStamp(i)-OnTStamp(i);
    if i>1
        locoInterval_timeSinceLast(i) = OnTStamp(i)-OffTStamp(i-1);
    end
end

% remove intervals that don't meet criteria
toDelIdx = find([locoInterval_meanSpeed]<3 | [locoInterval_duration]<5);

OnTStamp(toDelIdx)=[];
OffTStamp(toDelIdx)=[];

% find trials that start and end within a locomotion epoch
% loop through each locomotion epoch and find valid trials
nEpochs = numel(OnTStamp);

allStartTimes = cat(1,tsd.PDstart)-0.2;
allEndTimes = cat(1,tsd.PDend)+0.8;

runIdx = [];
% must start and finish within the time range

OnTStamp = OnTStamp+0.5; % remove smaall buffer from start and end of epochs
OffTStamp = OffTStamp-0.5;

for iepoch = 1:nEpochs
    runIdx_temp = find(allStartTimes>OnTStamp(iepoch) & allStartTimes<OffTStamp(iepoch)...
        & allEndTimes>OnTStamp(iepoch) & allEndTimes<OffTStamp(iepoch));

    runIdx=cat(1,runIdx,runIdx_temp);
end
        run_idx = runIdx;
        % use stricter stationary trial criteria
        stat_idx = find(cellfun(@(x) prop(x<0.5)>=1 & mean(x)<0.5, {tsd.WheelSpeed}));
end


[tsd.runFlag] = deal([nan]);

[tsd(stat_idx).runFlag] = deal([0]);
[tsd(run_idx).runFlag] = deal([1]);


tsd = tsd(~isnan([tsd.runFlag]));

% dot-motion trials
dmtrials = tsd([tsd.numDots1]==573);

for itrial =1 :numel(dmtrials)
    dmtrials(itrial).start_time = dmtrials(itrial).PDstart;
    dmtrials(itrial).absVel = abs(dmtrials(itrial).VelX1);
end

for iunit = 1:numel(units)
    units(iunit).spiketimes = units(iunit).spike_times;
end


%% get binned spike counts
options.intervalStart = -0.2;
options.intervalEnd = 1.8;
options.binSpacing=0.01; % 10ms bin spacing.

[anUnits, cond] = getBinnedSpikeCounts(dmtrials, units, {'absVel', 'runFlag'}, options);

for iunit = 1:numel(units)
    units(iunit).allSpikes = anUnits(iunit).allSpikes;
end

% blank trials
blankTrials  = tsd([tsd.numDots1]==0);

[anUnits, cond] = getBinnedSpikeCounts(blankTrials, units, {'runFlag'}, options);

for iunit = 1:numel(units)
    units(iunit).blankSpikes = anUnits(iunit).allSpikes;
end


%% downsample trials and calculate mean rates
% equal trial # for each condition

minTrial = min(cellfun(@(x) size(x,2), units(1).allSpikes(:)));

for iunit = 1:numel(units)
    units(iunit).allSpikes = cellfun(@(x) x(:,1:minTrial), units(iunit).allSpikes, 'UniformOutput', false);
end

minBlankTrial = min(cellfun(@(x) size(x,2), units(1).blankSpikes(:)));

for iunit = 1:numel(units)
    units(iunit).blankSpikes = cellfun(@(x) x(:,1:minBlankTrial), units(iunit).blankSpikes, 'UniformOutput', false);
end


% filter units for mean firing rate
for iunit =1:numel(units)
    units(iunit).dmFR = mean(cellfun(@(x) mean(x,'all'), units(iunit).allSpikes),1).*(1/options.binSpacing);
    units(iunit).blankFR = mean(cellfun(@(x) mean(x,'all'), units(iunit).blankSpikes),1).*(1/options.binSpacing);
end


%% get good units
minFR = 1;

goodUnits = units([units.isi_viol]<=0.1...
    & [units.amplitude_cutoff]<=0.1 & [units.amplitude]>=50 & all(cat(1,units.dmFR)>minFR,2)'); %% & [units.dmFR]>=1);

nTrials = size(goodUnits(1).allSpikes{1,1},2);
nInts = size(goodUnits(1).allSpikes{1,1},1);
nStim = numel(goodUnits(1).allSpikes(:,1));

%% get total, signal and noise correlations

windowSize = 0.2./options.binSpacing;

for istate = 1:2

    nca = nan(numel(goodUnits), numel(goodUnits), nInts-(windowSize-1));
    sca = nca;
    tca = nca;

    for iint = 1:nInts-(windowSize-1)

        % pre-process data
        spikeCountArray = nan(nTrials,nStim,numel(goodUnits));

        for iunit=1:numel(goodUnits)
            tempSpikeCounts = cellfun(@(x) sum(x(iint:iint+19,:),1), goodUnits(iunit).allSpikes(:,istate), 'UniformOutput', false);
            %spikeCountArray(:,iunit) = [tempSpikeCounts{:}]';
            spikeCountArray(:,:,iunit) = cell2mat(tempSpikeCounts)';
        end

        % get total correlation
        sca_reshape = reshape(spikeCountArray,nTrials*nStim,numel(goodUnits));
        totalCorrArray = corr(sca_reshape,'type','Pearson');

        % get signal correlation by shuffling within-condition trials
        for ishuffle = 1:10
            clear dataStruct
            sca_shuffle = spikeCountArray(randperm(nTrials),:,:);
            sca_shuffle_reshape = reshape(sca_shuffle,nTrials*nStim,numel(goodUnits));
            shuffleCorrArray(:,:,ishuffle) = corr(sca_reshape, sca_shuffle_reshape, 'type', 'Pearson');

        end

        signalCorrArray = mean(shuffleCorrArray,3); % final signal correlation estimate is mean of shuffles
        noiseCorrArray = totalCorrArray - signalCorrArray; % noise correlation estimate is total - signal

        tca(:,:,iint) = totalCorrArray;
        sca(:,:,iint) = signalCorrArray;
        nca(:,:,iint) = noiseCorrArray;


    end

    if istate == 1
        stat.noiseCorrArray = nca;
        stat.signalCorrArray = sca;
        stat.totalCorrArray = tca;

    elseif istate == 2
        run.noiseCorrArray = nca;
        run.signalCorrArray = sca;
        run.totalCorrArray = tca;

    end
end


%% Get correlations for each speed

windowSize = 0.2./options.binSpacing;

for istate = 1:2
clear stim
    for istim = 1:nStim

    % filter units for stim mean firing rate
%     for iunit =1:numel(units)
%     units(iunit).dmFR = mean(cellfun(@(x) mean(x,'all'), units(iunit).allSpikes{istim,:}),1).*(1/options.binSpacing);
%     end

    tca = nan(numel(goodUnits), numel(goodUnits), nInts-(windowSize-1));

    for iint = 1:nInts-(windowSize-1)

        % pre-process data
        spikeCountArray = nan(nTrials,numel(goodUnits));

        for iunit=1:numel(goodUnits)
            tempSpikeCounts = cellfun(@(x) sum(x(iint:iint+(windowSize-1),:),1), goodUnits(iunit).allSpikes(istim,istate), 'UniformOutput', false);
            %spikeCountArray(:,iunit) = [tempSpikeCounts{:}]';
            spikeCountArray(:, iunit) = cell2mat(tempSpikeCounts)';
        end

        % get total correlation
        totalCorrArray = corr(spikeCountArray,'type','Pearson');

        tca(:,:,iint) = totalCorrArray;

    end

    stim(istim).tca = tca;

    end

    if istate == 1
        stat.stim = stim;
    elseif istate == 2
        run.stim = stim;
    end


end




%% save data

session.stat = stat;
session.run = run;
session.gU = goodUnits;

save(fullfile(dataDir,outputFileName),'session')


end