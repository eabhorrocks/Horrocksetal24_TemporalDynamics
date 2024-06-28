function processSession_decoding_multPopSizes(inputFileName,outputFileName,dataDir)

%% Decoding - training and testing in the same window

%% load data

wheel=struct; % to prevent matlab func issue

load(fullfile(dataDir,inputFileName))


%% pre-process data - classify trials by behavioural state

% remove any trials that didnt run for correcgt amount of time
trialDur = 1;
tolerance = 0.05; % within 50ms
tsd = trials.Speed2D;
trialdurs = [tsd.PDend]-[tsd.PDstart];
invalidDurs_idx = abs(trialdurs-1)>tolerance;
tsd(invalidDurs_idx)=[];

tsd = tsd([tsd.Contrast1]==1 & [tsd.numDots1]==573); % only using full contrast trials


stateTrialType = 'normal'; % 'normal', 'strict', 'changepoints'

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



%% get binned spikes for all trials
dmtrials = tsd;

for itrial =1 :numel(dmtrials)
    dmtrials(itrial).start_time = dmtrials(itrial).PDstart;
    dmtrials(itrial).absVel = abs(dmtrials(itrial).VelX1);
end

for iunit = 1:numel(units)
    units(iunit).spiketimes = units(iunit).spike_times;
end

options.intervalStart = -0.2;
options.intervalEnd = 1.8;
options.binSpacing=0.01; % 10ms bin spacing.

[anUnits, cond] = getBinnedSpikeCounts(dmtrials, units, {'absVel', 'runFlag'}, options);

for iunit = 1:numel(units)
    units(iunit).allSpikes = anUnits(iunit).allSpikes;
end


%% downsample to min # of available trials
minTrial = min(cellfun(@(x) size(x,2), units(1).allSpikes(:)));

for iunit = 1:numel(units)
    units(iunit).allSpikes = cellfun(@(x) x(:,1:minTrial), units(iunit).allSpikes, 'UniformOutput', false);
end


%% get 'good' units

for iunit = 1:numel(units)
    allSpikes = vertcat(units(iunit).allSpikes{:});
    units(iunit).dm_fr = mean(allSpikes(:)).*100;
end

goodUnits = units([units.isi_viol]<=0.1...
    & [units.amplitude_cutoff]<=0.1 & [units.amplitude]>=50 & [units.dm_fr]>=1);


%% get unit spike counts for each trial and time bin
% formatting spike counts for use with decoder functions
clear stat run
stat.D=struct;
run.D=struct;
ii=0;
for ispeed = 1:6
    for itrial = 1:size(goodUnits(1).allSpikes{1},2)
        ii = ii+1;
        for iunit = 1:numel(goodUnits)
            stat.D(ii).data(:,iunit) = goodUnits(iunit).allSpikes{ispeed,1}(:,itrial);
            stat.D(ii).condition = num2str(ispeed);

            run.D(ii).data(:,iunit) = goodUnits(iunit).allSpikes{ispeed,2}(:,itrial);
            run.D(ii).condition = num2str(ispeed);
        end
    end
end

% re-arrange data into cond(ispeed).data(time,unit,trial)
for ispeed = 1:6
    idx = find(strcmp({stat.D.condition},num2str(ispeed)));
    stat.cond(ispeed).catData_sc = cat(3,stat.D(idx).data);

    idx = find(strcmp({run.D.condition},num2str(ispeed)));
    run.cond(ispeed).catData_sc = cat(3,run.D(idx).data);
end


%% decoding coarser time periods of stimulus w/ diff pop sizes

options.kfold = 3;
options.DiscrimType='linear';
nPerms = 3;
options.OptimizeHyperparameters = {'delta','gamma'};
options.Gamma=[];
options.useParallel = false;
nWorkers = 18;
% windowSize=10;

[nTime, nUnits, nTrials] = size(stat.cond(1).catData_sc);
tic

% pop sizes
popSizeVector = [10 20 40 80 120];
popSizeVector(popSizeVector>nUnits)=[];
% popSizeVector(end+1)=nUnits;
timebins = [21:30;...
    31:40;...
    41:50;...
    51:60;...
    61:70;...
    71:80;...
    81:90;...
    91:100;...
    101:110;...
    111:120;];

tic
for ipop = 1:numel(popSizeVector)
    ipop
    thisPopSize = popSizeVector(ipop);
    nReps = ceil((floor(nUnits/thisPopSize)*5)./nWorkers)*nWorkers;
    clear rep

    parfor irep = 1:nReps

        units2use = randperm(nUnits,thisPopSize);
        dataStruct = [];


        for iperm = 1:5
            iperm

            trainidx = randperm(nTrials,nTrain);
            testidx = find(~ismember(1:nTrials,trainidx));%             iperm



            %%% stat %%%
            %             disp('stat')
            for iint = 1:size(timebins,1)
                bin2use = timebins(iint,:);
                for ispeed = 1:6
                    dataStruct(ispeed).data = stat.cond(ispeed).catData_sc(bin2use,units2use,:); %
                end
                [rep(irep).stat.perm(iperm).meanPerf(iint), ~,~,...
                    ~] = do_cvDA(dataStruct,options);
            end

            %%% run %%%
            %             disp('run')
            for iint = 1:size(timebins,1)
                bin2use = timebins(iint,:);
                dataStruct = [];
                for ispeed = 1:6
                    dataStruct(ispeed).data = run.cond(ispeed).catData_sc(bin2use,units2use,:); %
                end
                [rep(irep).run.perm(iperm).meanPerf(iint), ~,~,...
                    ~] = do_cvDA(dataStruct,options);
            end

            stat_shufDataStruct=struct;
            run_shufDataStruct=struct;

            for ispeed = 1:6
                shufOrder = [];
                for iunit = 1:nUnits
                    shufOrder(iunit,:) = randperm(nTrials); % shuffle order for each unit for this speed
                end

                %%% stat %%%
                dataStruct(ispeed).data = stat.cond(ispeed).catData_sc(:,:,:); %
                for iunit = 1:nUnits
                    stat_shufDataStruct(ispeed).data(:,iunit,:) = dataStruct(ispeed).data(:,iunit,shufOrder(iunit,:));
                end

                %%% run %%%
                dataStruct(ispeed).data = run.cond(ispeed).catData_sc(:,:,:); %
                for iunit = 1:nUnits
                    run_shufDataStruct(ispeed).data(:,iunit,:) = dataStruct(ispeed).data(:,iunit,shufOrder(iunit,:));
                end

            end

            %%% do the shuffled-trial decoding %%%
            %             disp('stat shuf')
            for iint = 1:size(timebins,1)
                bin2use = timebins(iint,:);
                tempDataStruct=struct;
                for ispeed = 1:6
                    tempDataStruct(ispeed).data = stat_shufDataStruct(ispeed).data(bin2use,units2use,:); %
                end

                [rep(irep).statShuf.perm(iperm).meanPerf(iint), ~,~,...
                    ~]  = do_cvDA(tempDataStruct,options);
                %                 clear tempDataStruct
            end

            %             disp('run shuf')
            for iint = 1:size(timebins,1)
                bin2use = timebins(iint,:);
                tempDataStruct=struct;

                for ispeed = 1:6
                    tempDataStruct(ispeed).data = run_shufDataStruct(ispeed).data(bin2use,units2use,:); %
                end

                [rep(irep).runShuf.perm(iperm).meanPerf(iint), ~,~,...
                    ~]  = do_cvDA(tempDataStruct,options);
                %                 clear tempDataStruct

            end

        end


    end

    popSize(ipop).rep = rep;
    popSize(ipop).nUnits = thisPopSize;

end
toc


%% save data

session.stat = stat;
session.run = run;
session.goodUnits = goodUnits;
session.popSize = popSize;

dummyVar=[];

save(fullfile(dataDir,outputFileName),'dummyVar', 'session', '-v7.3')

end








