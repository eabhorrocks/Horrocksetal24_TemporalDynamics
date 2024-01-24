function processSession_decoding(inputFileName,outputFileName,dataDir)

%% Decoding - training and testing in the same window

%% load data

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

%% Decoding options

options.kfold = 3;
options.DiscrimType='linear';
nPerms = 5;
windowSize=10;

[~, nUnits, nTrials] = size(stat.cond(1).catData_sc);
tic

%% decoding with gamma optimised for each cross-val
% disp('decoding with gamma optimised')
% %params = hyperparameters('fitcdiscr', nan,nan);
% options.OptimizeHyperparameters = {'delta', 'gamma'};
% options.Gamma=[];
%
% for iperm = 1:nPerms
%     %iperm
%     %%% stat %%%
%
%     for iint = 1:200-(windowSize-1)
%         bin2use = iint:(iint+(windowSize-1));
%         dataStruct = [];
%         for ispeed = 1:6
%             dataStruct(ispeed).data = stat.cond(ispeed).catData_sc(bin2use,:,:); %
%         end
%         [stat.GammaOpt.perm(iperm).meanPerf(iint), ~,~,...
%             ~] = do_cvDA(dataStruct,options);
%     end
%
%     %%% run %%%
%
%     for iint = 1:200-(windowSize-1)
%         bin2use = iint:(iint+(windowSize-1));
%         dataStruct = [];
%         for ispeed = 1:6
%             dataStruct(ispeed).data = run.cond(ispeed).catData_sc(bin2use,:,:); %
%         end
%         [run.GammaOpt.perm(iperm).meanPerf(iint), ~,~,...
%             ~] = do_cvDA(dataStruct,options);
%     end
%
% end

%% decoding with gamma=1 (diagonal covariance matrix)
disp('decoding with gamma = 1')

options.OptimizeHyperparameters =  {'delta'};
options.Gamma=1;

for iperm = 1:nPerms
    %iperm
    %%% stat %%%

    for iint = 1:200-(windowSize-1)
        bin2use = iint:(iint+(windowSize-1));
        dataStruct = [];
        for ispeed = 1:6
            dataStruct(ispeed).data = stat.cond(ispeed).catData_sc(bin2use,:,:); %
        end
        [stat.GammaOne.perm(iperm).meanPerf(iint), ~,~,...
            ~] = do_cvDA(dataStruct,options);
    end

    %%% run %%%

    for iint = 1:200-(windowSize-1)
        bin2use = iint:(iint+(windowSize-1));
        dataStruct = [];
        for ispeed = 1:6
            dataStruct(ispeed).data = run.cond(ispeed).catData_sc(bin2use,:,:); %
        end
        [run.GammaOne.perm(iperm).meanPerf(iint), ~,~,...
            ~] = do_cvDA(dataStruct,options);
    end

end

%% decoding with gamma=0 (normal covariance matrix)
% disp('decoding with gamma = 0')
%
% options.OptimizeHyperparameters =  {'delta'};
% options.Gamma=0;
%
% for iperm = 1:nPerms
%     %iperm
%     %%% stat %%%
%
%     for iint = 1:200-(windowSize-1)
%         bin2use = iint:(iint+(windowSize-1));
%         dataStruct = [];
%         for ispeed = 1:6
%             dataStruct(ispeed).data = stat.cond(ispeed).catData_sc(bin2use,:,:); %
%         end
%         [stat.GammaZero.perm(iperm).meanPerf(iint), ~,~,...
%             ~] = do_cvDA(dataStruct,options);
%     end
%
%     %%% run %%%
%
%     for iint = 1:200-(windowSize-1)
%         bin2use = iint:(iint+(windowSize-1));
%         dataStruct = [];
%         for ispeed = 1:6
%             dataStruct(ispeed).data = run.cond(ispeed).catData_sc(bin2use,:,:); %
%         end
%         [run.GammaZero.perm(iperm).meanPerf(iint), ~,~,...
%             ~] = do_cvDA(dataStruct,options);
%     end
%
% end

%% Doing decoding on shuffled trials to disrupt noise correlations


%% decoding shuffled trials with optimised gamma
% disp('decoding with shuffled trials and gamma optimised')
%
%
% options.OptimizeHyperparameters = {'delta', 'gamma'};
% options.Gamma=[];
% %%% shuffle trials within stimulus conditions %%%
%
% for iperm =1:nPerms
%
%     % shuffle trials for each speed/unit
%
%     for ispeed = 1:6
%         shufOrder = [];
%         for iunit = 1:nUnits
%             shufOrder(iunit,:) = randperm(nTrials); % shuffle order for each unit for this speed
%         end
%
%         %%% stat %%%
%         dataStruct(ispeed).data = stat.cond(ispeed).catData_sc(:,:,:); %
%         for iunit = 1:nUnits
%             stat_shufDataStruct(ispeed).data(:,iunit,:) = dataStruct(ispeed).data(:,iunit,shufOrder(iunit,:));
%         end
%
%         %%% run %%%
%         dataStruct(ispeed).data = run.cond(ispeed).catData_sc(:,:,:); %
%         for iunit = 1:nUnits
%             run_shufDataStruct(ispeed).data(:,iunit,:) = dataStruct(ispeed).data(:,iunit,shufOrder(iunit,:));
%         end
%
%     end
%
%
%
%     %%% do the shuffled-trial decoding %%%
%     for iint = 1:200-(windowSize-1)
%         bin2use = iint:(iint+(windowSize-1));
%         for ispeed = 1:6
%             tempDataStruct(ispeed).data = stat_shufDataStruct(ispeed).data(bin2use,:,:); %
%         end
%
%         [stat.shuffleGammaOpt.perm(iperm).meanPerf(iint), ~,~,...
%             ~] = do_cvDA(tempDataStruct,options);
%         clear tempDataStruct
%     end
%
%
%     for iint = 1:200-(windowSize-1)
%         bin2use = iint:(iint+(windowSize-1));
%         for ispeed = 1:6
%             tempDataStruct(ispeed).data = run_shufDataStruct(ispeed).data(bin2use,:,:); %
%         end
%
%         [run.shuffleGammaOpt.perm(iperm).meanPerf(iint), ~,~,...
%             ~] = do_cvDA(tempDataStruct,options);
%         clear tempDataStruct
%     end
% end


%% decoding shuffled trials with gamma=0

% disp('decoding with shuffled trials and gamma=0')
%
% options.OptimizeHyperparameters = {'delta'};
% options.Gamma=0;
% %%% shuffle trials within stimulus conditions %%%
%
% for iperm =1:nPerms
%
%     % shuffle trials for each speed/unit
%
%     for ispeed = 1:6
%         shufOrder = [];
%         for iunit = 1:nUnits
%             shufOrder(iunit,:) = randperm(nTrials); % shuffle order for each unit for this speed
%         end
%
%         %%% stat %%%
%         dataStruct(ispeed).data = stat.cond(ispeed).catData_sc(:,:,:); %
%         for iunit = 1:nUnits
%             stat_shufDataStruct(ispeed).data(:,iunit,:) = dataStruct(ispeed).data(:,iunit,shufOrder(iunit,:));
%         end
%
%         %%% run %%%
%         dataStruct(ispeed).data = run.cond(ispeed).catData_sc(:,:,:); %
%         for iunit = 1:nUnits
%             run_shufDataStruct(ispeed).data(:,iunit,:) = dataStruct(ispeed).data(:,iunit,shufOrder(iunit,:));
%         end
%
%     end
%
%
%
%     %%% do the shuffled-trial decoding %%%
%     for iint = 1:200-(windowSize-1)
%         bin2use = iint:(iint+(windowSize-1));
%         for ispeed = 1:6
%             tempDataStruct(ispeed).data = stat_shufDataStruct(ispeed).data(bin2use,:,:); %
%         end
%
%         [stat.shuffleGammaOne.perm(iperm).meanPerf(iint), ~,~,...
%             ~] = do_cvDA(tempDataStruct,options);
%         clear tempDataStruct
%     end
%
%
%     for iint = 1:200-(windowSize-1)
%         bin2use = iint:(iint+(windowSize-1));
%         for ispeed = 1:6
%             tempDataStruct(ispeed).data = run_shufDataStruct(ispeed).data(bin2use,:,:); %
%         end
%
%         [run.shuffleGammaOne.perm(iperm).meanPerf(iint), ~,~,...
%             ~] = do_cvDA(tempDataStruct,options);
%         clear tempDataStruct
%     end
% end
%

%% decoding with quadratic discriminant analysis (covariance varies between stimuli)

% disp('decoding with QDA')
%
% options.OptimizeHyperparameters = {'gamma', 'delta'};
% options.Gamma=[];
% options.discrimType = 'quadratic';
%
% for iperm = 1:nPerms
%     %iperm
%     %%% stat %%%
%
%     for iint = 1%:200-(windowSize-1)
%         bin2use = iint:(iint+(windowSize-1));
%         dataStruct = [];
%         for ispeed = 1:6
%             dataStruct(ispeed).data = stat.cond(ispeed).catData_sc(bin2use,:,:); %
%         end
%         [stat.quadratic.perm(iperm).meanPerf(iint), ~,~,...
%             ~] = do_cvDA(dataStruct,options);
%     end
%
%     %%% run %%%
%
%     for iint = 1:200-(windowSize-1)
%         bin2use = iint:(iint+(windowSize-1));
%         dataStruct = [];
%         for ispeed = 1:6
%             dataStruct(ispeed).data = run.cond(ispeed).catData_sc(bin2use,:,:); %
%         end
%         [run.quadratic.perm(iperm).meanPerf(iint), ~,~,...
%             ~] = do_cvDA(dataStruct,options);
%     end
%
% end

%% decoding with quadratic discriminant analysis (covariance varies between stimuli)

% disp('decoding with diagQDA')
%
% options.OptimizeHyperparameters = {'delta'}
% options.Gamma=[];
% options.discrimType = 'diagquadratic';
%
% for iperm = 1:nPerms
%     %iperm
%     %%% stat %%%
%
%     for iint = 1:200-(windowSize-1)
%         bin2use = iint:(iint+(windowSize-1));
%         dataStruct = [];
%         for ispeed = 1:6
%             dataStruct(ispeed).data = stat.cond(ispeed).catData_sc(bin2use,:,:); %
%         end
%         [stat.diagquadratic.perm(iperm).meanPerf(iint), ~,~,...
%            ~] = do_cvDA(dataStruct,options);
%     end
%
%     %%% run %%%
%
%     for iint = 1:200-(windowSize-1)
%         bin2use = iint:(iint+(windowSize-1));
%         dataStruct = [];
%         for ispeed = 1:6
%             dataStruct(ispeed).data = run.cond(ispeed).catData_sc(bin2use,:,:); %
%         end
%         [run.diagquadratic.perm(iperm).meanPerf(iint), ~,~,...
%             ~] = do_cvDA(dataStruct,options);
%     end

% end

%% save data

session.stat = stat;
session.run = run;
session.goodUnits = goodUnits;

dummyVar=[];
save(fullfile(dataDir,outputFileName),'dummyVar', 'session', '-v7.3')

end








