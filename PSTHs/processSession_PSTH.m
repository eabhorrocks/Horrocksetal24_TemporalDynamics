function processSession_PSTH(inputFileName,outputFileName,dataDir)


%% load the **_basic.mat file to process
load(fullfile(dataDir,inputFileName))


%% process trials
% filter trials to the ones we want to analyse and sort according to
% behavioural state

% exclude any trials that did not run for appropriate time (sometime 1st
% trial)
trialDur = 1;
tolerance = 0.05; % within 50ms
tsd = trials.Speed2D;
trialdurs = [tsd.PDend]-[tsd.PDstart];
invalidDurs_idx = abs(trialdurs-1)>tolerance;
tsd(invalidDurs_idx)=[];

% only using full contrast trials for this analysis
tsd = tsd([tsd.Contrast1]==1);

%% split trials according to behavioural state

stateTrialType = 'FastSlowLoco'; % 'normal', 'strict', 'changepoints'

switch stateTrialType
    case 'normal'

        % normal paper criteria
        statTrials = tsd(cellfun(@(x) prop(x<3)>=0.75 & mean(x)<0.5, {tsd.WheelSpeed}));
        runTrials = tsd(cellfun(@(x) prop(x>0.5)>=0.75 & mean(x)>3, {tsd.WheelSpeed}));

    case 'strict'
        % stricter criteria
        statTrials = tsd(cellfun(@(x) prop(x<0.5)>=1 & mean(x)<0.5, {tsd.WheelSpeed}));
        runTrials = tsd(cellfun(@(x) prop(x>0.5)>=0.9 & mean(x)>3, {tsd.WheelSpeed}));

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
        runTrials = tsd(runIdx);
        % use stricter stationary trial criteria
        statTrials = tsd(cellfun(@(x) prop(x<0.5)>=1 & mean(x)<0.5, {tsd.WheelSpeed}));


    case 'FastSlowLoco' % split run ttrials according to quantiled run speed
        tempTrials = tsd(cellfun(@(x) prop(x>0.5)>=0.75 & mean(x)>3, {tsd.WheelSpeed}));
        tempTrials = tempTrials([tempTrials.MeanStimRunSpeed]>3);

            q = quantile([tempTrials.MeanStimRunSpeed],3); % split into quartiles, discard middle
            statTrials = tempTrials([tempTrials.MeanStimRunSpeed]<=q(1)); % fake 'stat' name
            runTrials = tempTrials([tempTrials.MeanStimRunSpeed]>=q(3));


end


%% optionally detect and remove any trials w/ clear eye movements
removeEyeMovementTrials = false;
if removeEyeMovementTrials

    for itrial = 1:numel(statTrials)
        nTimePoints = numel(statTrials(itrial).EyeX);
        dist=[];
        for itime = 1:nTimePoints-1
            dx = statTrials(itrial).EyeX(itime+1)-statTrials(itrial).EyeX(itime);
            dy = statTrials(itrial).EyeY(itime+1)-statTrials(itrial).EyeY(itime);
            dist(itime) = sqrt(dx^2+dy^2);
        end

        statTrials(itrial).meanEyeDist = mean(dist,'omitnan');
        statTrials(itrial).maxEyeDist = max(dist);

    end

    for itrial = 1:size(runTrials)
        nTimePoints = numel(runTrials(itrial).EyeX);
        dist=[];
        for itime = 1:nTimePoints-1
            dx = runTrials(itrial).EyeX(itime+1)-runTrials(itrial).EyeX(itime);
            dy = runTrials(itrial).EyeY(itime+1)-runTrials(itrial).EyeY(itime);
            dist(itime) = sqrt(dx^2+dy^2);
        end

        runTrials(itrial).meanEyeDist = mean(dist,'omitnan');
        runTrials(itrial).maxEyeDist = max(dist);
    end


    allEyeDists = cat(2,[statTrials.maxEyeDist], [runTrials.maxEyeDist]);

    [TF,L,U,C] = isoutlier(allEyeDists,'median','ThresholdFactor',2.5);

    % U is the upper limit for max eye movement.
    statTrials = statTrials([statTrials.maxEyeDist]<=U);
    runTrials = runTrials([runTrials.maxEyeDist]<=U);

end


%% split stationary trials according to pupil data - bottom/top tertiles
quants = quantile([statTrials.MeanStimEyeArea],2); % exclude middle 1/3
statSmallTrials = statTrials([statTrials.MeanStimEyeArea]<quants(1));
statBigTrials = statTrials([statTrials.MeanStimEyeArea]>quants(2));

% blank trials where no stim shown
statTrials_blank = statTrials([statTrials.numDots1]==0);
runTrials_blank = runTrials([runTrials.numDots1]==0);
statSmallTrials_blank = statSmallTrials([statSmallTrials.numDots1]==0);
statBigTrials_blank = statBigTrials([statBigTrials.numDots1]==0);

% set start_time field as PDstart value for blank trials
for itrial = 1:numel(statTrials_blank)
    statTrials_blank(itrial).start_time = statTrials_blank(itrial).PDstart;
end
for itrial = 1:numel(runTrials_blank)
    runTrials_blank(itrial).start_time = runTrials_blank(itrial).PDstart;
end


%% get trial intervals for PSTH function
% this finds the appropriate time intervals when different stimuli were
% shown and formats for use with the 'makePSTH' function.

speeds = [0, -16, -32, -64, -128, -256];

for ispeed = 1:numel(speeds)
    % stationary/locomotion trials
    stimTrials = statTrials([statTrials.VelX1]==speeds(ispeed));
    stat_intervalsCell{ispeed}= [vertcat(stimTrials.PDstart), vertcat(stimTrials.PDend)].*1000;

    stimTrials = runTrials([runTrials.VelX1]==speeds(ispeed));
    run_intervalsCell{ispeed}= [vertcat(stimTrials.PDstart), vertcat(stimTrials.PDend)].*1000;

    % stationary trials with small and large pupil size
    stimTrials = statSmallTrials([statSmallTrials.VelX1]==speeds(ispeed));
    statSmall_intervalsCell{ispeed}= [vertcat(stimTrials.PDstart), vertcat(stimTrials.PDend)].*1000;

    stimTrials = statBigTrials([statBigTrials.VelX1]==speeds(ispeed));
    statBig_intervalsCell{ispeed}= [vertcat(stimTrials.PDstart), vertcat(stimTrials.PDend)].*1000;

end

stat_blankIntervals = [vertcat(statTrials_blank.PDstart), vertcat(statTrials_blank.PDend)].*1000;
run_blankIntervals = [vertcat(runTrials_blank.PDstart), vertcat(runTrials_blank.PDend)].*1000;

statSmall_blankIntervals = [vertcat(statSmallTrials_blank.PDstart), vertcat(statSmallTrials_blank.PDend)].*1000;
statBig_blankIntervals = [vertcat(statBigTrials_blank.PDstart), vertcat(statBigTrials_blank.PDend)].*1000;


%% downsample nTrials
% downsample trials so that each condition has equal trial counts
% taking 1:minTrial (repeatable).
minTrial = min([min(cellfun(@(x) size(x,1), stat_intervalsCell)),...
    min(cellfun(@(x) size(x,1), run_intervalsCell))]);


%stat_intervalsCell = cellfun(@(x) x(1:minTrial,:), stat_intervalsCell,'UniformOutput',false);
%run_intervalsCell = cellfun(@(x) x(1:minTrial,:), run_intervalsCell,'UniformOutput',false);


%% generate PSTHs

% options for PSTH function
options.binWidth = 10;
options.preTime = 200;
options.postTime = 800;
options.smoothType = 'gaussian';
options.smoothWidth = 175;

options.getReliability = true;
options.kfold = 3; % set to number>nTrials to do Leave-one-out.
options.distType = 'correlation';
options.nReliPerms = 25;
options.nShuffle = 40;

options.plot = false;
options.cols = inferno(7); options.cols(1,:) = [];
options.plotFRrange = [];
options.plot95CI = false;

for iunit =1:numel(units)
    %[psth, blank, bmp] = makePSTH(spikeTimes,intervalsCell,blankIntervals,options);
    [units(iunit).stat_psth, units(iunit).stat_blank, ~]            = makePSTH(units(iunit).spike_times*1000,stat_intervalsCell,stat_blankIntervals,options);
    [units(iunit).run_psth, units(iunit).run_blank, bmp]            = makePSTH(units(iunit).spike_times*1000,run_intervalsCell,run_blankIntervals,options);

 %   [units(iunit).statSmall_psth, units(iunit).statSmall_blank, ~]  = makePSTH(units(iunit).spike_times*1000,statSmall_intervalsCell,statSmall_blankIntervals,options);
 %   [units(iunit).statBig_psth, units(iunit).statBig_blank, ~]      = makePSTH(units(iunit).spike_times*1000,statBig_intervalsCell,statBig_blankIntervals,options);
end


%% get psth metrics

% input args for PSTH metrics
zThresh = 3.29; % 99.9% CI
minRespEnd = 1000; % set earliest response 'end' to stimulus offset (for sustainedness index metric)

for iunit = 1:numel(units)

    %%% stat psth %%%
    stat_baseline = mean(units(iunit).stat_blank.psth);
    clear metrics
    % generate stuct array of psth metrics
    for ispeed = 1:6
        metrics(ispeed) = calc_PSTHmetrics(units(iunit).stat_psth(ispeed).psth,...
            units(iunit).stat_psth(ispeed).zpsth, bmp, stat_baseline, zThresh, minRespEnd);
    end
    % copy metrics fields to psth structs
    units(iunit).stat_psth = copyStructFields(metrics,units(iunit).stat_psth);


    %%% run psth %%%
    run_baseline = mean(units(iunit).run_blank.psth);
    clear metrics
    % generate stuct array of psth metrics
    for ispeed = 1:6
        metrics(ispeed) = calc_PSTHmetrics(units(iunit).run_psth(ispeed).psth,...
            units(iunit).run_psth(ispeed).zpsth, bmp, run_baseline, zThresh, minRespEnd);
    end

    % copy metrics fields to psth structs
    units(iunit).run_psth = copyStructFields(metrics,units(iunit).run_psth);

    %%% stat small pupil psth %%%
%     statSmall_baseline = mean(units(iunit).statSmall_blank.psth);
%     clear metrics
%     % generate stuct array of psth metrics
%     for ispeed = 1:6
%         metrics(ispeed) = calc_PSTHmetrics(units(iunit).statSmall_psth(ispeed).psth,...
%             units(iunit).statSmall_psth(ispeed).zpsth, bmp, statSmall_baseline, zThresh, minRespEnd);
%     end
%     % copy metrics fields to psth structs
%     units(iunit).statSmall_psth = copyStructFields(metrics,units(iunit).statSmall_psth);
% 
%     %%% stat big pupil psth %%%
%     statBig_baseline = mean(units(iunit).statBig_blank.psth);
%     clear metrics
%     % generate stuct array of psth metrics
%     for ispeed = 1:6
%         metrics(ispeed) = calc_PSTHmetrics(units(iunit).statBig_psth(ispeed).psth,...
%             units(iunit).statBig_psth(ispeed).zpsth, bmp, statBig_baseline, zThresh, minRespEnd);
%     end
%     % copy metrics fields to psth structs
%     units(iunit).statBig_psth = copyStructFields(metrics,units(iunit).statBig_psth);

end



%% fit descriptive functions to PSTHs
% fit descriptive functions to classify stimulus onset and offset features

rangeThresh = 3;
reliThresh =-1.645;
zThresh = 3.29;

for iunit= 1:numel(units)

    for ispeed = 1:6

        % stat
        response = normalize(units(iunit).stat_psth(ispeed).psth,'range');
        resp_onset = response(21:51); % t=0 to t=300ms
        resp_offset = response(121:151); % t=1000ms to t=1300ms

        [units(iunit).stat_psth(ispeed).onset_bestParams,...
            units(iunit).stat_psth(ispeed).onset_char,...
            units(iunit).stat_psth(ispeed).onset_bestR2] = fitGaussianTemplates_meanPSTH(resp_onset,5, false);

        [units(iunit).stat_psth(ispeed).offset_bestParams,...
            units(iunit).stat_psth(ispeed).offset_char,...
            units(iunit).stat_psth(ispeed).offset_bestR2] = fitGaussianTemplates_meanPSTH(resp_offset,5, false);

        % run
        response = normalize(units(iunit).run_psth(ispeed).psth,'range');
        resp_onset = response(21:51);
        resp_offset = response(121:151);

        [units(iunit).run_psth(ispeed).onset_bestParams,...
            units(iunit).run_psth(ispeed).onset_char,...
            units(iunit).run_psth(ispeed).onset_bestR2] = fitGaussianTemplates_meanPSTH(resp_onset,5, false);

        [units(iunit).run_psth(ispeed).offset_bestParams,...
            units(iunit).run_psth(ispeed).offset_char,...
            units(iunit).run_psth(ispeed).offset_bestR2] = fitGaussianTemplates_meanPSTH(resp_offset,5, false);


        % stat small
%         response = normalize(units(iunit).statSmall_psth(ispeed).psth,'range');
%         resp_onset = response(21:51);
%         resp_offset = response(121:151);
% 
%         [units(iunit).statSmall_psth(ispeed).onset_bestParams,...
%             units(iunit).statSmall_psth(ispeed).onset_char,...
%             units(iunit).statSmall_psth(ispeed).onset_bestR2] = fitGaussianTemplates_meanPSTH(resp_onset,5, false);
% 
%         [units(iunit).statSmall_psth(ispeed).offset_bestParams,...
%             units(iunit).statSmall_psth(ispeed).offset_char,...
%             units(iunit).statSmall_psth(ispeed).offset_bestR2] = fitGaussianTemplates_meanPSTH(resp_offset,5, false);
% 
% 
%         % stat big
%         response = normalize(units(iunit).statBig_psth(ispeed).psth,'range');
%         resp_onset = response(21:51);
%         resp_offset = response(121:151);
% 
%         [units(iunit).statBig_psth(ispeed).onset_bestParams,...
%             units(iunit).statBig_psth(ispeed).onset_char,...
%             units(iunit).statBig_psth(ispeed).onset_bestR2] = fitGaussianTemplates_meanPSTH(resp_onset,5, false);
% 
%         [units(iunit).statBig_psth(ispeed).offset_bestParams,...
%             units(iunit).statBig_psth(ispeed).offset_char,...
%             units(iunit).statBig_psth(ispeed).offset_bestR2] = fitGaussianTemplates_meanPSTH(resp_offset,5, false);

    end
end






%% save data
%saveName = [sessionTags{isession,1},'_', sessionTags{isession,2},'_psth3rds.mat'];
save(fullfile(dataDir,outputFileName),'trials', 'wheel', 'units')

end


