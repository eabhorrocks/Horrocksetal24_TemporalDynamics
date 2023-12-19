function processSession_PSTH(inputFileName,outputFileName,dataDir)


load(fullfile(dataDir,inputFileName))



%% process trials
% fitler trials to the ones we want to analyse and sort according to
% behavioural state

% exclude any trials that did not run for appropriate time (sometime 1st
% trial)
trialDur = 1;
tolerance = 0.05; % within 50ms
tsd = trials.Speed2D;
trialdurs = [tsd.PDend]-[tsd.PDstart];
invalidDurs_idx = abs(trialdurs-1)>tolerance;
tsd(invalidDurs_idx)=[];

tsd = tsd([tsd.Contrast1]==1); % only using full contrast trials

% split trials by state according to wheel data
runTrials = tsd(cellfun(@(x) prop(x>0.5)>=0.75 & mean(x)>3, {tsd.WheelSpeed}));
statTrials = tsd(cellfun(@(x) prop(x<3)>=0.75 & mean(x)<0.5, {tsd.WheelSpeed}));

% split stat trials according to pupil data - bottom/top tertiles
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

stat_intervalsCell = cellfun(@(x) x(1:minTrial,:), stat_intervalsCell,'UniformOutput',false);
run_intervalsCell = cellfun(@(x) x(1:minTrial,:), run_intervalsCell,'UniformOutput',false);


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

    [units(iunit).statSmall_psth, units(iunit).statSmall_blank, ~]  = makePSTH(units(iunit).spike_times*1000,statSmall_intervalsCell,statSmall_blankIntervals,options);
    [units(iunit).statBig_psth, units(iunit).statBig_blank, ~]      = makePSTH(units(iunit).spike_times*1000,statBig_intervalsCell,statBig_blankIntervals,options);
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
    statSmall_baseline = mean(units(iunit).statSmall_blank.psth);
    clear metrics
    % generate stuct array of psth metrics
    for ispeed = 1:6
        metrics(ispeed) = calc_PSTHmetrics(units(iunit).statSmall_psth(ispeed).psth,...
            units(iunit).statSmall_psth(ispeed).zpsth, bmp, statSmall_baseline, zThresh, minRespEnd);
    end
    % copy metrics fields to psth structs
    units(iunit).statSmall_psth = copyStructFields(metrics,units(iunit).statSmall_psth);

    %%% stat big pupil psth %%%
    statBig_baseline = mean(units(iunit).statBig_blank.psth);
    clear metrics
    % generate stuct array of psth metrics
    for ispeed = 1:6
        metrics(ispeed) = calc_PSTHmetrics(units(iunit).statBig_psth(ispeed).psth,...
            units(iunit).statBig_psth(ispeed).zpsth, bmp, statBig_baseline, zThresh, minRespEnd);
    end
    % copy metrics fields to psth structs
    units(iunit).statBig_psth = copyStructFields(metrics,units(iunit).statBig_psth);

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
        response = normalize(units(iunit).statSmall_psth(ispeed).psth,'range');
        resp_onset = response(21:51);
        resp_offset = response(121:151);

        [units(iunit).statSmall_psth(ispeed).onset_bestParams,...
            units(iunit).statSmall_psth(ispeed).onset_char,...
            units(iunit).statSmall_psth(ispeed).onset_bestR2] = fitGaussianTemplates_meanPSTH(resp_onset,5, false);

        [units(iunit).statSmall_psth(ispeed).offset_bestParams,...
            units(iunit).statSmall_psth(ispeed).offset_char,...
            units(iunit).statSmall_psth(ispeed).offset_bestR2] = fitGaussianTemplates_meanPSTH(resp_offset,5, false);


        % stat big
        response = normalize(units(iunit).statBig_psth(ispeed).psth,'range');
        resp_onset = response(21:51);
        resp_offset = response(121:151);

        [units(iunit).statBig_psth(ispeed).onset_bestParams,...
            units(iunit).statBig_psth(ispeed).onset_char,...
            units(iunit).statBig_psth(ispeed).onset_bestR2] = fitGaussianTemplates_meanPSTH(resp_onset,5, false);

        [units(iunit).statBig_psth(ispeed).offset_bestParams,...
            units(iunit).statBig_psth(ispeed).offset_char,...
            units(iunit).statBig_psth(ispeed).offset_bestR2] = fitGaussianTemplates_meanPSTH(resp_offset,5, false);

    end
end






%% save data
%saveName = [sessionTags{isession,1},'_', sessionTags{isession,2},'_psth3rds.mat'];
save(fullfile(dataDir,outputFileName),'trials', 'wheel', 'units')

end


