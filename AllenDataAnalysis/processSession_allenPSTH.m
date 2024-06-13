%% analysis of PSTH responses from allen's 'Functional Connectivity' dataset

dataDirectory = 'E:\AllenInstituteAnalyses\AllenDataAnalysis\Data\rawMatFiles';
dataFiles = dir(fullfile(dataDirectory, 'session*'));

tic
for isession = 1:numel(dataFiles)
    isession
    clear units trials runInfo pupilInfo sessionInfo
    sessionFile = fullfile([dataDirectory,'\' dataFiles(isession).name]);


    requiredStims = {'dot_motion'};

    % default options
    options = struct;

    [units, trials, runInfo, pupilInfo, sessionInfo] = preprocessAllenData(sessionFile, requiredStims, options);


    %% basic info about trials and areas of interest
    areas = {'VISp', 'VISl', 'VISal', 'VISrl', 'VISam', 'VISpm', 'LGd', 'LP'};
    dmtrials = trials(strcmp([trials.stimulus_name],'dot_motion'));
    speeds = unique([dmtrials.Speed]);
    alldirs = unique([dmtrials.Dir]);

    %% process trials (stationary vs running)

    % only use trials with durations ~1s
    trialDur = 1;
    tolerance = 0.05; % within 50ms
    trialdurs = [dmtrials.stop_time]-[dmtrials.start_time];
    invalidDurs_idx = abs(trialdurs-1)>tolerance;
    dmtrials(invalidDurs_idx)=[];

    % split trials by state according to wheel data
    runTrials = dmtrials(cellfun(@(x) prop(x>0.5)>=0.75 & mean(x)>3, {dmtrials.runTrace}));
    statTrials = dmtrials(cellfun(@(x) prop(x<3)>=0.75 & mean(x)<0.5, {dmtrials.runTrace}));

    %% PSTH options

    % options
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
    options.cols = inferno(7);
    options.plotFRrange = [];
    options.plot95CI = false;


    %% initialise struct fields for PSTHs

    nUnits = numel(units);

    for iunit = 1:nUnits
        units(iunit).statPSTH = [];
        units(iunit).runPSTH = [];
    end


    %% generate PSTHs
    bmp = -195:10:1795;
    % STAT PSTHS
    for idir = 1:4
        for ispeed = 1:numel(speeds)
            t = statTrials([statTrials.Speed]==speeds(ispeed) & [statTrials.Dir]==alldirs(idir));
            sints{idir,ispeed}= [vertcat(t.start_time), vertcat(t.stop_time)].*1000;
        end
    end
    parfor iunit = 1:nUnits
        if ismember(units(iunit).ecephys_structure_acronym, areas)
            [units(iunit).statPSTH,  ~, bmp] =...
                makePSTH(units(iunit).spiketimes*1000,sints,[],options);
           
        end
    end


    clear sints

    % RUN PSTHS
    for idir = 1:4
        for ispeed = 1:numel(speeds)
            t = runTrials([runTrials.Speed]==speeds(ispeed) & [runTrials.Dir]==alldirs(idir));
            sints{idir,ispeed}= [vertcat(t.start_time), vertcat(t.stop_time)].*1000;
        end
    end
    parfor iunit = 1:nUnits
        if ismember(units(iunit).ecephys_structure_acronym, areas)
            [units(iunit).runPSTH,  ~, bmp] =...
                makePSTH(units(iunit).spiketimes*1000,sints,[],options);
        end
    end

    %% get basic metrics

    % sorry this is an ugly adaptation of previous code
    bmp = -195:10:1795;
    % input args for PSTH metrics
    zThresh = 3.29; % 99.9% CI
    minRespEnd = 1000;
    minTrials = 10;

    for iunit = 1:nUnits
        if ismember(units(iunit).ecephys_structure_acronym, areas)

            %%% stat psth %%%
            clear metrics
            % generate stuct array of psth metrics
            parfor idir = 1:4
                for ispeed = 1:7

                    if units(iunit).statPSTH(idir,ispeed).nTrials>=minTrials

                        metrics(idir,ispeed) = calc_PSTHmetrics(units(iunit).statPSTH(idir,ispeed).psth,...
                            units(iunit).statPSTH(idir,ispeed).zpsth, bmp, mean(units(iunit).statPSTH(idir,ispeed).psth(1:20)), zThresh, minRespEnd);
                    else
                        metrics(idir,ispeed).responsive = nan;
                        metrics(idir,ispeed).excited = nan;
                        metrics(idir,ispeed).suppressed = nan;
                        metrics(idir,ispeed).onsetLatency = nan;
                        metrics(idir,ispeed).offsetLatency = nan;
                        metrics(idir,ispeed).peakLatency = nan;

                        metrics(idir,ispeed).peakSign = nan;
                        metrics(idir,ispeed).peakAbsZ = nan;
                        metrics(idir,ispeed).peak90Latency = nan;
                        metrics(idir,ispeed).meanZ = nan;
                        metrics(idir,ispeed).minZ = nan;
                        metrics(idir,ispeed).maxZ = nan;
                        metrics(idir,ispeed).rangeZ = nan;
                        metrics(idir,ispeed).minZLat = nan;
                        metrics(idir,ispeed).maxZLat = nan;

                        metrics(idir,ispeed).maxFR = nan;
                        metrics(idir,ispeed).minFR = nan;
                        metrics(idir,ispeed).meanFR = nan;
                        metrics(idir,ispeed).rangeFR = nan;
                        metrics(idir,ispeed).susIdx = nan;
                    end
                end
            end


            % copy metrics fields to psth structs
            units(iunit).statPSTH = copyStructFields(metrics,units(iunit).statPSTH);
            clear metrics
        end

        if ismember(units(iunit).ecephys_structure_acronym, areas)

            %%% run psth %%%

            % generate stuct array of psth metrics
            parfor idir = 1:4
                for ispeed = 1:7

                    if units(iunit).runPSTH(idir,ispeed).nTrials>=minTrials

                        metrics(idir,ispeed) = calc_PSTHmetrics(units(iunit).runPSTH(idir,ispeed).psth,...
                            units(iunit).runPSTH(idir,ispeed).zpsth, bmp, mean(units(iunit).runPSTH(idir,ispeed).psth(1:20)), zThresh, minRespEnd);
                    else
                        metrics(idir,ispeed).responsive = nan;
                        metrics(idir,ispeed).excited = nan;
                        metrics(idir,ispeed).suppressed = nan;
                        metrics(idir,ispeed).onsetLatency = nan;
                        metrics(idir,ispeed).offsetLatency = nan;
                        metrics(idir,ispeed).peakLatency = nan;

                        metrics(idir,ispeed).peakSign = nan;
                        metrics(idir,ispeed).peakAbsZ = nan;
                        metrics(idir,ispeed).peak90Latency = nan;
                        metrics(idir,ispeed).meanZ = nan;
                        metrics(idir,ispeed).minZ = nan;
                        metrics(idir,ispeed).maxZ = nan;
                        metrics(idir,ispeed).rangeZ = nan;
                        metrics(idir,ispeed).minZLat = nan;
                        metrics(idir,ispeed).maxZLat = nan;

                        metrics(idir,ispeed).maxFR = nan;
                        metrics(idir,ispeed).minFR = nan;
                        metrics(idir,ispeed).meanFR = nan;
                        metrics(idir,ispeed).rangeFR = nan;
                        metrics(idir,ispeed).susIdx = nan;
                    end
                end
            end


            % copy metrics fields to psth structs
            units(iunit).runPSTH = copyStructFields(metrics,units(iunit).runPSTH);
            clear metrics
        end
    end




    %% fit descriptive functions to onset and offset periods
    %
    rangeThresh = 3;
    reliThresh =-1.645;
    zThresh = 3.29;
    bmp = -195:10:1795;
    parfor iunit= 1:numel(units)
        if ismember(units(iunit).ecephys_structure_acronym, areas)
            for idir = 1:4
                for ispeed = 1:7

                    if units(iunit).statPSTH(idir,ispeed).nTrials>=minTrials

                        % stat
                        response = normalize(units(iunit).statPSTH(idir,ispeed).psth,'range');
                        resp_onset = response(21:51);
                        resp_offset = response(121:151);

                        [units(iunit).statPSTH(idir,ispeed).onset_bestParams,...
                            units(iunit).statPSTH(idir,ispeed).onset_char,...
                            units(iunit).statPSTH(idir,ispeed).onset_bestR2] = fitGaussianTemplates_meanPSTH(resp_onset,5, false);

                        [units(iunit).statPSTH(idir,ispeed).offset_bestParams,...
                            units(iunit).statPSTH(idir,ispeed).offset_char,...
                            units(iunit).statPSTH(idir,ispeed).offset_bestR2] = fitGaussianTemplates_meanPSTH(resp_offset,5, false);

                    else
                        units(iunit).statPSTH(idir,ispeed).onset_bestParams = nan;
                        units(iunit).statPSTH(idir,ispeed).onset_char = nan;
                        units(iunit).statPSTH(idir,ispeed).onset_bestR2 = nan;

                        units(iunit).statPSTH(idir,ispeed).offset_bestParams = nan;
                        units(iunit).statPSTH(idir,ispeed).offset_char = nan;
                        units(iunit).statPSTH(idir,ispeed).offset_bestR2 = nan;
                    end
                    % run
                    if units(iunit).runPSTH(idir,ispeed).nTrials>=minTrials
                        response = normalize(units(iunit).runPSTH(idir,ispeed).psth,'range');
                        resp_onset = response(21:51);
                        resp_offset = response(121:151);

                        [units(iunit).runPSTH(idir,ispeed).onset_bestParams,...
                            units(iunit).runPSTH(idir,ispeed).onset_char,...
                            units(iunit).runPSTH(idir,ispeed).onset_bestR2] = fitGaussianTemplates_meanPSTH(resp_onset,5, false);

                        [units(iunit).runPSTH(idir,ispeed).offset_bestParams,...
                            units(iunit).runPSTH(idir,ispeed).offset_char,...
                            units(iunit).runPSTH(idir,ispeed).offset_bestR2] = fitGaussianTemplates_meanPSTH(resp_offset,5, false);

                    else
                        units(iunit).runPSTH(idir,ispeed).onset_bestParams = nan;
                        units(iunit).runPSTH(idir,ispeed).onset_char = nan;
                        units(iunit).runPSTH(idir,ispeed).onset_bestR2 = nan;

                        units(iunit).runPSTH(idir,ispeed).offset_bestParams = nan;
                        units(iunit).runPSTH(idir,ispeed).offset_char = nan;
                        units(iunit).runPSTH(idir,ispeed).offset_bestR2 = nan;
                    end

                end
            end
        end
    end



    outputFileName =  [dataDirectory,'\psth_' dataFiles(isession).name];

    save(outputFileName, 'units', 'trials', 'runInfo', 'pupilInfo', 'sessionInfo')

end
toc
