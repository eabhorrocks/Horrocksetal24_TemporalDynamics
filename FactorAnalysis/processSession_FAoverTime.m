function processSession_FAoverTime(inputFileName,outputFileName,dataDir)

load(fullfile(dataDir,inputFileName))
%% get binned spike counts
tic
trialDur = 1;
tolerance = 0.05; % within 50ms
tsd = trials.Speed2D;
trialdurs = [tsd.PDend]-[tsd.PDstart];
invalidDurs_idx = abs(trialdurs-1)>tolerance;
tsd(invalidDurs_idx)=[];

tsd = tsd([tsd.Contrast1]==1); % only using full contrast trials

% split trials by state according to wheel data
run_idx = find(cellfun(@(x) prop(x>0.5)>=0.75 & mean(x)>3, {tsd.WheelSpeed}));
stat_idx = find(cellfun(@(x) prop(x<3)>=0.75 & mean(x)<0.5, {tsd.WheelSpeed}));
[tsd.runFlag] = deal([nan]);

[tsd(stat_idx).runFlag] = deal([0]);
[tsd(run_idx).runFlag] = deal([1]);

runTrials = tsd(cellfun(@(x) prop(x>0.5)>=0.75 & mean(x)>3, {tsd.WheelSpeed}));
statTrials = tsd(cellfun(@(x) prop(x<3)>=0.75 & mean(x)<0.5, {tsd.WheelSpeed}));

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

options.intervalStart = -0.2;
options.intervalEnd = 1.8;
options.binSpacing=0.01; %

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

minTrial = min(cellfun(@(x) size(x,2), units(1).allSpikes(:)));

for iunit = 1:numel(units)
    units(iunit).allSpikes = cellfun(@(x) x(:,1:minTrial), units(iunit).allSpikes, 'UniformOutput', false);
    units(iunit).allSpikesShuffled = cellfun(@(x) x(:,randperm(size(x,2))), units(iunit).allSpikes, 'UniformOutput', false);
end

minBlankTrial = min(cellfun(@(x) size(x,2), units(1).blankSpikes(:)));

for iunit = 1:numel(units)
    units(iunit).blankSpikes = cellfun(@(x) x(:,1:minBlankTrial), units(iunit).blankSpikes, 'UniformOutput', false);
end
toc

% filter units for mean firing rate
for iunit =1:numel(units)
    units(iunit).dmFR = mean(cellfun(@(x) mean(x,'all'), units(iunit).allSpikes),1).*(1/options.binSpacing);
    units(iunit).blankFR = mean(cellfun(@(x) mean(x,'all'), units(iunit).blankSpikes),1).*(1/options.binSpacing);
end





%% get trials struct for stat and run

minFR = 1;

tempUnits = units([units.isi_viol]<=0.1...
    & [units.amplitude_cutoff]<=0.1 & [units.amplitude]>=50 & all(cat(1,units.dmFR)>minFR,2)');

%% factor analysis over time. 

windowSize=20;

handles=[];
minFR = 1;

for iint = 1:200-(windowSize-1)
    bin2use = iint:(iint+(windowSize-1));
    %bin2use = intervals(iint,:);

    %%% get spike counts %%%
    D_stat=struct;
    ii = 0;
    for ispeed = 1:6
        for itrial = 1:minTrial
            ii = ii+1;
            for iunit = 1:numel(tempUnits)
                D_stat(ii).data(iunit,:) = sum(tempUnits(iunit).allSpikes{ispeed,1}(bin2use,itrial),1).*5;
                D_stat(ii).condition = num2str(ispeed);
            end
        end
    end

    D_run=struct;
    ii = 0;
    for ispeed = 1:6
        for itrial = 1:minTrial
            ii = ii+1;
            for iunit = 1:numel(tempUnits)
                D_run(ii).data(iunit,:) = sum(tempUnits(iunit).allSpikes{ispeed,2}(bin2use,itrial),1).*5;
                D_run(ii).condition = num2str(ispeed);
            end
        end
    end

     %%% shuffled %%%
    D_statShuf=struct;
    ii = 0;
    for ispeed = 1:6
        for itrial = 1:minTrial
            ii = ii+1;
            for iunit = 1:numel(tempUnits)
                D_statShuf(ii).data(iunit,:) = sum(tempUnits(iunit).allSpikesShuffled{ispeed,1}(bin2use,itrial),1).*5;
                D_statShuf(ii).condition = num2str(ispeed);
            end
        end
    end

    D_runShuf=struct;
    ii = 0;
    for ispeed = 1:6
        for itrial = 1:minTrial
            ii = ii+1;
            for iunit = 1:numel(tempUnits)
                D_runShuf(ii).data(iunit,:) = sum(tempUnits(iunit).allSpikesShuffled{ispeed,2}(bin2use,itrial),1).*5;
                D_runShuf(ii).condition = num2str(ispeed);
            end
        end
    end



    %%% do cv to get qOpt %%%
    D = D_stat;
    idx2remove= find(mean(cat(2,D.data),2)<minFR);
    for itrial = 1:numel(D)
        D(itrial).data(idx2remove) = [];
    end
    % estimate dimensionality using FA
    candidateDims = 1:10; %
    alg = 3; % 1PCA, 2PPCA, 3FA, 4 LDA, 5 GPFA
    [projs, mse, like] = cvreducedims_edd(D, alg, candidateDims, handles);
    [~, idx] = max(like); % find q that maximises likelihood of data
    q = candidateDims(idx);
    [newD, C, lat, explained, params] = reducedims_EH(D,alg, q, handles); % do FA using q dims

    statDims(iint) = compute_dshared(params);
    statLoadingSim{iint} = compute_load_sim(params);
    statSV(iint) = compute_perc_shared(params);


    %%% do cv to get qOpt %%%
    D = D_run;
    idx2remove= find(mean(cat(2,D.data),2)<minFR);
    for itrial = 1:numel(D)
        D(itrial).data(idx2remove) = [];
    end
    % estimate dimensionality using FA
    candidateDims = 1:15; %
    alg = 3; % 1PCA, 2PPCA, 3FA, 4 LDA, 5 GPFA
    [projs, mse, like] = cvreducedims_edd(D, alg, candidateDims, handles);
    [~, idx] = max(like); % find q that maximises likelihood of data
    q = candidateDims(idx);
    [newD, C, lat, explained, params] = reducedims_EH(D,alg, q, handles); % do FA using q dims


    runDims(iint) = compute_dshared(params);
    runLoadingSim{iint} = compute_load_sim(params);
    runSV(iint) = compute_perc_shared(params);


    %%% shuffled %%%
    D = D_statShuf;
    idx2remove= find(mean(cat(2,D.data),2)<minFR);
    for itrial = 1:numel(D)
        D(itrial).data(idx2remove) = [];
    end
    % estimate dimensionality using FA
    candidateDims = 1:10; %
    alg = 3; % 1PCA, 2PPCA, 3FA, 4 LDA, 5 GPFA
    [projs, mse, like] = cvreducedims_edd(D, alg, candidateDims, handles);
    [~, idx] = max(like); % find q that maximises likelihood of data
    q = candidateDims(idx);
    [newD, C, lat, explained, params] = reducedims_EH(D,alg, q, handles); % do FA using q dims

    statShufDims(iint) = compute_dshared(params);
    statShufLoadingSim{iint} = compute_load_sim(params);
    statShufSV(iint) = compute_perc_shared(params);


    %%% do cv to get qOpt %%%
    D = D_runShuf;
    idx2remove= find(mean(cat(2,D.data),2)<minFR);
    for itrial = 1:numel(D)
        D(itrial).data(idx2remove) = [];
    end
    % estimate dimensionality using FA
    candidateDims = 1:15; %
    alg = 3; % 1PCA, 2PPCA, 3FA, 4 LDA, 5 GPFA
    [projs, mse, like] = cvreducedims_edd(D, alg, candidateDims, handles);
    [~, idx] = max(like); % find q that maximises likelihood of data
    q = candidateDims(idx);
    [newD, C, lat, explained, params] = reducedims_EH(D,alg, q, handles); % do FA using q dims


    runShufDims(iint) = compute_dshared(params);
    runShufLoadingSim{iint} = compute_load_sim(params);
    runShufSV(iint) = compute_perc_shared(params);

end


%% save data


session.statDims = statDims;
session.statLoadingSim = statLoadingSim;
session.statSV = statSV;

session.statShufDims = statShufDims;
session.statShufLoadingSim = statShufLoadingSim;
session.statShufSV = statShufSV;

session.runDims = runDims;
session.runLoadingSim = runLoadingSim;
session.runSV = runSV;

session.runShufDims = runShufDims;
session.runShufLoadingSim = runShufLoadingSim;
session.runShufSV = runShufSV;

save(fullfile(dataDir,outputFileName),'session')

end
