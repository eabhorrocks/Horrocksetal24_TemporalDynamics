function processSession_popAnalysis_fitSep(inputFileName,outputFileName,dataDir)

load(fullfile(dataDir,inputFileName))

%% process trials
% filter trials to the ones we want to analyse and sort according to
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

%% pre-proces binned spike counts for FA

% split stat/run trials
for iunit = 1:numel(tempUnits)
    tempUnits(iunit).dmSC_10ms_stat = tempUnits(iunit).allSpikes(:,1)';
    tempUnits(iunit).dmSC_10ms_run = tempUnits(iunit).allSpikes(:,2)';
end

% use only good units
units = tempUnits([units.isi_viol]<=0.1...
    & [tempUnits.amplitude_cutoff]<=0.1 & [tempUnits.amplitude]>=50);
%% preprocess data

minFR = 1;
smoothWidth = 175; % in ms
bin_N = 1; % bin every N elements
binSize = 10; % bin size in ms
binVector = -200:binSize*bin_N:(1800-binSize*bin_N);
nConds = 12;
nDecodingBins = 5; %  n*10 ms decoding window (e.g. 5 =50ms)

plotFlag = false;

options.kfold = 9999;

vispu_new = struct;
scName = 'dmSC_10ms_stat';
for iunit =1:numel(units)

    % optionally bin data
    %binnedData = cellfun(@(x) sumEveryN(sqrt(x), bin_N, 1)', units(iunit).(scName), 'UniformOutput', false);
    binnedData = units(iunit).(scName)';
    % sqrt and smooth data for FA
    processedData = cellfun(@(x) smoothdata(sqrt(x), 1, 'gaussian', smoothWidth/(bin_N*binSize)), binnedData,'UniformOutput', false);
    % get average for calculting mean firing rate
    average_smooth = cellfun(@(x) smoothdata(mean(x,2),'gaussian', smoothWidth/(bin_N*binSize)), binnedData, 'UniformOutput', false);

    mu = mean(vertcat(average_smooth{:}));
    mu = mu*(1000/(binSize*bin_N));

    % zscore data
    tempData = horzcat(processedData{:});
    zmean = mean(tempData(:));
    zstd = std(tempData(:));
    processedData = cellfun(@(x) arrayfun(@(x) (x-zmean)/(zstd), x), processedData,'UniformOutput',false);
    vispu_new(iunit).processedData = processedData(:);

    % for gpfa
    vispu_new(iunit).binnedData = binnedData(:);

    vispu_new(iunit).mu = mu;
    vispu_new(iunit).ID = string(units(iunit).cluster_id);
end

% remove low fr units
vispu_new([vispu_new.mu]<minFR) = [];

% get number fo available trials
nTrialsByCond = cellfun(@(x) size(x,2), vispu_new(1).processedData);


% generate input data


cols = inferno(7);


% get D struct for stat and run separately
% input data for Factor Analysis (pre-smoothed))
D_stat = struct;
ii = 0;
for ispeed = 1:6
    for itrial = 1:size(vispu_new(1).processedData{ispeed},2)
        ii = ii+1;
        for iunit = 1:numel(vispu_new)
            D_stat(ii).data(iunit,:) = vispu_new(iunit).processedData{ispeed}(:,itrial);
            D_stat(ii).condition = num2str(ispeed);
            %                D(ii).epochColors = cols(ispeed,:);
        end
    end
end




%% DIM REDUCTION (DataHigh toolbox)

handles = []; % no longer needed


% stat 
candidateDims = 10:37; % 1:50, 
alg = 3; % 1PCA, 2PPCA, 3FA, 4 LDA, 5 GPFA
[projs, mse, like] = cvreducedims_edd(D_stat, alg, candidateDims, handles);
[~, idx] = max(like); % find q that maximises likelihood of data
q = candidateDims(idx);
[newD, C, lat, explained, params] = reducedims_EH(D_stat,alg, q, handles); % do FA using q dims

s.q = q;
s.qOpt = compute_dshared(params);
s.LoadingSim = compute_load_sim(params);
s.SV = compute_perc_shared(params);
s.params = params;
s.propSharedVariance = [lat(1); diff(lat)];
s.loadings = C;
s.nUnits = size(C,1);
s.D=D_stat;
s.explained = explained;

for itrial = 1:numel(newD)
    newD(itrial).y = D_stat(itrial).data;
end

loadings = C;

nUnits = numel(tempUnits);
s.nUnits = nUnits;


nBins = 200/bin_N;
cond = struct; nConds = 6;
for icond = 1:nConds
    idx = find(strcmp({newD.condition},num2str(icond)));
    cond(icond).catData = cat(3,newD(idx).data);
    cond(icond).catData = permute(cond(icond).catData, [2, 1, 3]);
    cond(icond).meanTrajectory = mean(cond(icond).catData,3);
    Dout(icond).data = cond(icond).meanTrajectory;
    %cond(ispeed).catData_sc = cat(3,newD(idx).y); % input spike counts
    %cond(ispeed).catData_sc = permute(cond(ispeed).catData_sc, [2, 1, 3]);
end

s.cond= cond;

stat.s = s;
clear s


%% run 


vispu_new = struct;
scName = 'dmSC_10ms_run';
for iunit =1:numel(units)

    % optionally bin data
    %binnedData = cellfun(@(x) sumEveryN(sqrt(x), bin_N, 1)', units(iunit).(scName), 'UniformOutput', false);
    binnedData = units(iunit).(scName)';
    % sqrt and smooth data for FA
    processedData = cellfun(@(x) smoothdata(sqrt(x), 1, 'gaussian', smoothWidth/(bin_N*binSize)), binnedData,'UniformOutput', false);
    % get average for calculting mean firing rate
    average_smooth = cellfun(@(x) smoothdata(mean(x,2),'gaussian', smoothWidth/(bin_N*binSize)), binnedData, 'UniformOutput', false);

    mu = mean(vertcat(average_smooth{:}));
    mu = mu*(1000/(binSize*bin_N));

    % zscore data
    tempData = horzcat(processedData{:});
    zmean = mean(tempData(:));
    zstd = std(tempData(:));
    processedData = cellfun(@(x) arrayfun(@(x) (x-zmean)/(zstd), x), processedData,'UniformOutput',false);
    vispu_new(iunit).processedData = processedData(:);

    % for gpfa
    vispu_new(iunit).binnedData = binnedData(:);

    vispu_new(iunit).mu = mu;
    vispu_new(iunit).ID = string(units(iunit).cluster_id);
end

% remove low fr units
vispu_new([vispu_new.mu]<minFR) = [];

% get number fo available trials
nTrialsByCond = cellfun(@(x) size(x,2), vispu_new(1).processedData);


D_run = struct;
ii = 0;
for ispeed = 1:6
    for itrial = 1:size(vispu_new(1).processedData{ispeed},2)
        ii = ii+1;
        for iunit = 1:numel(vispu_new)
            D_run(ii).data(iunit,:) = vispu_new(iunit).processedData{ispeed}(:,itrial);
            D_run(ii).condition = num2str(ispeed);
            %                D(ii).epochColors = cols(ispeed,:);
        end
    end
end




handles = []; % no longer needed
candidateDims = 10:50; % 1:50, 
alg = 3; % 1PCA, 2PPCA, 3FA, 4 LDA, 5 GPFA
[projs, mse, like] = cvreducedims_edd(D_run, alg, candidateDims, handles);
[~, idx] = max(like); % find q that maximises likelihood of data
q = candidateDims(idx);
[newD, C, lat, explained, params] = reducedims_EH(D_run,alg, q, handles); % do FA using q dims

s.q = q;
s.qOpt = compute_dshared(params);
s.LoadingSim = compute_load_sim(params);
s.SV = compute_perc_shared(params);
s.params = params;
s.propSharedVariance = [lat(1); diff(lat)];
s.loadings = C;
s.nUnits = size(C,1);
s.D=D_run;
s.explained = explained;

for itrial = 1:numel(newD)
    newD(itrial).y = D_run(itrial).data;
end

loadings = C;

nUnits = numel(tempUnits);
s.nUnits = nUnits;


nBins = 200/bin_N;
cond = struct; nConds = 6;
for icond = 1:nConds
    idx = find(strcmp({newD.condition},num2str(icond)));
    cond(icond).catData = cat(3,newD(idx).data);
    cond(icond).catData = permute(cond(icond).catData, [2, 1, 3]);
    cond(icond).meanTrajectory = mean(cond(icond).catData,3);
    Dout(icond).data = cond(icond).meanTrajectory;
    %cond(ispeed).catData_sc = cat(3,newD(idx).y); % input spike counts
    %cond(ispeed).catData_sc = permute(cond(ispeed).catData_sc, [2, 1, 3]);
end

s.cond= cond;

run.s = s;


%% save session data


session.stat = stat;
session.run = run;

try
save(fullfile(dataDir,outputFileName),'session')
catch
    dummyVar=1;
save(fullfile(dataDir,outputFileName),'dummyVar','session', '-v7.3')
end


end

