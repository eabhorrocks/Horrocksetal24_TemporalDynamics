%%


%% i/o args
% output args:
% PSTH
% z-scored PSTH
% ?reliability
% plotHandle


% input args:
% stimulus trial event times (on + off times)
% blank trial event times
% spike times


spikeTimes = unitex.spiketimes.*1000;
statTrials = tsd(cellfun(@(x) prop(x<3)>=0.8 & mean(x)<0.5, {tsd.WheelSpeed}));
statTrials_blank = statTrials([statTrials.numDots1]==0);

speeds = [0, -16, -32, -64, -128, -256];

for ispeed = 1:numel(speeds)
    stimTrials = statTrials([statTrials.Contrast1]==1 & [statTrials.VelX1]==speeds(ispeed));
    intervalsCell{ispeed}= [vertcat(stimTrials.PDstart), vertcat(stimTrials.PDend)].*1000;
end

blankIntervals = [vertcat(statTrials_blank.PDstart), vertcat(statTrials_blank.PDend)].*1000;

%%
options.binWidth = 10;
options.preTime = 200;
options.postTime = 800;
options.smoothType = 'gaussian';
options.smoothWidth = 175;


options.getReliability = true;
options.kfold = 1000; % set to number>nTrials to do Leave-one-out. % smaller the k-fold, higher the reliability score (more training trials)
options.distType = @(x, Y)dtw(x,Y,20,'euclidean'); % any pdist argument, e.g. 'correlation', 'euclidean'
%options.distType = 'correlation';

options.nReliPerms = 10;
options.nShuffle = 10;


options.plot = true;
options.cols = inferno(7); options.cols(1,:) = [];
options.plotFRrange = [0 60];
options.plot95CI = false;

figure
[psth, blank, bmp] = makePSTH(spikeTimes,intervalsCell,blankIntervals,options);

%% get metrics for PSTH

% input args for check responsive
blankResp = mean(smoothdata(unitex.blankTrials_stat{1}.*100,1,'gaussian', 17.5),2);
distType = 'euclidean';

% input args for PSTH metrics
zThresh = 3.29; % 99.9% CI
minRespEnd = 1000;
baselineFR = mean(blankResp(:));

for ipsth = 1:numel(psth)
    
    % smoothed stim trials
    stimTrials = smoothdata(unitex.trials_stat{ispeed}.*100,1,'gaussian', 17.5);
    [pval, stats, distRatio] = checkResponsive(stimTrials, blankResp, distType);

    
    psth(ipsth).resp_pval = checkResponsive(stimTrials, blankResp, distType);
    metrics(ipsth) = calc_PSTHmetrics(psth(ipsth).psth, psth(ipsth).zpsth, bmp, baselineFR, zThresh, minRespEnd);
end

psth = copyStructFields(metrics,psth);

%%




