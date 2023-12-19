%% Cross-time decoding analysis


sessionTags =...
    {'M22027', '20220517';...
    'M22029', '20220607';...
    'M22031', '20220524';...
    'M22032', '20220621';...
    'M22033', '20220706'};

for isession = 1:5

    clear stat run

fname = [sessionTags{isession,1},'_', sessionTags{isession,2},'_basic.mat'];

load(fname)

isession



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

% downsample trials and calculate mean rates

minTrial = min(cellfun(@(x) size(x,2), units(1).allSpikes(:)));

for iunit = 1:numel(units)
    units(iunit).allSpikes = cellfun(@(x) x(:,1:minTrial), units(iunit).allSpikes, 'UniformOutput', false);
end

%% get good units

for iunit = 1:numel(units)
    allSpikes = vertcat(units(iunit).allSpikes{:});
    units(iunit).dm_fr = mean(allSpikes(:)).*100;
end

goodUnits = units([units.isi_viol]<=0.1...
    & [units.amplitude_cutoff]<=0.1 & [units.amplitude]>=50 & [units.dm_fr]>=1);


%% get data into decoding format


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



%% train/test in different windwos

do=true;
if do
clear perf

kfold = 2;
options.OptimizeHyperparameters = {'gamma', 'delta'};
options.Gamma=[];
windowSize = 10;
intervals = reshape(1:200,10,20)';
nTest = floor(nTrials/kfold);
nTrain = nTrials-nTest;
options.DiscrimType = 'Linear';

for iperm = 1:5
    iperm

trainidx = randperm(nTrials,nTrain);
testidx = find(~ismember(1:nTrials,trainidx));

% stat
tic
for itrain_int = 1:size(intervals,1)
    itrain_int
    bin2use = intervals(itrain_int,:);
    clear dataStructTrain
    for icond = 1:6
        dataStructTrain(icond).data = stat.cond(icond).catData_sc(bin2use,:,trainidx);
    end
    trainingData = squeeze(mean(cat(3,dataStructTrain.data),1))'; % itrial x iunit
    trainingLabels = repelem(1:6, 1, nTrain)';

     idx2remove = find(var(trainingData,1)<eps);
     trainingData(:,idx2remove)=[];


    Mdl = fitcdiscr(trainingData,trainingLabels,'DiscrimType', ...
        options.DiscrimType,...
        'OptimizeHyperparameters',options.OptimizeHyperparameters,...
        'Gamma',options.Gamma,...
        'HyperparameterOptimizationOptions',...
        struct('ShowPlots',false,...
        'AcquisitionFunctionName','expected-improvement-plus',...
        'UseParallel',true,'Verbose',0));


    % testing
    for itest_int = 1:size(intervals,1)
        bin2use = intervals(itest_int,:);
        clear dataStructTest;
        for icond = 1:6
            dataStructTest(icond).data = stat.cond(icond).catData_sc(bin2use,:,testidx);
        end

        testingData = squeeze(mean(cat(3,dataStructTest.data),1))'; % itrial x iunit
        testingData(:,idx2remove)=[];
        testingLabels = repelem(1:6, 1, nTest)';

        predictions = Mdl.predict(testingData);
        perf(itrain_int, itest_int) = prop(predictions==testingLabels);

    end
end


stat.perm(iperm).perf=perf;
clear perf


% run
tic
for itrain_int = 1:size(intervals,1)
    itrain_int
    bin2use = intervals(itrain_int,:);
    clear dataStructTrain
    for icond = 1:6
        dataStructTrain(icond).data = run.cond(icond).catData_sc(bin2use,:,trainidx);
    end
    trainingData = squeeze(mean(cat(3,dataStructTrain.data),1))'; % itrial x iunit
    trainingLabels = repelem(1:6, 1, nTrain)';

     idx2remove = find(var(trainingData,1)<eps);
     trainingData(:,idx2remove)=[];


    Mdl = fitcdiscr(trainingData,trainingLabels,'DiscrimType', ...
        options.DiscrimType,...
        'OptimizeHyperparameters',options.OptimizeHyperparameters,...
        'Gamma',options.Gamma,...
        'HyperparameterOptimizationOptions',...
        struct('ShowPlots',false,...
        'AcquisitionFunctionName','expected-improvement-plus',...
        'UseParallel',true,'Verbose',0));


    % testing
    for itest_int = 1:size(intervals,1)
        bin2use = intervals(itest_int,:);
        clear dataStructTest;
        for icond = 1:6
            dataStructTest(icond).data = run.cond(icond).catData_sc(bin2use,:,testidx);
        end

        testingData = squeeze(mean(cat(3,dataStructTest.data),1))'; % itrial x iunit
        testingData(:,idx2remove)=[];
        testingLabels = repelem(1:6, 1, nTest)';

        predictions = Mdl.predict(testingData);
        perf(itrain_int, itest_int) = prop(predictions==testingLabels);

    end
end

run.perm(iperm).perf = perf;


end


session(isession).stat = stat;
session(isession).run = run;


end

end



%% plots

for isession = 1:5
    session(isession).meanStat = mean(cat(3,session(isession).stat.perm.perf),3);
    session(isession).meanRun = mean(cat(3,session(isession).run.perm.perf),3);
end



%% plot each session

for isession = 1:5
    figure
    subplot(121)
    imagesc(session(isession).meanStat)
    colorbar
    colormap(turbo)
    axis xy
    subplot(122)
    imagesc(session(isession).meanRun)
    colormap(turbo)
    colorbar
    axis xy
end


%% plot average of sessions

allStat = cat(3,session.meanStat);
allRun = cat(3,session.meanRun);

figure
subplot(121), hold on
imagesc(mean(allStat,3)');
xlabel('Testing window')
ylabel('Training window')
caxis([1/6 0.5])
axis xy
colorbar
plot([3 3], [0.5 20.5], 'r')
plot([13 13], [0.5 20.5], 'r')
plot([0.5 20.5], [3 3], 'r')
plot([0.5 20.5], [13 13],'r')
xlim([0.5 20.5]), ylim([0.5 20.5])
ax = gca; ax.XTick = 1:20; ax.YTick = 1:20;
defaultAxesProperties(gca, false)

subplot(122), hold on
imagesc(mean(allRun,3)');
xlabel('Testing window')
ylabel('Training window')
caxis([1/6 0.7])
axis xy
colorbar
plot([3 3], [0.5 20.5], 'r')
plot([13 13], [0.5 20.5], 'r')
plot([0.5 20.5], [3 3], 'r')
plot([0.5 20.5], [13 13],'r')
xlim([0.5 20.5]), ylim([0.5 20.5])
ax = gca; ax.XTick = 1:20; ax.YTick = 1:20;
defaultAxesProperties(gca, false)

%% plot mean performance for each train/test window 
allStat = cat(3,session.meanStat);
allRun = cat(3,session.meanRun);
statDiag = diag(mean(allStat,3));
runDiag = diag(mean(allRun,3));
binVector = -150:100:1750;

%figure
for train_idx=1:20
figure, hold on
    title(['t = ', num2str(binVector(train_idx))])
shadedErrorBar(1:20, mean(allStat(train_idx,:,:),3), sem(allStat(train_idx,:,:),3),'lineProps', 'k')
shadedErrorBar(1:20, mean(allRun(train_idx,:,:),3), sem(allRun(train_idx,:,:),3),'lineProps', 'r')
plot([train_idx, train_idx], [0.1 0.75], 'm')
ylim([0.1 0.75])
ax = gca; ax.YTick = 0.1:0.1:0.7;
ax.XTick = 1:20; ax.XTick = 1.5:20.5;
ax.XTickLabel = -100:100:1800; % ax.XTickLabel = binVector;
 plot(1:20, statDiag,'k:')
 plot(1:20, runDiag, 'r:')
 xlim([2.5 20.5])
 defaultAxesProperties(gca, false)
end



%% plot fractional performance of different train test windows

allStat = cat(3,session.meanStat);
allRun = cat(3,session.meanRun);
statDiag = diag(mean(allStat,3));
runDiag = diag(mean(allRun,3));

%$figure
for train_idx=1:20
    figure, hold on
        title(['t = ', num2str(binVector(train_idx))])

    statDiff = (mean(allStat(train_idx,:,:),3)-1/6)' ./ (statDiag-1/6);
    runDiff = (mean(allRun(train_idx,:,:),3)-1/6)' ./ (runDiag-1/6);
    plot(statDiff,'k')
    plot(runDiff,'r')
    ylim([0 1.3])
    plot([train_idx, train_idx], [0 1.3], 'm')
     xlim([0.5 20.5])

end



%% plot fractional performance of different train test windows by sesssion

% add checks for: diag performance > chance 
% train/test performance > chance

% if diag performance is < chance, then all values must be nan for that
% if train/test is < chance but diag is > chance, then values should be 0.

threshPerf = 0.2083;

for isession = 1:5
    figure

    statDiag = diag(session(isession).meanStat);
    runDiag = diag(session(isession).meanRun);
    statDiag(statDiag<=threshPerf) = nan;
    runDiag(runDiag<=threshPerf) = nan;


    for train_idx=1:20
    subplot(4,5,train_idx), hold on
    
    thisStatPerf = session(isession).meanStat(train_idx,:); thisStatPerf(thisStatPerf<=threshPerf)=1/6;
    thisRunPerf = session(isession).meanRun(train_idx,:); thisRunPerf(thisRunPerf<=threshPerf)=1/6;

    
    session(isession).statDiff(train_idx,:) = (thisStatPerf'-1/6)./(statDiag-1/6);
    session(isession).runDiff(train_idx,:) = (thisRunPerf'-1/6)./(runDiag-1/6);
    plot(session(isession).statDiff(train_idx,:),'k')
    plot(session(isession).runDiff(train_idx,:),'r')
    plot([train_idx, train_idx], [0 1.3], 'm')

    ylim([0 1.3])
    end
  

end

allStatDiff = cat(3,session.statDiff);
allRunDiff = cat(3, session.runDiff);
% 
allStatDiff(allStatDiff>1) = 1;
allRunDiff(allRunDiff>1) = 1;


figure
for train_idx=1:20
    figure, hold on
            title(['t = ', num2str(binVector(train_idx))])

    shadedErrorBar(1:20, nanmean(allStatDiff(train_idx,:,:),3),nansem(allStatDiff(train_idx,:,:),3), 'lineProps','k.-');
   shadedErrorBar(1:20, nanmean(allRunDiff(train_idx,:,:),3),nansem(allRunDiff(train_idx,:,:),3), 'lineProps','r.-');
           plot([train_idx, train_idx], [0 1.3], 'm')
 
   ylim([0 1])
    xlim([2.5 20.5])
end


    %%

    isession = 2;
allStat = session(isession).meanStat;
allRun = session(isession).meanRun;

statDiag = diag(allStat);
runDiag = diag(allRun);

figure
for train_idx=1:20
    subplot(4,5,train_idx), hold on
    title(['t = ', num2str(binVector(train_idx)+50)])
shadedErrorBar(1:20, mean(allStat(train_idx,:,:),3), sem(allStat(train_idx,:,:),3),'lineProps', 'k')
shadedErrorBar(1:20, mean(allRun(train_idx,:,:),3), sem(allRun(train_idx,:,:),3),'lineProps', 'r')
plot([train_idx, train_idx], [0.1 0.7], 'm')
ylim([0.1 0.7])
ax = gca; ax.YTick = 0.1:0.1:0.7;
ax.XTick = 1:20; ax.XTickLabel = binVector;
 plot(1:20, statDiag,'k:')
 plot(1:20, runDiag, 'r:')
 xlim([0.5 20.5])
end

figure
for train_idx=1:20
    subplot(4,5,train_idx), hold on
    statDiff = (mean(allStat(train_idx,:,:),3)-1/6)' ./ (statDiag-1/6);
    runDiff = (mean(allRun(train_idx,:,:),3)-1/6)' ./ (runDiag-1/6);
    plot(statDiff,'k')
    plot(runDiff,'r')
    ylim([0 1.3])
end
