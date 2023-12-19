%% analyse train/test decoding

load('rLDA_trainTest_withConsensus.mat')

%% get mean of permutations for each session

for isession = 1:5
    session(isession).meanStat = mean(cat(3,session(isession).stat.perm.perf),3);
    session(isession).meanRun = mean(cat(3,session(isession).run.perm.perf),3);
end

%% plot consensus decoding

allStat = cat(3,session.meanStat);
allRun = cat(3,session.meanRun);
statDiag = diag(mean(allStat,3));
runDiag = diag(mean(allRun,3));


figure
shadedErrorBar(1:191,nanmean(statCons,1),nansem(statCons,1))
hold on
shadedErrorBar(1:191,nanmean(runCons,1),nansem(runCons,1),'lineProps','r')
plot(1:10:191, statDiag,'k')
plot(1:10:191, runDiag, 'r')

title('Consensus decoding (shaded) vs normal decoding (lines)')



%% plot train x test performance for each session

for isession = 1:5
    
%     subplot(5,2,isession*2-1)
subplot(121)
    imagesc(session(isession).meanStat)
    colorbar
    colormap(turbo)
    axis xy
%     subplot(5,2,isession*2)
subplot(122)
    imagesc(session(isession).meanRun)
    colormap(turbo)
    colorbar
    axis xy
end


%% plot train x test performance average of sessions

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
for train_idx=[4 5 6]
figure, hold on
    title(['t = ', num2str(binVector(train_idx))])
shadedErrorBar(1:20, mean(allStat(train_idx,:,:),3), sem(allStat(train_idx,:,:),3),'lineProps', 'k')
shadedErrorBar(1:20, mean(allRun(train_idx,:,:),3), sem(allRun(train_idx,:,:),3),'lineProps', 'r')
plot([train_idx, train_idx], [0.1 0.75], 'm')
ylim([0.1 0.75])
ax = gca; ax.YTick = 0.1:0.1:0.7;
ax.XTick = 1:20; ax.XTick = 0.5:20.5;
ax.XTickLabel = -200:100:1800; % ax.XTickLabel = binVector;
plot(1:20, statDiag,'k:')
 plot([1 20], [1/6 1/6], 'k')
  plot(1:20, runDiag, 'r:')
 xlim([0.5 20.5])
 defaultAxesProperties(gca, true)
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
ax=gca;
    ax.XTick = 1:20; ax.XTick = 1.5:20.5;
ax.XTickLabel = -100:100:1800; % ax.XTickLabel = binVector;
defaultAxesProperties(gca,false)
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



shadedErrorBar(-150:10:1750, nanmean(runCons,1), nansem(runCons,1),'lineProps','r')


%% plot example trial rasters

for itrial = 1:6
figure
data = run.cond(3).catData_sc(:,:,itrial); % time, unit, trial

imagesc(data'), colormap(1-gray), caxis([0 1])
box on
ax = gca;
ax.XTick=[];
ax.YTick=[];
end

%% get normalised performance for each training window using AUC-type measure

bv = -0.15:0.1:1.75;
timeIdx2Use =3:12;
trainIdx=4;

for isession=1:5

    session(isession).normPerf_stat=[];
    session(isession).normPerf_run=[];


diagPerf_stat = diag(session(isession).meanStat);
diagPerf_stat = diagPerf_stat-1/6; 
diagPerf_stat(diagPerf_stat<0)=0;
diagPerf_stat = diagPerf_stat(timeIdx2Use);

diagPerf_run = diag(session(isession).meanRun);
diagPerf_run = diagPerf_run-1/6; 
diagPerf_run(diagPerf_run<0)=0;
diagPerf_run = diagPerf_run(timeIdx2Use);

for itime = 1:timeIdx2Use(end)
    tempPerf = session(isession).meanStat(itime,:)-1/6; tempPerf(tempPerf<0)=0;
    tempPerf = tempPerf(timeIdx2Use);
    session(isession).normPerf_stat(itime) = sum(tempPerf)/sum(diagPerf_stat);

    tempPerf = session(isession).meanRun(itime,:)-1/6; tempPerf(tempPerf<0)=0;
    tempPerf = tempPerf(timeIdx2Use);
    session(isession).normPerf_run(itime) = sum(tempPerf)/sum(diagPerf_run);
end


end

sesh2Use=1:5;

allStat = cat(1,session(sesh2Use).normPerf_stat); 
allRun = cat(1,session(sesh2Use).normPerf_run);  

% relative performance of different training windows
figure, hold on
errorbar(bv(1:timeIdx2Use(end)), mean(allStat,1), sem(allStat,1),'k-')
errorbar(bv(1:timeIdx2Use(end)), mean(allRun,1), sem(allRun,1), 'r-')
ax = gca;
ax.XTick = -0.2:0.1:1;
% ylim([0, 1.1])
xlabel('Training window')
ylabel('Relative Performance to train=test')
defaultAxesProperties(gca, true)

%% plot fractional difference in relataive performance

fracDiffVals  = allRun./allStat

figure, hold on
shadedErrorBar(4:12, mean(fracDiffVals(:,4:12)), sem(fracDiffVals(:,4:12)),'lineProps',{'m'})
plot([4 12], [ 1 1], 'k:')
ax = gca; ax.XTick = 4:12; ax.XTickLabel = 150:100:950;
defaultAxesProperties(gca, true)
%% bar chart of performance for decoder trained 100-200ms
figure, hold on
bar([1 2], [mean(allStat(:,trainIdx)), mean(allRun(:,trainIdx))])
% for isesh = 1:5
%     plot([1 2], [allStat(isesh,trainIdx), allRun(isesh,trainIdx)],'k')
% end
errorbar(1:2, [mean(allStat(:,trainIdx)), mean(allRun(:,trainIdx))], ...
    [sem(allStat(:,trainIdx)), sem(allRun(:,trainIdx))], 'LineStyle', 'none')






%% rm anova for entire stimulus period, with plot

time = 3:12;
statVals = allStat(:,time)'; runVals = allRun(:,time)';
ffVec = cat(1,statVals(:), runVals(:));
stateVec = (cat(1,repelem(1, numel(statVals(:)))', repelem(2,numel(runVals(:)))'));
timeVec = cat(1,repmat(1:numel(time), 1,size(statVals,2))',repmat(1:numel(time), 1,size(runVals,2))');
subjVec = cat(1,repelem(1:size(statVals,2),1, numel(time))', repelem(1:size(runVals,2),1, numel(time))');

[p,tbl,stats,terms] = anovan(ffVec,{timeVec,stateVec,subjVec},'model',2,'random',3,'varnames',{'Time','State','Subj'});


% rel perf over diff training windows

figure, hold on
shadedErrorBar(bv(1:12), mean(allStat), sem(allStat),'lineProps',{'k'})
hold on
shadedErrorBar(bv(1:12), mean(allRun), sem(allRun),'lineProps',{'r'})



%% t-test for early time window (t=100-200)

timeIdx=10;

[h, p] = adtest(allStat(:,timeIdx)-allRun(:,timeIdx));

[h p] = ttest(allStat(:,timeIdx), allRun(:,timeIdx))
