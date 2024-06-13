%% analyse cross-time decoding
sessionTags = {'M22027', '20220517';...
    'M22029', '20220607';...
    'M22031', '20220524';
    'M22032', '20220621';
    'M22033', '20220706'};

dataDir = 'C:\Users\edward.horrocks\Documents\Code\V1Dynamics\Data\basic_111022';


for isession = 1:size(sessionTags,1)
    fname = [sessionTags{isession,1},'_', sessionTags{isession,2},'_crossTime_SigNoise.mat']; %PSTH_noEye %psth3rds

    load(fullfile(dataDir,fname))
    [units.session] = deal(isession);

    s(isession) = session;



end

%%
for ipop = 1:4

allStat=[];
allRun=[];
allStatShuf=[];
allRunShuf=[];

for isession = 1:5
    s(isession).allStat = [];
    s(isession).allRun = [];
    s(isession).allStatShuf = [];
    s(isession).allRunShuf = [];
    try
        for irep =1:numel(s(isession).popSize(ipop).rep)
            thisStatPerf = mean(cat(3,s(isession).popSize(ipop).rep(irep).stat.perm.perf),3);
            thisRunPerf = mean(cat(3,s(isession).popSize(ipop).rep(irep).run.perm.perf),3);

            s(isession).allStat = cat(3,s(isession).allStat,thisStatPerf);
            s(isession).allRun = cat(3,s(isession).allRun,thisRunPerf);

            thisStatShufPerf = mean(cat(3,s(isession).popSize(ipop).rep(irep).statShuf.perm.perf),3);
            thisRunShufPerf = mean(cat(3,s(isession).popSize(ipop).rep(irep).runShuf.perm.perf),3);

            s(isession).allStatShuf = cat(3,s(isession).allStatShuf,thisStatShufPerf);
            s(isession).allRunShuf = cat(3,s(isession).allRunShuf,thisRunShufPerf);

        end
    catch
    end

%     figure
%     subplot(121)
%     imagesc(mean(s(isession).allStat,3))
%     colormap(turbo)
%     colorbar
%     axis xy
% 
%     subplot(122)
%     imagesc(mean(s(isession).allRun,3))
%     colormap(turbo)
%     colorbar
%     axis xy


end






%% plot average of sessions
% 
% allStat = cat(3,s.allStat);
% allRun = cat(3,s.allRun);
% 
% figure
% subplot(121), hold on
% imagesc(mean(allStat,3)');
% xlabel('Testing window')
% ylabel('Training window')
% caxis([1/6 0.4])
% axis xy
% colorbar
% plot([3 3], [0.5 20.5], 'r')
% plot([13 13], [0.5 20.5], 'r')
% plot([0.5 20.5], [3 3], 'r')
% plot([0.5 20.5], [13 13],'r')
% xlim([0.5 20.5]), ylim([0.5 20.5])
% ax = gca; ax.XTick = 1:20; ax.YTick = 1:20;
% defaultAxesProperties(gca, false)
% 
% subplot(122), hold on
% imagesc(mean(allRun,3)');
% xlabel('Testing window')
% ylabel('Training window')
% caxis([1/6 0.6])
% axis xy
% colorbar
% plot([3 3], [0.5 20.5], 'r')
% plot([13 13], [0.5 20.5], 'r')
% plot([0.5 20.5], [3 3], 'r')
% plot([0.5 20.5], [13 13],'r')
% xlim([0.5 20.5]), ylim([0.5 20.5])
% ax = gca; ax.XTick = 1:20; ax.YTick = 1:20;
% defaultAxesProperties(gca, false)
% title(ipop)





%% plot mean performance for each train/test window
allStat = cat(3,s.allStat);
allRun = cat(3,s.allRun);
statDiag = diag(mean(allStat,3));
runDiag = diag(mean(allRun,3));
binVector = -150:100:1750;

%figure
for train_idx=[4]% 5 6]
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
title(ipop)

%% plot fractional performance of different train test windows by sesssion

% add checks for: diag performance > chance
% train/test performance > chance

% if diag performance is < chance, then all values must be nan for that
% if train/test is < chance but diag is > chance, then values should be 0.

timeIdx2Use =3:12;
trainIdx=4;

for isession = 1:5
    s(isession).statDiff=[];
    s(isession).runDiff=[];

    s(isession).statShufDiff=[];
    s(isession).runShufDiff=[];

    if numel(s(isession).popSize)>=ipop

        for irep =1:size(s(isession).allStat,3)

            thisStatPerf = s(isession).allStat(:,:,irep);
            thisRunPerf = s(isession).allRun(:,:,irep);

            thisStatShufPerf = s(isession).allStatShuf(:,:,irep);
            thisRunShufPerf = s(isession).allRunShuf(:,:,irep);


            diagPerf_stat = diag(thisStatPerf);
            diagPerf_stat = diagPerf_stat-1/6;
            diagPerf_stat(diagPerf_stat<0)=0;
            diagPerf_stat = diagPerf_stat(timeIdx2Use);

            diagPerf_run = diag(thisRunPerf);
            diagPerf_run = diagPerf_run-1/6;
            diagPerf_run(diagPerf_run<0)=0;
            diagPerf_run = diagPerf_run(timeIdx2Use);

            diagPerf_statShuf = diag(thisStatShufPerf);
            diagPerf_statShuf = diagPerf_statShuf-1/6;
            diagPerf_statShuf(diagPerf_statShuf<0)=0;
            diagPerf_statShuf = diagPerf_statShuf(timeIdx2Use);

            diagPerf_runShuf = diag(thisRunShufPerf);
            diagPerf_runShuf = diagPerf_runShuf-1/6;
            diagPerf_runShuf(diagPerf_runShuf<0)=0;
            diagPerf_runShuf = diagPerf_runShuf(timeIdx2Use);

            %         thisMeanStatPerf = thisMeanStatPerf(timeIdx2Use);
            %         thisMeanRunPerf = thisMeanRunPerf(timeIdx2Use);

            for itime = 1:timeIdx2Use(end)
                tempPerf = thisStatPerf(itime,:)-1/6; tempPerf(tempPerf<0)=0;
                tempPerf = tempPerf(timeIdx2Use);
                thisDiffStatPerf(itime) = sum(tempPerf)/sum(diagPerf_stat);

                tempPerf = thisRunPerf(itime,:)-1/6; tempPerf(tempPerf<0)=0;
                tempPerf = tempPerf(timeIdx2Use);
                thisDiffRunPerf(itime) = sum(tempPerf)/sum(diagPerf_run);

                tempPerf = thisStatShufPerf(itime,:)-1/6; tempPerf(tempPerf<0)=0;
                tempPerf = tempPerf(timeIdx2Use);
                thisDiffStatShufPerf(itime) = sum(tempPerf)/sum(diagPerf_stat);

                tempPerf = thisRunShufPerf(itime,:)-1/6; tempPerf(tempPerf<0)=0;
                tempPerf = tempPerf(timeIdx2Use);
                thisDiffShufRunPerf(itime) = sum(tempPerf)/sum(diagPerf_run);
            end

            s(isession).statDiff=cat(1,s(isession).statDiff,thisDiffStatPerf);
            s(isession).runDiff=cat(1,s(isession).runDiff,thisDiffRunPerf);

            s(isession).statShufDiff=cat(1,s(isession).statShufDiff,thisDiffStatShufPerf);
            s(isession).runShufDiff=cat(1,s(isession).runShufDiff,thisDiffShufRunPerf);
        end

    end


end

%% plot normalised performance for all reps

popSizeVector = [s(1).popSize.nUnits];

allStatDiff = cat(1,s.statDiff);
allRunDiff = cat(1,s.runDiff);

allStatDiff(allStatDiff>1) = 1;
allRunDiff(allRunDiff>1) = 1;

allStatShufDiff = cat(1,s.statShufDiff);
allRunShufDiff = cat(1,s.runShufDiff);

allStatShufDiff(allStatShufDiff>1) = 1;
allRunShufDiff(allRunShufDiff>1) = 1;

trainIdx = 4;
% figure, hold on
% histogram(allStatDiff(:,trainIdx),'FaceColor','k','BinEdges',[0:0.05:1])
% histogram(allRunDiff(:,trainIdx),'Facecolor','r','BinEdges',[0:0.05:1])
% plot([0 1],[0 1],'k:')
% 
% plot(allStatDiff(:,trainIdx),allRunDiff(:,trainIdx),'k.')
% title(num2str(popSizeVector(ipop)))


%% plot median rel perf for all time windows


for itime = 1:12
[statDiff_CI_lower(itime), statDiff_CI_upper(itime)] = medianCI(allStatDiff(:,itime), 1);
[runDiff_CI_lower(itime), runDiff_CI_upper(itime)] = medianCI(allRunDiff(:,itime), 1);

[statShufDiff_CI_lower(itime), statShufDiff_CI_upper(itime)] = medianCI(allStatShufDiff(:,itime), 1);
[runShufDiff_CI_lower(itime), runShufDiff_CI_upper(itime)] = medianCI(allRunShufDiff(:,itime), 1);
end

figure(122), hold on
title(num2str(popSizeVector(ipop)))
errorbar(1:12, median(allStatDiff,1), median(allStatDiff,1)-statDiff_CI_lower, statDiff_CI_upper-median(allStatDiff,1),'k')
errorbar(1:12, median(allRunDiff,1), median(allRunDiff,1)-runDiff_CI_lower, runDiff_CI_upper-median(allRunDiff,1),'r')

%errorbar(1:12, median(allStatShufDiff,1), median(allStatShufDiff,1)-statShufDiff_CI_lower, statShufDiff_CI_upper-median(allStatShufDiff,1),'k:')
%errorbar(1:12, median(allRunShufDiff,1), median(allRunShufDiff,1)-runShufDiff_CI_lower, runShufDiff_CI_upper-median(allRunShufDiff,1),'r:')



figure, hold on
plot([0 1],[0 1],'k:')
plot(allStatDiff(:,trainIdx),allStatShufDiff(:,trainIdx),'k.')
plot(allRunDiff(:,trainIdx),allRunShufDiff(:,trainIdx),'r.')
title(num2str(popSizeVector(ipop)))



figure(999), hold on

vals = {allStatDiff(:,trainIdx), allRunDiff(:,trainIdx)};

hold on
distributionPlot(vals,'histOpt',1,'globalNorm',3,'showMM',6,'divFactor',3,...
    'xValues',[ipop-0.2, ipop+0.2],'color',{'k','r'})
plot([0 4],[0 0],'k:')
ylabel('Relative Performanace')
xlabel('pop size')


end


%%


