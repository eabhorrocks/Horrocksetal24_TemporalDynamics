%% analysing decoding with diff pop sizes


sessionTags = {'M22027', '20220517';...
    'M22029', '20220607';...
    'M22031', '20220524';
    'M22032', '20220621';
    'M22033', '20220706'};

dataDir = 'C:\Users\edward.horrocks\Documents\Code\V1Dynamics\Data\basic_111022';


for isession = 1:size(sessionTags,1)
    fname = [sessionTags{isession,1},'_', sessionTags{isession,2},'_decodePopSize100ms.mat']; %PSTH_noEye %psth3rds

    load(fullfile(dataDir,fname))
    [units.session] = deal(isession);

    s(isession) = session;



end

%% some basic processing


for itime = 1

    for isession = 1:5
        popSize = s(isession).popSize;


        for ipop = 1:numel(popSize)
            popSize(ipop).stat.allPerf = [];
            popSize(ipop).run.allPerf = [];
            popSize(ipop).statShuf.allPerf = [];
            popSize(ipop).runShuf.allPerf = [];

            popSize(ipop).stat.allPerfTimeAv=[];
            popSize(ipop).run.allPerfTimeAv = [];
            popSize(ipop).statShuf.allPerfTimeAv = [];
            popSize(ipop).runShuf.allPerfTimeAv = [];

            for irep = 1:numel(popSize(ipop).rep)

                % average of perms
                popSize(ipop).rep(irep).stat.meanPerf = mean(cat(1,popSize(ipop).rep(irep).stat.perm.meanPerf),1);
                popSize(ipop).rep(irep).run.meanPerf = mean(cat(1,popSize(ipop).rep(irep).run.perm.meanPerf),1);
                popSize(ipop).rep(irep).statShuf.meanPerf = mean(cat(1,popSize(ipop).rep(irep).statShuf.perm.meanPerf),1);
                popSize(ipop).rep(irep).runShuf.meanPerf = mean(cat(1,popSize(ipop).rep(irep).runShuf.perm.meanPerf),1);

                % get aaverrage over timettime
                time_idx = 1:10;
                popSize(ipop).rep(irep).stat.timeAv = mean(popSize(ipop).rep(irep).stat.meanPerf(time_idx));
                popSize(ipop).rep(irep).run.timeAv = mean(popSize(ipop).rep(irep).run.meanPerf(time_idx));
                popSize(ipop).rep(irep).statShuf.timeAv = mean(popSize(ipop).rep(irep).statShuf.meanPerf(time_idx));
                popSize(ipop).rep(irep).runShuf.timeAv = mean(popSize(ipop).rep(irep).runShuf.meanPerf(time_idx));


                % add rep to allPerf feld
                popSize(ipop).stat.allPerf = cat(1, popSize(ipop).stat.allPerf, popSize(ipop).rep(irep).stat.meanPerf);
                popSize(ipop).run.allPerf = cat(1, popSize(ipop).run.allPerf, popSize(ipop).rep(irep).run.meanPerf);
                popSize(ipop).statShuf.allPerf = cat(1, popSize(ipop).statShuf.allPerf, popSize(ipop).rep(irep).statShuf.meanPerf);
                popSize(ipop).runShuf.allPerf = cat(1, popSize(ipop).runShuf.allPerf, popSize(ipop).rep(irep).runShuf.meanPerf);


                % gt all ttime-averagaed perf
                popSize(ipop).stat.allPerfTimeAv = cat(1, popSize(ipop).stat.allPerfTimeAv, popSize(ipop).rep(irep).stat.timeAv);
                popSize(ipop).run.allPerfTimeAv = cat(1, popSize(ipop).run.allPerfTimeAv, popSize(ipop).rep(irep).run.timeAv);
                popSize(ipop).statShuf.allPerfTimeAv = cat(1, popSize(ipop).statShuf.allPerfTimeAv, popSize(ipop).rep(irep).statShuf.timeAv);
                popSize(ipop).runShuf.allPerfTimeAv = cat(1, popSize(ipop).runShuf.allPerfTimeAv, popSize(ipop).rep(irep).runShuf.timeAv);

            end



        end

        s(isession).popSize = popSize;
    end


    %% time-averaged combine across sessions

    popSizeVector = [10 20 40 80 120];

    for ipop=1:5
        ps(ipop).allStat=[];
        ps(ipop).allRun=[];
        ps(ipop).allStatShuf=[];
        ps(ipop).allRunShuf=[];

        for isession = 1:5
            popSize = s(isession).popSize;

            if numel(popSize)>=ipop
                ps(ipop).allStat=cat(1, ps(ipop).allStat,popSize(ipop).stat.allPerfTimeAv);
                ps(ipop).allRun=cat(1, ps(ipop).allRun,popSize(ipop).run.allPerfTimeAv);
                ps(ipop).allStatShuf=cat(1,ps(ipop).allStatShuf, popSize(ipop).statShuf.allPerfTimeAv);
                ps(ipop).allRunShuf=cat(1, ps(ipop).allRunShuf, popSize(ipop).runShuf.allPerfTimeAv);

            end
        end

        ps(ipop).allStatNorm = log2(ps(ipop).allStat./ps(ipop).allStatShuf);
        ps(ipop).allRunNorm = log2(ps(ipop).allRun./ps(ipop).allRunShuf);


    end


    figure, hold on
    errorbar(popSizeVector, cellfun(@mean, {ps.allStat}), cellfun(@sem, {ps.allStat}), 'k')
    errorbar(popSizeVector, cellfun(@mean, {ps.allRun}), cellfun(@sem, {ps.allRun}), 'r')
    errorbar(popSizeVector, cellfun(@mean, {ps.allStatShuf}), cellfun(@sem, {ps.allStatShuf}), 'k:')
    errorbar(popSizeVector, cellfun(@mean, {ps.allRunShuf}), cellfun(@sem, {ps.allRunShuf}), 'r:')
    ylim([0.18, 0.65])
    title('time averaged perf')
    xlabel('pop size')

    % figure, hold on
    % errorbar(popSizeVector, cellfun(@mean, {ps.allStatNorm}), cellfun(@sem, {ps.allStatNorm}), 'k')
    % errorbar(popSizeVector, cellfun(@mean, {ps.allRunNorm}), cellfun(@sem, {ps.allRunNorm}), 'r')
    % errorbar(popSizeVector, cellfun(@mean, {ps.allStatShuf}), cellfun(@sem, {ps.allStatShuf}), 'k:')
    % errorbar(popSizeVector, cellfun(@mean, {ps.allRunShuf}), cellfun(@sem, {ps.allRunShuf}), 'r:')
    % ylim([0.18, 0.65])


    %title([num2str(itime*100-100), ' to ' num2str(itime*100) 'ms'])
    xlim([0 121])

end


%% plot mean performance over ttime for diff popsizes

popSizeVector = [10 20 40 80 120];
figure
for ipop=1:5
    ps(ipop).allStat=[];
    ps(ipop).allRun=[];
    ps(ipop).allStatShuf=[];
    ps(ipop).allRunShuf=[];

    for isession = 1:5
        popSize = s(isession).popSize;

        if numel(popSize)>=ipop
            ps(ipop).allStat=cat(1, ps(ipop).allStat,popSize(ipop).stat.allPerf);
            ps(ipop).allRun=cat(1, ps(ipop).allRun,popSize(ipop).run.allPerf);
            ps(ipop).allStatShuf=cat(1,ps(ipop).allStatShuf, popSize(ipop).statShuf.allPerf);
            ps(ipop).allRunShuf=cat(1, ps(ipop).allRunShuf, popSize(ipop).runShuf.allPerf);

        end
    end

    subplot(1,5,ipop)
    shadedErrorBar(1:10, mean(ps(ipop).allStat), sem(ps(ipop).allStat),'lineProps','k')
    shadedErrorBar(1:10, mean(ps(ipop).allRun), sem(ps(ipop).allRun),'lineProps','r')
    shadedErrorBar(1:10, mean(ps(ipop).allStatShuf), sem(ps(ipop).allStatShuf),'lineProps','k:')
    shadedErrorBar(1:10, mean(ps(ipop).allRunShuf), sem(ps(ipop).allRunShuf),'lineProps','r:')
    title(popSizeVector(ipop))
    ylim([0.1 0.65])
end

%% plot diff in perf bettween stat/run for diff timebins




%% scaattter plots
figure


ipop=4;

% for itime = 1:10
%     subplot(2,5,itime)
%     title(itime)
%     hold on
%     plot(ps(ipop).allStat(:,itime), ps(ipop).allStatShuf(:,itime),'k.')
%     plot(ps(ipop).allRun(:,itime), ps(ipop).allRunShuf(:,itime),'r.')
%     plot([0 1], [0 1],'k:')
% end

% hsotograms
figure
for itime = 1:10
    subplot(2,5,itime)
    title(itime)
    hold on
    histogram(real(log2((ps(ipop).allStatShuf(:,itime)-1/6)./(ps(ipop).allStat(:,itime)-1/6))),'BinEdges',[-1:0.1:1])
    histogram(real(log2((ps(ipop).allRunShuf(:,itime)-1/6)./(ps(ipop).allRun(:,itime)-1/6))),'BinEdges',[-1:0.1:1])
end


%% plot average proportional change in performaance over time (shuf vs normal)
figure
for ipop=1:4;


    statDiff = [];
    runDiff = [];
    statDiff_CI_lower = [];
    statDiff_CI_upper = [];
    runDiff_CI_lower = [];
    runDiff_CI_upper = [];

    for itime = 2:10
        statDiff(:,itime-1) = real(log2((ps(ipop).allStatShuf(:,itime)-1/6)./(ps(ipop).allStat(:,itime)-1/6)));
        runDiff(:,itime-1) = real(log2((ps(ipop).allRunShuf(:,itime)-1/6)./(ps(ipop).allRun(:,itime)-1/6)));

        [statDiff_CI_lower(itime-1), statDiff_CI_upper(itime-1)] = medianCI(statDiff(:,itime-1), 1);
        [runDiff_CI_lower(itime-1), runDiff_CI_upper(itime-1)] = medianCI(runDiff(:,itime-1), 1);

    end

    figure, hold on
    title(num2str(popSizeVector(ipop)))


    errorbar(2:10, median(statDiff,1), median(statDiff,1)-statDiff_CI_lower, statDiff_CI_upper-median(statDiff,1),'k')
    errorbar(2:10, median(runDiff,1), median(runDiff,1)-runDiff_CI_lower, runDiff_CI_upper-median(runDiff,1),'r')
    plot([2 10], [0 0], 'k:')
    ylim([-0.15 0.25])
    ax = gca; ax.XTick = 2:10; ax.XTickLabel = ax.XTick*100-50;
    xlabel('100ms window midpoint (ms)')
    ylabel('log2 performance difference (shuffled/normal)')
    defaultAxesProperties(gca, true)
end


%% try plotting perf diff

%% scatter plot of time-avaeragaed performaance

ipop = 4;

allStat =[];
allRun =[];
allStatShuf =[];
allRunShuf =[];


for isession = 1:5
    popSize = s(isession).popSize;

    if numel(popSize)>=ipop
        for irep = 1:numel(popSize(ipop).rep)
            allStat(end+1) = popSize(ipop).rep(irep).stat.timeAv;
            allRun(end+1) = popSize(ipop).rep(irep).run.timeAv;
            allStatShuf(end+1) = popSize(ipop).rep(irep).statShuf.timeAv;
            allRunShuf(end+1) = popSize(ipop).rep(irep).runShuf.timeAv;
        end
    end

end


figure
subplot(1,3,1), hold on
plot([0 1], [0 1], 'k')
plot(allStat, allRun,'ko')
xlabel('stat'), ylabel('run')

subplot(1,3,2), hold on
plot([0 1], [0 1], 'k')
plot(allStat, allStatShuf,'ko')
xlabel('stat'), ylabel('statShuf')
% ylim([0.2 0.6])
% xlim([0.2 0.6])


subplot(1,3,3), hold on
plot([0 1], [0 1], 'k')
plot(allRun, allRunShuf,'ko')
xlabel('run'), ylabel('runShuf')
% % ylim([0.2 0.6])
% % xlim([0.2 0.6])


%% for each session, compare perf for each pop size

for isession = 1:5
    figure, hold on
    popSize = s(isession).popSize;

    for ipop = 1:numel(popSize)
        errorbar(ipop, mean(popSize(ipop).stat.allPerfTimeAv), sem(popSize(ipop).stat.allPerfTimeAv),'k')
        errorbar(ipop, mean(popSize(ipop).run.allPerfTimeAv), sem(popSize(ipop).run.allPerfTimeAv),'r')

        errorbar(ipop, mean(popSize(ipop).statShuf.allPerfTimeAv), sem(popSize(ipop).statShuf.allPerfTimeAv),'b')
        errorbar(ipop, mean(popSize(ipop).runShuf.allPerfTimeAv), sem(popSize(ipop).runShuf.allPerfTimeAv),'g')

    end
end







%% plot
for isession = 1:5
    popSize = s(isession).popSize;


    figure
    for ipop = 1:numel(popSize)
        subplot(2,5,ipop)
        shadedErrorBar(50:100:950, mean(popSize(ipop).stat.allPerf,1),sem(popSize(ipop).stat.allPerf,1),'lineProps','k')
        shadedErrorBar(50:100:950, mean(popSize(ipop).statShuf.allPerf,1),sem(popSize(ipop).statShuf.allPerf,1),'lineProps','k:')
        %shadedErrorBar(100:100:900, mean(popSize(ipop).run.allPerf,1),sem(popSize(ipop).run.allPerf,1),'lineProps','r')
        %shadedErrorBar(100:100:900, mean(popSize(ipop).runShuf.allPerf,1),sem(popSize(ipop).runShuf.allPerf,1),'lineProps','r:')
        ylim([0 0.8])
    end

    for ipop = 1:numel(popSize)
        subplot(2,5,ipop+5)
        shadedErrorBar(50:100:950, mean(popSize(ipop).run.allPerf,1),sem(popSize(ipop).run.allPerf,1),'lineProps','r')
        shadedErrorBar(50:100:950, mean(popSize(ipop).runShuf.allPerf,1),sem(popSize(ipop).runShuf.allPerf,1),'lineProps','r:')
        ylim([0 0.8])
    end


    % pause
end



%% exaample session looking at mult pop sizes


for isession =2
    figure
    for itime = 1:10
        clear pstat prun pstatshuf prunshuf

        for ipop = 1:numel(s(isession).popSize)

            pstat{ipop} = s(isession).popSize(ipop).stat.allPerf(:,itime);
            prun{ipop} = s(isession).popSize(ipop).run.allPerf(:,itime);

            pstatshuf{ipop} = s(isession).popSize(ipop).statShuf.allPerf(:,itime);
            prunshuf{ipop} = s(isession).popSize(ipop).runShuf.allPerf(:,itime);

        end

        pSizes = [10 20 40 80 120];
        subplot(2,5,itime)
        shadedErrorBar(pSizes(1:numel(pstat)), cellfun(@mean, pstat), cellfun(@sem, pstat))
        shadedErrorBar(pSizes(1:numel(pstat)), cellfun(@mean, prun), cellfun(@sem, prun),'lineProps','r')

        shadedErrorBar(pSizes(1:numel(pstat)), cellfun(@mean, pstatshuf), cellfun(@sem, pstatshuf),'lineProps','k:')
        shadedErrorBar(pSizes(1:numel(pstat)), cellfun(@mean, prunshuf), cellfun(@sem, prunshuf),'lineProps','r:')
        title(['time: ', num2str(itime*100) ' ms'])
        ylim([0 0.7])
        xlim([0 pSizes(end)])
    end
end



%% compare changs in performance with changing pop sizes

% stat
figure
for itime = 1:10;
            clear diffStat diffRun

    for ipop = 1:3

        for isesh = 1:5
            if numel(s(isesh).popSize)>ipop

                diffStat(isesh,ipop) = (mean(s(isesh).popSize(ipop+1).stat.allPerf(:,itime))-1/6)./(mean(s(isesh).popSize(ipop).stat.allPerf(:,itime))-1/6);
                diffRun(isesh,ipop) = (mean(s(isesh).popSize(ipop+1).run.allPerf(:,itime))-1/6)./(mean(s(isesh).popSize(ipop).run.allPerf(:,itime))-1/6);

                diffStatShuf(isesh,ipop) = (mean(s(isesh).popSize(ipop+1).statShuf.allPerf(:,itime))-1/6)./(mean(s(isesh).popSize(ipop).statShuf.allPerf(:,itime))-1/6);
                diffRunShuf(isesh,ipop) = (mean(s(isesh).popSize(ipop+1).runShuf.allPerf(:,itime))-1/6)./(mean(s(isesh).popSize(ipop).runShuf.allPerf(:,itime))-1/6);

            else
                diffStat(isesh,ipop) = nan;
                diffRun(isesh,ipop) = nan;
                diffStatShuf(isesh,ipop) = nan;
                diffRunShuf(isesh,ipop) = nan;
            end

        end
    end
    subplot(2,5,itime), hold on
    %errorbar(1:3, nanmean(diffStat,1),nansem(diffStat,1),'k')
    errorbar(1:3, nanmean(diffRun,1),nansem(diffRun,1),'r')
    % errorbar(1:3, nanmean(diffStatShuf,1),nansem(diffStatShuf,1),'k:')
     errorbar(1:3, nanmean(diffRunShuf,1),nansem(diffRunShuf,1),'r:')
    ylim([1 1.6])
    ax = gca; ax.XTick = 1:3;
    ylabel('Fractional change in perf of increasing pop size')
    title(itime)
defaultAxesProperties(gca, true)

end