%% Analysis of neural tangling

%% load dataa
sessionTags = {'M22027', '20220517';...
    'M22029', '20220607';...
    'M22031', '20220524';
    'M22032', '20220621';
    'M22033', '20220706'};

dataDir = 'C:\Users\edward.horrocks\Documents\Code\V1Dynamics\Data\basic_111022';

s = struct;

for isession = 1:5
    tic
    fname = [sessionTags{isession,1},'_', sessionTags{isession,2},'_FAsmooth_Mar23.mat'];

    load(fullfile(dataDir,fname))

    s(isession).session = session;
end

%% calculate tanagling
clear sesh Qstat Qrun
for isesh = 1:5
    qOpt =s(isesh).session.s.qOpt;
    %qOpt = find(cumsum(s(isesh).session.s.propSharedVariance)>=0.85,1,'first');%6;% s(isesh).session.s.qOpt;
    for ispeed = 1:6
        idim = 1:qOpt;
        % ispeed = 1;
        itime = 1:200;
        statTraj = s(isesh).session.cond(ispeed).meanTrajectory(itime,idim);
        runTraj = s(isesh).session.cond(ispeed+6).meanTrajectory(itime,idim);
        weights = s(isesh).session.s.propSharedVariance(1:qOpt);
        deltaT = 10;
        epsVal = 0.1;

        Qstat=[];
        for it = 2:numel(itime)
            for it2 = 2:numel(itime)

                xt = statTraj(itime(it),:);
                xt2 = statTraj(itime(it2),:);

                xderiv_t = (statTraj(itime(it),:) - statTraj(itime(it-1),:))/deltaT;
                xderiv_t2 = (statTraj(itime(it2),:) - statTraj(itime(it2-1),:))/deltaT;

                Qstat(it,it2) = (norm(xderiv_t-xderiv_t2)^2) / (norm(xt-xt2)^2 + epsVal);

            end
        end



        Qrun=[];
        for it = 2:numel(itime)
            for it2 = 2:numel(itime)

                xt = runTraj(itime(it),:);
                xt2 = runTraj(itime(it2),:);
                %
                xderiv_t = (runTraj(itime(it),:) - runTraj(itime(it-1),:))/deltaT;
                xderiv_t2 = (runTraj(itime(it2),:) - runTraj(itime(it2-1),:))/deltaT;

                Qrun(it,it2) = (norm(xderiv_t-xderiv_t2)^2) / (norm(xt-xt2)^2 + epsVal);

            end
        end

        plotFlag=false;
        if plotFlag
            figure,
            subplot(131)
            imagesc(Qstat)
            subplot(132)
            imagesc(Qrun)

            maxVal =max(cat(1,Qstat(:), Qrun(:)));

            ax1 = subplot(131), hold on
            caxis([0 maxVal])
            plot([1 200], [21 21], 'r')
            plot([1 200], [121 121], 'r')
            plot([21 21], [1, 200], 'r')
            plot([121 121],[1, 200], 'r')
            colorbar
            ax = gca
            ax.XTick = 0:10:200
            ax.YTick = 0:10:200
            ax.XTickLabel = -200:100:1800;
            ax.YTickLabel = -200:100:1800;
            grid on

            ax2 = subplot(132), hold on
            plot([1 200], [21 21], 'r')
            plot([1 200], [121 121], 'r')
            plot([21 21], [1, 200], 'r')
            plot([121 121],[1, 200], 'r')
            caxis([0 maxVal])
            colormap(1-gray)
            colorbar
            ax = gca
            ax.XTick = 0:10:200
            ax.YTick = 0:10:200
            ax.XTickLabel = -200:100:1800;
            ax.YTickLabel = -200:100:1800;
            grid on

            ax3 = subplot(133)
            imagesc(Qrun-Qstat)
        end

        sesh(isesh).speed(ispeed).Qstat = Qstat;
        sesh(isesh).speed(ispeed).Qrun = Qrun;

    end
end



%% get temporally local tangling

binSize = 10;
windowSize = 200/binSize;

for isesh = 1:5
    for ispeed = 1:6
        this_Qstat = sesh(isesh).speed(ispeed).Qstat;
        this_Qrun = sesh(isesh).speed(ispeed).Qrun;
        for itime = 1:180
            tqs = this_Qstat(itime:itime+windowSize, itime:itime+windowSize);
            tqr = this_Qrun(itime:itime+windowSize, itime:itime+windowSize);


            sesh(isesh).speed(ispeed).statLocalTangling(itime) = ...
                prctile(tqs(1+windowSize/2,:),90,'all');

            sesh(isesh).speed(ispeed).runLocalTangling(itime) = ...
                prctile(tqr(1+windowSize/2,:),90,'all');
        end
    end
end

%% plot overall average local tangling
allStat=[];
allRun=[];
for isesh=1:5
    allStat = cat(1,allStat,sesh(isesh).speed.statLocalTangling);
    allRun = cat(1,allRun,sesh(isesh).speed.runLocalTangling);
end

figure, hold on
shadedErrorBar(-95:10:1695, mean(allStat,1), sem(allStat,1),'lineProps','k')
shadedErrorBar(-95:10:1695, mean(allRun,1), sem(allRun,1),'lineProps','r')
xlim([-200 1800])
ax=gca; ax.XTick = -200:200:1800;
defaultAxesProperties(gca, true)
ylabel('Local Tangling')
xlabel('Time (ms)')

figure

% plot by speed
subplot(211)
for ispeed = 1:6
    allStat=[];
    allRun=[];
    for isesh=1:5
        allStat = cat(1,allStat,sesh(isesh).speed(ispeed).statLocalTangling);
        allRun = cat(1,allRun,sesh(isesh).speed(ispeed).runLocalTangling);
    end
    subplot(211), hold on
    plot(-95:10:1695, mean(allStat,1),'Color', speedcols(ispeed,:))
    subplot(212), hold on
    plot(-95:10:1695, mean(allRun,1),'Color', speedcols(ispeed,:))
end

subplot(211)
xlim([-200 1800])
ylim([0 0.008])
defaultAxesProperties(gca, true)

subplot(212)
xlim([-200 1800])
ylim([0 0.008])
defaultAxesProperties(gca, true)




%% RM-anova for tangling

time = 1:180;
statVals = allStat';
runVals = allRun';
allVals = cat(1,statVals(:),runVals(:));

timeVec = categorical(repmat(time,1,size(statVals,2)*2)');
subjVec = categorical([repelem(1:5, 1, 6*180)]); subjVec=[subjVec,subjVec]';
speedVec = categorical([repmat(1:6,1,5*180), repmat(1:6,1,5*180)]');
stateVec = categorical(repelem(1:2,1,numel(statVals))');

% [p,tbl,stats,terms] = anovan(allVals,{timeVec,stateVec,subjVec},'model',2,'random',3,'varnames',{'Time','State','Subj'});


[p,tbl,stats,terms] = anovan(allVals,{timeVec,stateVec,subjVec,speedVec},'model','full','varnames',{'Time','State','Subj','Speed'});


%% scatter of mean tangling

allStat_mean = mean(allStat(:,21:180),2);
allRun_mean = mean(allRun(:,21:180),2);

figure, hold on
plot(allStat_mean, allRun_mean,'ko')

xlim([0 0.0035])
ylim([0 0.0035])
plot([0 0.0035],[0 0.0035],'k:')

figure, hold on
plot([0 0.006],[0 0.006],'k:')

for ispeed = 1:6
    allStat=[];
    allRun=[];
    for isesh=1:5
        allStat = cat(1,allStat,sesh(isesh).speed(ispeed).statLocalTangling);
        allRun = cat(1,allRun,sesh(isesh).speed(ispeed).runLocalTangling);
    end
    allStat_mean = mean(allStat(:,21:71),2);
    allRun_mean = mean(allRun(:,21:71),2);

    plot(allStat_mean, allRun_mean,'o','MarkerEdgeColor','none','MarkerFaceColor',speedcols(ispeed,:))

end

xlim([0 0.006])
ylim([0 0.006])