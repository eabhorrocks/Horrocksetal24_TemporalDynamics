%% analyse FA - fit separately

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
    fname = [sessionTags{isession,1},'_', sessionTags{isession,2},'_FA_fitSep.mat'];

    load(fullfile(dataDir,fname))

    s(isession).session = session;
end


%% recreate cond

for isesh=1:5
for icond = 1:6
    s(isesh).session.s.cond(icond) = s(isesh).session.stat.s.cond(icond);
end
for icond = 7:12
    s(isesh).session.s.cond(icond) = s(isesh).session.run.s.cond(icond-6);
end
end

%% calculate tangling
clear sesh Qstat Qrun
for isesh = 1:5
    qOpt_stat = s(isesh).session.stat.s.qOpt;
    qOpt_run = s(isesh).session.run.s.qOpt;
    qOpt = min([qOpt_stat, qOpt_run]);
    %qOpt = find(cumsum(s(isesh).session.s.propSharedVariance)>=0.85,1,'first');%6;% s(isesh).session.s.qOpt;
    for ispeed = 1:6
        % ispeed = 1;
        itime = 1:200;
        statTraj = s(isesh).session.s.cond(ispeed).meanTrajectory(itime,1:qOpt);
        runTraj = s(isesh).session.s.cond(ispeed+6).meanTrajectory(itime,1:qOpt);
%         weights = s(isesh).session.s.propSharedVariance(1:qOpt);
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

speedcols = inferno(6);
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
ylim([0 0.004])
defaultAxesProperties(gca, true)

subplot(212)
xlim([-200 1800])
ylim([0 0.004])
defaultAxesProperties(gca, true)


%  plot speed-wise local tangling
figure

for ispeed = 1:6
    allStat =[];
    allRun=[];
    for isesh = 1:5
        allStat = cat(1,allStat,sesh(isesh).speed(ispeed).statLocalTangling);
        allRun = cat(1,allRun,sesh(isesh).speed(ispeed).runLocalTangling);
    end
    subplot(2,3,ispeed), hold on
    shadedErrorBar(-95:10:1695, mean(allStat,1), sem(allStat,1))
    shadedErrorBar(-95:10:1695, mean(allRun,1), sem(allRun,1),...
        'lineProps','r')

end



%% RM-anova for tangling

% time = 1:180;
% statVals = allStat';
% runVals = allRun';
% allVals = cat(1,statVals(:),runVals(:));
% 
% timeVec = categorical(repmat(time,1,size(statVals,2)*2)');
% subjVec = categorical([repelem(1:5, 1, 6*180)]); subjVec=[subjVec,subjVec]';
% speedVec = categorical([repmat(1:6,1,5*180), repmat(1:6,1,5*180)]');
% stateVec = categorical(repelem(1:2,1,numel(statVals))');
% 
% % [p,tbl,stats,terms] = anovan(allVals,{timeVec,stateVec,subjVec},'model',2,'random',3,'varnames',{'Time','State','Subj'});
% 
% 
% [p,tbl,stats,terms] = anovan(allVals,{timeVec,stateVec,subjVec,speedVec},'model','full','varnames',{'Time','State','Subj','Speed'});


%% scatter of mean tangling

allStat_mean = mean(allStat(:,21:180),2);
allRun_mean = mean(allRun(:,21:180),2);

figure, hold on
plot(allStat_mean, allRun_mean,'ko')


plot([0 0.0035],[0 0.0035],'k:')
xlim([0 0.0005])
ylim([0 0.0005])

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

xlim([0 0.002])
ylim([0 0.002])


%% Speed and acceleration of tarjectories

allStat=[];
allRun=[];

allStatAcc = [];
allRunAcc = [];

allStatJerk = [];
allRunJerk = [];


for isesh = 1:5
%     weights = s(isesh).session.s.propSharedVariance';
    qOpt_stat = s(isesh).session.stat.s.qOpt;
    qOpt_run = s(isesh).session.run.s.qOpt;
    qOpt = min([qOpt_stat, qOpt_run]);
    weights = repelem(1, 1, qOpt);
    nFactors = qOpt;

    for icond = 1:12
        for itime = 1:200-1
            s(isesh).session.s.cond(icond).speed(itime) = ....
                (sum(abs(s(isesh).session.s.cond(icond).meanTrajectory(itime+1,1:qOpt) - ...
                s(isesh).session.s.cond(icond).meanTrajectory(itime,1:qOpt)) .* weights) *(1000/10)) / nFactors;
        end

        s(isesh).session.s.cond(icond).maxSpeedOnset = max(s(isesh).session.s.cond(icond).speed(1:61));
        s(isesh).session.s.cond(icond).maxSpeedOffset = max(s(isesh).session.s.cond(icond).speed(100:161));


        s(isesh).session.s.cond(icond).speedAcc = diff( s(isesh).session.s.cond(icond).speed);
        [s(isesh).session.s.cond(icond).accPower,f] = pspectrum(s(isesh).session.s.cond(icond).speedAcc,100);
        s(isesh).session.s.cond(icond).speedJerk = diff(s(isesh).session.s.cond(icond).speedAcc);
    end
    s(isesh).session.s.meanStatSpeed = cat(1,s(isesh).session.s.cond(1:6).speed);
    s(isesh).session.s.meanRunSpeed = cat(1,s(isesh).session.s.cond(7:12).speed);

    s(isesh).session.s.meanStatSpeedAcc = cat(1,s(isesh).session.s.cond(1:6).speedAcc);
    s(isesh).session.s.meanRunSpeedAcc = cat(1,s(isesh).session.s.cond(7:12).speedAcc);



    s(isesh).session.s.meanStatJerk = cat(1,s(isesh).session.s.cond(1:6).speedJerk);
    s(isesh).session.s.meanRunJerk = cat(1,s(isesh).session.s.cond(7:12).speedJerk);


    allStat = cat(1,allStat, s(isesh).session.s.meanStatSpeed);
    allRun = cat(1, allRun, s(isesh).session.s.meanRunSpeed);

    allStatAcc = cat(1, allStatAcc, s(isesh).session.s.meanStatSpeedAcc);
    allRunAcc = cat(1, allRunAcc, s(isesh).session.s.meanRunSpeedAcc);

    allStatJerk = cat(1, allStatJerk, s(isesh).session.s.meanStatJerk);
    allRunJerk = cat(1, allRunJerk, s(isesh).session.s.meanRunJerk);
end

figure, hold on
shadedErrorBar(-190:10:1795, mean(allRun,1), sem(allRun,1),'lineProps', 'r')
shadedErrorBar(-190:10:1795, mean(allStat,1), sem(allStat,1))

ylim([0 4])
xlim([-200 1800])
ax=gca; ax.XTick = -200:200:1800;
defaultAxesProperties(gca, true)
title('Speed')

figure, hold on
shadedErrorBar(-180:10:1795, mean(allRunAcc.*100,1), sem(allRunAcc.*100,1),'lineProps', 'r')
shadedErrorBar(-180:10:1795, mean(allStatAcc.*100,1), sem(allStatAcc.*100,1))

% ylim([0 500])
xlim([-200 1800])
ax=gca; ax.XTick = -200:200:1800;
defaultAxesProperties(gca, true)
title('Acceleration')



%% stats for speed

speedVec=[];
seshVec=[];
onsetMaxSpdVec=[];
offsetMaxSpdVec=[];
stateVec=[];

for isesh = 1:5
    for ispeed = 1:6
        speedVec=cat(1,speedVec,ispeed);
        seshVec=cat(1,seshVec,isesh);
        stateVec=cat(1,stateVec,1);
        onsetMaxSpdVec=cat(1,onsetMaxSpdVec,s(isesh).session.s.cond(ispeed).maxSpeedOnset);
        offsetMaxSpdVec=cat(1,offsetMaxSpdVec,s(isesh).session.s.cond(ispeed).maxSpeedOffset);
    end
end



for isesh = 1:5
    for ispeed = 7:12
        speedVec=cat(1,speedVec,ispeed-6);
        seshVec=cat(1,seshVec,isesh);
        stateVec=cat(1,stateVec,2);
        onsetMaxSpdVec=cat(1,onsetMaxSpdVec,s(isesh).session.s.cond(ispeed).maxSpeedOnset);
        offsetMaxSpdVec=cat(1,offsetMaxSpdVec,s(isesh).session.s.cond(ispeed).maxSpeedOffset);
    end
end


valsVec = onsetMaxSpdVec;

tbl = table(valsVec,speedVec,stateVec,seshVec,...
    'VariableNames',{'vals','speed','state','sesh'});

f = 'vals ~ state + (1|speed) + (1|sesh)';

lme = fitlme(tbl, f, 'DummyVarCoding', 'reference')


[mean(valsVec(1:30)), sem(valsVec(1:30)); mean(valsVec(31:60)), sem(valsVec(31:60))]



%% angle of trajectory vs sustained period


for isession = 1:5
    qOpt_stat = s(isession).session.stat.s.qOpt;
    qOpt_run = s(isession).session.run.s.qOpt;
    qOpt = min([qOpt_stat, qOpt_run]);

    for icond = 1:12
        bin_N = 1; % bin every N elements
        binSize = 10; % bin size in ms
        binVector = -200:binSize*bin_N:(1800-binSize*bin_N);
        timeBinVector = round(binVector,2); % time associated with each index (1st edge of bin)
        nWins = 10; % number of consecutive windows required to reach steady state

        % set parameters
        % weights = s.propSharedVariance(:)'; % weights are frac. of shared variance explained by each component
        weights = repelem(1, 1, qOpt);

        % onset trajectories
        initialSSIdx = 1:find(timeBinVector==0)-1;
        finalSSIdx = find(timeBinVector==500):find(timeBinVector==1000)-1;
        finalBaseIdx = find(timeBinVector==1500):find(timeBinVector==1790)-1;


        traj = s(isession).session.s.cond(icond).meanTrajectory(:, 1:qOpt);


        ss_mean = mean(traj(finalSSIdx,:),1);
        ini_mean = mean(traj(initialSSIdx,:),1);
        final_mean = mean(traj(finalBaseIdx,:),1);

        refAxis = ss_mean-ini_mean;
        refAxisOffset = final_mean-ss_mean;


        % onset
        angle=[];
        angleSub=[];
        angleDir=[];
        v1full=[];
        for itime = 1:200
            v1full(itime,:) = ss_mean-traj(itime,:);
            v1 = ss_mean-traj(itime,:); v1=v1(:);
            v2 = refAxis; v2=v2(:);
            angle(itime) = dot(v1 / norm(v1), v2 / norm(v2));
            distFromSS(itime) = sum(abs(traj(itime,:)-ini_mean));
            %      angleSub(itime) = subspace(v1,v2);


        end

        % offset
        angleOffset=[];
        distFromSSOffset=[];

        for itime = 1:200
            v1full(itime,:) = final_mean-traj(itime,:);
            v1 = final_mean-traj(itime,:); v1=v1(:);
            v2 = refAxisOffset; v2=v2(:);
            angleOffset(itime) = dot(v1 / norm(v1), v2 / norm(v2));
            distFromSSOffset(itime) = sum(abs(traj(itime,:)-ini_mean));
            %      angleSub(itime) = subspace(v1,v2);


        end


        s(isession).session.cond(icond).angle = acosd(angle);%mod(acosd(angle).*sign(angle),360);
        s(isession).session.cond(icond).angleOffset = acosd(angleOffset);%mod(acosd(angle).*sign(angle),360);

        s(isession).session.cond(icond).distFromStim = distFromSS;
        s(isession).session.cond(icond).distFromStim = distFromSSOffset;

        s(isession).session.cond(icond).sumDiffAngle = sum(abs(diff(s(isession).session.cond(icond).angle(1:111))));
        s(isession).session.cond(icond).sumDiffAngleOffset = sum(abs(diff(s(isession).session.cond(icond).angleOffset(111:end))));

        s(isession).session.cond(icond).v1full = v1full;
        % s(isession).session.cond(icond).angleSub = angleSub;



    end
end



%% sum change in angle

speedVec=[];
seshVec=[];
onsetSumAngle=[];
offsetSumAngle=[];
stateVec=[];

for isesh = 1:5
    for ispeed = 1:6
        speedVec=cat(1,speedVec,ispeed);
        seshVec=cat(1,seshVec,isesh);
        stateVec=cat(1,stateVec,1);
        onsetSumAngle=cat(1,onsetSumAngle,s(isesh).session.cond(ispeed).sumDiffAngle);
        offsetSumAngle=cat(1,offsetSumAngle,s(isesh).session.cond(ispeed).sumDiffAngleOffset);
    end
end



for isesh = 1:5
    for ispeed = 7:12
        speedVec=cat(1,speedVec,ispeed-6);
        seshVec=cat(1,seshVec,isesh);
        stateVec=cat(1,stateVec,2);
        onsetSumAngle=cat(1,onsetSumAngle,s(isesh).session.cond(ispeed).sumDiffAngle);
        offsetSumAngle=cat(1,offsetSumAngle,s(isesh).session.cond(ispeed).sumDiffAngleOffset);
    end
end


valsVec = offsetSumAngle;

tbl = table(valsVec,speedVec,stateVec,seshVec,...
    'VariableNames',{'vals','speed','state','sesh'});

    f = 'vals ~ state + (1|speed) + (1|sesh)';
      
    lme = fitlme(tbl, f, 'DummyVarCoding', 'reference')


    [mean(valsVec(1:30)), sem(valsVec(1:30)); mean(valsVec(31:60)), sem(valsVec(31:60))]


    
%% plot overall change in angle
figure
subplot(121), hold on
for ispeed = 1:6
    for isession = 1:5

        plot(s(isession).session.cond(ispeed).sumDiffAngle, s(isession).session.cond(ispeed+6).sumDiffAngle,'o', 'MarkerEdgeColor','none','MarkerFaceColor',speedcols(ispeed,:));
    end
end

subplot(122), hold on
for ispeed = 1:6
    for isession = 1:5

        plot(s(isession).session.cond(ispeed).sumDiffAngleOffset, s(isession).session.cond(ispeed+6).sumDiffAngleOffset,'o', 'MarkerEdgeColor','none','MarkerFaceColor',speedcols(ispeed,:));
    end
end



subplot(121)
plot([0 600],[0 600],'k:')
xlim([0 600])
ylim([0 600])
axis equal
defaultAxesProperties(gca, true)

subplot(122)
plot([0 600],[0 600],'k:')
xlim([0 600])
ylim([0 600])
axis equal
defaultAxesProperties(gca, true)

%% plot average angle over time

allStat=[];
allRun=[];

for isession = 1:5
    for ispeed = 1:6
        allStat=cat(1,allStat, s(isession).session.cond(ispeed).angle);
        allRun=cat(1,allRun, s(isession).session.cond(ispeed+6).angle);
    end
end

figure
shadedErrorBar(-195:10:1795,mean(allStat,1), sem(allStat,1), 'lineProps','k')
shadedErrorBar(-195:10:1795,mean(allRun,1), sem(allRun,1), 'lineProps','r')
defaultAxesProperties(gca, true)
title('Onset')

allStat=[];
allRun=[];

for isession = 1:5
    for ispeed = 1:6
        allStat=cat(1,allStat, s(isession).session.cond(ispeed).angleOffset);
        allRun=cat(1,allRun, s(isession).session.cond(ispeed+6).angleOffset);
    end
end

figure
shadedErrorBar(-195:10:1795,mean(allStat,1), sem(allStat,1), 'lineProps','k')
shadedErrorBar(-195:10:1795,mean(allRun,1), sem(allRun,1), 'lineProps','r')
ylim([0 140])
defaultAxesProperties(gca, true)
title('Offset')



figure
% plot by speed
for ispeed = 1:6
    allStat=[];
    allRun=[];
for isession=1:5
    allStat = cat(1,allStat,s(isession).session.cond(ispeed).angle);
    allRun = cat(1,allRun,s(isession).session.cond(ispeed+6).angle);
end
subplot(221), hold on
plot(-195:10:1795, mean(allStat,1),'Color', speedcols(ispeed,:))
subplot(223), hold on
plot(-195:10:1795, mean(allRun,1),'Color', speedcols(ispeed,:))
end



%% plot by speed
for ispeed = 1:6
    allStat=[];
    allRun=[];
for isession=1:5
    allStat = cat(1,allStat,s(isession).session.cond(ispeed).angle);
    allRun = cat(1,allRun,s(isession).session.cond(ispeed+6).angle);
end
subplot(221), hold on
plot(-195:10:1795, mean(allStat,1),'Color', speedcols(ispeed,:))
subplot(223), hold on
plot(-195:10:1795, mean(allRun,1),'Color', speedcols(ispeed,:))
end

subplot(221)
xlim([-200 1800])
ylim([0 160])
defaultAxesProperties(gca, true)

subplot(223)
xlim([-200 1800])
ylim([0 160])
defaultAxesProperties(gca, true)

subplot(222)
xlim([-200 1800])
ylim([0 160])
defaultAxesProperties(gca, true)

subplot(224)
xlim([-200 1800])
ylim([0 160])
defaultAxesProperties(gca, true)





%% get distance metrics with no weights

bin_N = 1; % bin every N elements
binSize = 10; % bin size in ms
binVector = -200:binSize*bin_N:(1800-binSize*bin_N);
timeBinVector = round(binVector,2); % time associated with each index (1st edge of bin)
nWins = 10; % number of consecutive windows required to reach steady state

for isession = 1:5

    % set parameters
    % weights = s.propSharedVariance(:)'; % weights are frac. of shared variance explained by each component
        qOpt_stat = s(isession).session.stat.s.qOpt;
    qOpt_run = s(isession).session.run.s.qOpt;
    qOpt = min([qOpt_stat, qOpt_run]);
    weights = repelem(1, 1, qOpt);
    cond = s(isession).session.s.cond;

    % onset trajectories
    initialSSIdx = 1:find(timeBinVector==0)-1;
    finalSSIdx = find(timeBinVector==500):find(timeBinVector==1000)-1;



    for icond = 1:numel(cond)
        neuralTrajectory = cond(icond).meanTrajectory(:,1:qOpt);
        output_onset(icond) = calcPopulationGeoMetrics(neuralTrajectory, weights, timeBinVector, initialSSIdx, finalSSIdx, nWins);
    end

    % offset trajectories
    initialSSIdx = find(timeBinVector==500):find(timeBinVector==1000)-1;
    finalSSIdx = find(timeBinVector==1500):find(timeBinVector==1800-binSize*bin_N);

    for icond = 1:numel(cond)
        neuralTrajectory = cond(icond).meanTrajectory(:,1:qOpt);
        output_offset(icond) = calcPopulationGeoMetrics(neuralTrajectory, weights, timeBinVector, initialSSIdx, finalSSIdx, nWins);
    end


    fnames = fieldnames(output_onset);
    for ifield = 1:numel(fnames)
        s(isession).session.stat.(['onset_', fnames{ifield}]) = {output_onset(1:6).(fnames{ifield})};
        s(isession).session.run.(['onset_', fnames{ifield}]) = {output_onset(7:12).(fnames{ifield})};
    end

    fnames = fieldnames(output_offset);
    for ifield = 1:numel(fnames)
        s(isession).session.stat.(['offset_', fnames{ifield}]) = {output_offset(1:6).(fnames{ifield})};
        s(isession).session.run.(['offset_', fnames{ifield}]) = {output_offset(7:12).(fnames{ifield})};
    end


end


%% get standard geometric values

statDR = [];
runDR  = [];

statDD = [];
runDD = [];

statCD = [];
runCD = [];

statDRoff = [];
runDRoff  = [];

statDDoff = [];
runDDoff = [];

statCDoff = [];
runCDoff = [];

statPerf = [];
runPerf = [];


for isession = 1:5

    statDR  =cat(1, statDR, s(isession).session.stat.onset_distanceRatio);
    runDR  =cat(1, runDR, s(isession).session.run.onset_distanceRatio);

    statDD  =cat(1, statDD, s(isession).session.stat.onset_directDistance);
    runDD  =cat(1, runDD, s(isession).session.run.onset_directDistance);

    statCD  =cat(1, statCD, s(isession).session.stat.onset_cumulativeDistanceTravelled);
    runCD  =cat(1, runCD, s(isession).session.run.onset_cumulativeDistanceTravelled);

    statDRoff  =cat(1, statDRoff, s(isession).session.stat.offset_distanceRatio);
    runDRoff  =cat(1, runDRoff, s(isession).session.run.offset_distanceRatio);

    statDDoff  =cat(1, statDDoff, s(isession).session.stat.offset_directDistance);
    runDDoff  =cat(1, runDDoff, s(isession).session.run.offset_directDistance);

    statCDoff  =cat(1, statCDoff, s(isession).session.stat.offset_cumulativeDistanceTravelled);
    runCDoff  =cat(1, runCDoff, s(isession).session.run.offset_cumulativeDistanceTravelled);
% 
%     statPerf  =cat(1, statPerf, s(isession).session.stat.meanPerf);
%     runPerf  =cat(1, runPerf, s(isession).session.run.meanPerf);

end


statDR = cell2mat(statDR);
runDR  = cell2mat(runDR);

statDD = cell2mat(statDD);
runDD = cell2mat(runDD);

statCD = cell2mat(statCD);
runCD = cell2mat(runCD);

statDRoff = cell2mat(statDRoff);
runDRoff  = cell2mat(runDRoff);

statDDoff = cell2mat(statDDoff);
runDDoff = cell2mat(runDDoff);

statCDoff = cell2mat(statCDoff);
runCDoff = cell2mat(runCDoff);

%% plot scatter plots of trajectory distance measures

speedcols = inferno(6);
areacols = tab10(5);
figure
subplot(131), hold on
plot([0 16], [0 16], 'k')
for ispeed = 1:6
    plot(statDD(:,ispeed), runDD(:,ispeed), 'o', 'MarkerFaceColor', speedcols(ispeed,:), 'MarkerEdgeColor', 'none')
end
axis equal

title('Direct distance')
xlabel('stat'), ylabel('run')
xlim([0 16]), ylim([0 16])
ax = gca; ax.XTick = 0:16; ax.YTick = 0:16;
defaultAxesProperties(gca, true)


subplot(132), hold on
plot([0 42], [0 42], 'k')
for ispeed = 1:6
    plot(statCD(:,ispeed), runCD(:,ispeed), 'o', 'MarkerFaceColor', speedcols(ispeed,:), 'MarkerEdgeColor', 'none')
end
% for isession = 1:5
%     plot(mean(statCD(isession,:)), mean(runCD(isession,:)), 'o', 'MarkerFaceColor', areacols(isession,:), 'MarkerEdgeColor', 'w')
% end
axis equal

title('Distance travelled')
xlabel('stat'), ylabel('run')
xlim([0 42]), ylim([0 42])
ax = gca; ax.XTick = 0:5:40; ax.YTick = 0:5:40;
defaultAxesProperties(gca, true)


subplot(133), hold on
plot([0 5], [0 5], 'k')
for ispeed = 1:6
    plot(log2(statDR(:,ispeed)), log2(runDR(:,ispeed)), 'o', 'MarkerFaceColor', speedcols(ispeed,:), 'MarkerEdgeColor', 'none')
end
% for isession = 1:5
%     plot(mean(statDR(isession,:)), mean(runDR(isession,:)), 'o', 'MarkerFaceColor', areacols(isession,:), 'MarkerEdgeColor', 'w')
% end
axis equal

title('Distance ratio')
xlabel('stat'), ylabel('run')
xlim([0 5]), ylim([0 5])
ax = gca; ax.XTick = log2([1 2 4 8 16 32]); ax.YTick = ax.XTick;
ax.XTickLabel = [1 2 4 8 16 32]; ax.YTickLabel = ax.XTickLabel;
defaultAxesProperties(gca, true)



%% offset period

speedcols = inferno(6);
areacols = tab10(5);
figure
subplot(131), hold on
plot([0 16], [0 16], 'k')
for ispeed = 1:6
    plot(statDDoff(:,ispeed), runDDoff(:,ispeed), 'o', 'MarkerFaceColor', speedcols(ispeed,:), 'MarkerEdgeColor', 'none')
end
axis equal

title('Direct distance')
xlabel('stat'), ylabel('run')
xlim([0 16]), ylim([0 16])
ax = gca; ax.XTick = 0:16; ax.YTick = 0:16;
defaultAxesProperties(gca, true)


subplot(132), hold on
plot([0 42], [0 42], 'k')
for ispeed = 1:6
    plot(statCDoff(:,ispeed), runCDoff(:,ispeed), 'o', 'MarkerFaceColor', speedcols(ispeed,:), 'MarkerEdgeColor', 'none')
end
% for isession = 1:5
%     plot(mean(statCD(isession,:)), mean(runCD(isession,:)), 'o', 'MarkerFaceColor', areacols(isession,:), 'MarkerEdgeColor', 'w')
% end
axis equal

title('Distance travelled')
xlabel('stat'), ylabel('run')
xlim([0 42]), ylim([0 42])
ax = gca; ax.XTick = 0:5:40; ax.YTick = 0:5:40;
defaultAxesProperties(gca, true)


subplot(133), hold on
plot([0 5], [0 5], 'k')
for ispeed = 1:6
    plot(log2(statDRoff(:,ispeed)), log2(runDRoff(:,ispeed)), 'o', 'MarkerFaceColor', speedcols(ispeed,:), 'MarkerEdgeColor', 'none')
end
% for isession = 1:5
%     plot(mean(statDR(isession,:)), mean(runDR(isession,:)), 'o', 'MarkerFaceColor', areacols(isession,:), 'MarkerEdgeColor', 'w')
% end
axis equal

title('Distance ratio')
xlabel('stat'), ylabel('run')
xlim([0 5]), ylim([0 5])
ax = gca; ax.XTick = log2([1 2 4 8 16 32]); ax.YTick = ax.XTick;
ax.XTickLabel = [1 2 4 8 16 32]; ax.YTickLabel = ax.XTickLabel;
defaultAxesProperties(gca, true)
