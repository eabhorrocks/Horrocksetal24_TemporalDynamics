%% plot correlation analyses

sessionTags = {'M22027', '20220517';...
    'M22029', '20220607';...
    'M22031', '20220524';
    'M22032', '20220621';
    'M22033', '20220706'};

dataDir = 'E:\V1Data\Data\v1_fromC24';


for isession = 1:size(sessionTags,1)
    isession
    fname = [sessionTags{isession,1},'_', sessionTags{isession,2},'_','corrs_longAxis.mat']; %_FAcorrs_dev.mat

    load(fullfile(dataDir,fname))


    session1(isession).stat = session.stat;
    session1(isession).run = session.run;
    session1(isession).units = session.gU;
    
% cell-type classification
%     for iunit = 1:numel(session1(isession).units)
%         [ccg, t] = CCG(session1(isession).units(iunit).spike_times, ones(size(session1(isession).units(iunit).spike_times)),...
%             'binSize', 0.0005, 'duration', 0.1,'norm', 'rate');
%         %         goodUnits(iunit).acg = ccg;
%         fit_params_out = fit_ACG(ccg,false);
% 
%         session1(isession).units(iunit).tau_rise = fit_params_out.acg_tau_rise;
%     end
% 
%     narrow_idx = find([session1(isession).units.duration]<=0.45);
%     wide_idx = find([session1(isession).units.duration]>0.45 & [session1(isession).units.tau_rise]>6);
%     pyr_idx = find(~ismember(1:numel(session1(isession).units), [narrow_idx,wide_idx]));
% 
%     [session1(isession).units.cellType]=deal(nan);
%     [session1(isession).units(pyr_idx).cellType]=deal(1);
%     [session1(isession).units(narrow_idx).cellType]=deal(2);
%     [session1(isession).units(wide_idx).cellType]=deal(3);



    %     session1(isession).stat.tl_noiseCorrArray = session.stat.tl_noiseCorrArray;
    %     session1(isession).run.tl_noiseCorrArray = session.run.tl_noiseCorrArray;
    %
    %     clear session
    %

end


session = session1;


%% plot correlations for each speed

speeds = [0 16 32 64 128 256];
nTimeBins = numel(1:200-19);
speedcols = inferno(6);


figure
for istim = 1:6
    allStat =[];
    allRun= [];


    for iint = 1:nTimeBins

        statTemp=[];
        runTemp=[];


        %subplot(1,6,istim), hold on, title(speeds(istim))
        for isession = 1:5

            temp = session(isession).stat.stim(istim).tca(:,:,iint);
            temp = setUpperTri2NaN(temp);

            statTemp = cat(1, statTemp, temp(:));

            temp = session(isession).run.stim(istim).tca(:,:,iint);
            temp = setUpperTri2NaN(temp);
            runTemp = cat(1, runTemp, temp(:));

        end

        allStat = cat(2,allStat,statTemp);
        allRun = cat(2,allRun, runTemp);

    end

    subplot(121), hold on,
    shadedErrorBar(-100:10:1700, mean(allStat,1,'omitnan'),nansem(allStat,1),'lineProps',{'Color', speedcols(istim,:)})
    subplot(122), hold on
    shadedErrorBar(-100:10:1700, mean(allRun,1,'omitnan'),nansem(allRun,1),'lineProps',{'Color', speedcols(istim,:)});


    %shadedErrorBar(-100:10:1700, mean(allStat,1,'omitnan'),sem(allStat,1),'lineProps','k')
    %shadedErrorBar(-100:10:1700, mean(allRun,1,'omitnan'),sem(allRun,1),'lineProps','r')
    %ax = gca;ax.XTick = -200:200:1800;
    %defaultAxesProperties(gca, true)
    %ylim([-0.02 0.27])
end

subplot(121)
ax=gca; ax.XTick = -200:200:1800; ylim([0 0.2]); xlim([-200 1800])
defaultAxesProperties(gca,true)

subplot(122)
ax=gca; ax.XTick = -200:200:1800; ylim([0 0.2]); xlim([-200 1800])
defaultAxesProperties(gca,true)

%% plot pair-average of total, signal and noise correlations

sesh2include=1:5;


nTimeBins = numel(1:200-19);

%%%%%%% TOTAL $$$$$$$
figure(123)

allStat = [];
allRun = [];


for iint = 1:nTimeBins
    tempStat = [];
    tempRun = [];
    for isession = sesh2include
        temp = session(isession).stat.totalCorrArray(:,:,iint);
        temp = setUpperTri2NaN(temp);

        tempStat = cat(1,tempStat, temp(:));

        temp = session(isession).run.totalCorrArray(:,:,iint);
        temp = setUpperTri2NaN(temp);

        tempRun = cat(1,tempRun, temp(:));
    end

    allStat(:,iint) = tempStat;
    allRun(:,iint) = tempRun;
end

allStatTotal = allStat;
allRunTotal = allRun;

subplot(3,2,1), hold on
shadedErrorBar(-100:10:1700, nanmean(allStat,1), nansem(allStat,1))
shadedErrorBar(-100:10:1700, nanmean(allRun,1), nansem(allRun,1), 'lineProps', {'color', 'r'})
ax = gca;
ax.XTick = -200:100:1800;
xlim([-200 1800])
title('total r_{sc} mean')
defaultAxesProperties(gca, false)


subplot(3,2,2), hold on
plot(-100:10:1700, nanstd(allStat,[],1),'k');
plot(-100:10:1700, nanstd(allRun,[],1),'r')
ax=gca;
ax.XTick = -200:100:1800;
xlim([-200 1800])
title('total r_{sc} std')
defaultAxesProperties(gca, false)


%%%%%%% SIGNAL $$$$$$$
figure(123)
allStat = [];
allRun = [];

for iint = 1:nTimeBins
    tempStat = [];
    tempRun = [];
    for isession = sesh2include
        temp = session(isession).stat.signalCorrArray(:,:,iint);
        temp = setUpperTri2NaN(temp);

        tempStat = cat(1,tempStat, temp(:));

        temp = session(isession).run.signalCorrArray(:,:,iint);
        temp = setUpperTri2NaN(temp);

        tempRun = cat(1,tempRun, temp(:));
    end

    allStat(:,iint) = tempStat;
    allRun(:,iint) = tempRun;
end

subplot(3,2,3), hold on
shadedErrorBar(-100:10:1700, nanmean(allStat,1), nansem(allStat,1))
shadedErrorBar(-100:10:1700, nanmean(allRun,1), nansem(allRun,1), 'lineProps', {'color', 'r'})
ax = gca;
ax.XTick = -200:100:1800;
xlim([-200 1800])
title('signal r_{sc} mean')
defaultAxesProperties(gca, false)


subplot(3,2,4), hold on
plot(-100:10:1700, nanstd(allStat,[],1),'k');
plot(-100:10:1700, nanstd(allRun,[],1),'r')
ax = gca;
ax.XTick = -200:100:1800;
xlim([-200 1800])
title('signal r_{sc} std')
defaultAxesProperties(gca, false)

% figure, hold on
% title('signal correlations')
% plot(nanmean(allStat),nanstd(allStat,[],1),'k.')
% plot(nanmean(allRun),nanstd(allRun,[],1),'r.')

allStatSig = allStat;
allRunSig = allRun;


%%%%%%% NOISE $$$$$$$
figure(123)
allStat = [];
allRun = [];


for iint = 1:nTimeBins
    tempStat = [];
    tempRun = [];
    for isession = sesh2include
        temp = session(isession).stat.noiseCorrArray(:,:,iint);
        temp = setUpperTri2NaN(temp);

        tempStat = cat(1,tempStat, temp(:));

        temp = session(isession).run.noiseCorrArray(:,:,iint);
        temp = setUpperTri2NaN(temp);

        tempRun = cat(1,tempRun, temp(:));
    end

    allStat(:,iint) = tempStat;
    allRun(:,iint) = tempRun;
end

subplot(3,2,5), hold on
shadedErrorBar(-100:10:1700, nanmean(allStat,1), nansem(allStat,1))
shadedErrorBar(-100:10:1700, nanmean(allRun,1), nansem(allRun,1), 'lineProps', {'color', 'r'})
ax = gca;
ax.XTick = -200:100:1800;
xlim([-200 1800])
title('noise r_{sc} mean')
defaultAxesProperties(gca, false)


subplot(3,2,6), hold on
plot(-100:10:1700, nanstd(allStat,[],1),'k');
plot(-100:10:1700, nanstd(allRun,[],1),'r')
ax = gca;
ax.XTick = -200:100:1800;
xlim([-200 1800])
title('noise r_{sc} std')
defaultAxesProperties(gca, false)

% figure, hold on
% title('noise correlations')
% plot(nanmean(allStat),nanstd(allStat,[],1),'k')
% plot(nanmean(allRun),nanstd(allRun,[],1),'r')

allStatNoise = allStat;
allRunNoise = allRun;

%% plot abs corrs
figure
subplot(311)
shadedErrorBar(-100:10:1700, nanmean(abs(allStatTotal),1), nansem(abs(allStatTotal),1),'lineProps','k')
shadedErrorBar(-100:10:1700, nanmean(abs(allRunTotal),1), nansem(abs(allRunTotal),1),'lineProps','r')
ax = gca;
ax.XTick = -200:100:1800;
xlim([-200 1800])
title('total')
defaultAxesProperties(gca, false)

subplot(312)
shadedErrorBar(-100:10:1700, nanmean(abs(allStatSig),1), nansem(abs(allStatSig),1),'lineProps','k')
shadedErrorBar(-100:10:1700, nanmean(abs(allRunSig),1), nansem(abs(allRunSig),1),'lineProps','r')
ax = gca;
ax.XTick = -200:100:1800;
xlim([-200 1800])
title('signal')
defaultAxesProperties(gca, false)

subplot(313)
shadedErrorBar(-100:10:1700, nanmean(abs(allStatNoise),1), nansem(abs(allStatNoise),1),'lineProps','k')
shadedErrorBar(-100:10:1700, nanmean(abs(allRunNoise),1), nansem(abs(allRunNoise),1),'lineProps','r')
ax = gca;
ax.XTick = -200:100:1800;
xlim([-200 1800])
title('noise')
defaultAxesProperties(gca, false)



figure

subplot(121)
shadedErrorBar(-100:10:1700, nanmean(abs(allStatSig),1), nansem(abs(allStatSig),1),'lineProps','k')
shadedErrorBar(-100:10:1700, nanmean(abs(allRunSig),1), nansem(abs(allRunSig),1),'lineProps','r')
ax = gca;
ax.XTick = -200:100:1800;
xlim([-200 1800])
% ylim([0, 0.16])
title('signal')
defaultAxesProperties(gca, false)

subplot(122)
shadedErrorBar(-100:10:1700, nanmean(abs(allStatNoise),1), nansem(abs(allStatNoise),1),'lineProps','k')
shadedErrorBar(-100:10:1700, nanmean(abs(allRunNoise),1), nansem(abs(allRunNoise),1),'lineProps','r')
ax = gca;
ax.XTick = -200:100:1800;
xlim([-200 1800])
% ylim([0, 0.16])
title('noise')
defaultAxesProperties(gca, false)




%% stats for mangitude of correlations
%noise

for isession = 1:5

    tempStat = [];
    tempRun = [];
        for iint = 1:181

        temp = session(isession).stat.noiseCorrArray(:,:,iint);
        temp = setUpperTri2NaN(temp);

        tempStat = cat(2,tempStat, temp(:));

        temp = session(isession).run.noiseCorrArray(:,:,iint);
        temp = setUpperTri2NaN(temp);

        tempRun = cat(2,tempRun, temp(:));
        end

        sesh(isession).statNoise = removeNANrows(nanmean(tempStat(:,11:111),2));
        sesh(isession).runNoise = removeNANrows(nanmean(tempRun(:,11:111),2));

end


% stats on magnitude of noise
allStat = abs(cat(1,sesh.statNoise));
allRun = abs(cat(1,sesh.runNoise));

allVals = cat(1,allStat,allRun);
nSesh = cellfun(@numel, {sesh.statNoise});

sessionVec = repelem(1:5,1,nSesh)';

valsVec = allVals;
pairVec = categorical(cat(2,1:numel(allStat), 1:numel(allRun))');
stateVec = categorical(cat(1,repelem(1,numel(allStat),1), repelem(2,numel(allRun),1)));
seshVec = categorical(cat(1,sessionVec, sessionVec));

tbl = table(valsVec,pairVec,stateVec,seshVec,...
    'VariableNames',{'vals','pair','state','sesh'});

    f = 'vals ~ state + (1|pair) + (1|sesh)';
      
    lme = fitlme(tbl, f, 'DummyVarCoding', 'reference')

%% signal

for isession = 1:5

    tempStat = [];
    tempRun = [];
        for iint = 1:181

        temp = session(isession).stat.signalCorrArray(:,:,iint);
        temp = setUpperTri2NaN(temp);

        tempStat = cat(2,tempStat, temp(:));

        temp = session(isession).run.signalCorrArray(:,:,iint);
        temp = setUpperTri2NaN(temp);

        tempRun = cat(2,tempRun, temp(:));
        end

        sesh(isession).statSignal = removeNANrows(nanmean(tempStat(:,11:111),2));
        sesh(isession).runSignal = removeNANrows(nanmean(tempRun(:,11:111),2));

end


% stats on magnitude of noise
allStat = abs(cat(1,sesh.statSignal));
allRun = abs(cat(1,sesh.runSignal));

allVals = cat(1,allStat,allRun);
nSesh = cellfun(@numel, {sesh.statSignal});

sessionVec = repelem(1:5,1,nSesh)';

valsVec = allVals;
pairVec = categorical(cat(2,1:numel(allStat), 1:numel(allRun))');
stateVec = categorical(cat(1,repelem(1,numel(allStat),1), repelem(2,numel(allRun),1)));
seshVec = categorical(cat(1,sessionVec, sessionVec));

tbl = table(valsVec,pairVec,stateVec,seshVec,...
    'VariableNames',{'vals','pair','state','sesh'});

    f = 'vals ~ state + (1|pair) + (1|sesh)';
      
    lme = fitlme(tbl, f, 'DummyVarCoding', 'reference')

%% rergression between signal and noise correlations over time


for isession = 1:5
    statRho = [];
    statCosTheta = [];
    for iint = 1:181
        temp = session(isession).stat.signalCorrArray(:,:,iint);
        temp = setUpperTri2NaN(temp);
        tempStatSignal = temp;

        temp = session(isession).stat.noiseCorrArray(:,:,iint);
        temp = setUpperTri2NaN(temp);
        tempStatNoise = temp;

        allVals = [tempStatSignal(:), tempStatNoise(:)];
        allVals = removeNANrows(allVals);


        [statRho(iint,:), statS(iint)] = polyfit(allVals(:,1),allVals(:,2),1);
        statCosTheta(iint) = calcCosineTheta(allVals(:,1),allVals(:,2));
    end

    tempS(isession).statRho = statRho(:,1);
    tempS(isession).statCosTheta = statCosTheta;
    
end



for isession = 1:5
    runRho = [];
    runCosTheta = [];
    for iint = 1:181
        temp = session(isession).run.signalCorrArray(:,:,iint);
        temp = setUpperTri2NaN(temp);
        tempRunSignal = temp;

        temp = session(isession).run.noiseCorrArray(:,:,iint);
        temp = setUpperTri2NaN(temp);
        tempRunNoise = temp;

        allVals = [tempRunSignal(:), tempRunNoise(:)];
        allVals = removeNANrows(allVals);


        [runRho(iint,:), runS(iint)] = polyfit(allVals(:,1),allVals(:,2),1);
        runCosTheta(iint) = calcCosineTheta(allVals(:,1),allVals(:,2));

    end

    tempS(isession).runRho = runRho(:,1);
    tempS(isession).runCosTheta = runCosTheta;
end

idx2use = 15:120;
bv = -100:10:1700;

allStat = cat(2,tempS.statRho)';
allRun = cat(2,tempS.runRho)';


figure, hold on
shadedErrorBar(bv(idx2use), mean(allStat(:,idx2use),1), sem(allStat(:,idx2use),1),'lineProps','k')
shadedErrorBar(bv(idx2use), mean(allRun(:,idx2use),1), sem(allRun(:,idx2use),1),'lineProps','r')
xlim([-200 1800])
ax = gca; ax.XTick = -200:200:1800;
plot([-200 1800], [0 0], 'k:')
ylim([-.2 1.3])
ax.XTick = 0:200:1800;
xlim([0 1200])
defaultAxesProperties(gca,true)


% rm-anova on slope

statVals = allStat(:,idx2use)';
runVals  = allRun(:,idx2use)';


slopeVec = cat(1,statVals(:), runVals(:));
stateVec = categorical(cat(1,repelem(1, numel(statVals(:)))', repelem(2,numel(runVals(:)))'));
timeVec = categorical(cat(1,repmat(1:size(statVals,1), 1,size(statVals,2))',repmat(1:size(statVals,1), 1,size(runVals,2))'));
subjVec = categorical(cat(1,repelem(1:size(statVals,2),1, size(statVals,1))', repelem(1:size(runVals,2),1, size(statVals,1))'));


[p,tbl,stats,terms] = anovan(slopeVec,{timeVec,stateVec,subjVec},'model',2,'random',3,'varnames',{'Time','State','Subj'});

%%
%%% cell-type %%%



%% Cell-type total, signal and noise correlations

nTimeBins = numel(1:200-19);
sesh2include=1:5;

clear tempStat tempRun stat run

%%%%%%% TOTAL $$$$$$$

for itype1=1:3
    for itype2=1:3
        stat(itype1,itype2).corr = [];
        run(itype1,itype2).corr = [];
    end
end

for iint = 1:nTimeBins

    for itype1=1:3
        for itype2=1:3
            tempStat(itype1,itype2).corr = [];
            tempRun(itype1,itype2).corr = [];
        end
    end

    for isession = sesh2include

        for itype1 = 1:3
            for itype2 = 1:3

                t1idx = [session(isession).units.cellType]==itype1;
                t2idx = [session(isession).units.cellType]==itype2;

                temp = session(isession).stat.totalCorrArray(:,:,iint);
%                 temp = setUpperTri2NaN(temp);

                thisTemp = temp(t1idx,t2idx);

                tempStat(itype1,itype2).corr = cat(1,tempStat(itype1,itype2).corr, thisTemp(:));

                temp = session(isession).run.totalCorrArray(:,:,iint);
%                 temp = setUpperTri2NaN(temp);

                thisTemp = temp(t1idx,t2idx);

                tempRun(itype1,itype2).corr = cat(1,tempRun(itype1,itype2).corr, thisTemp(:));


            end
        end
    end

        for itype1 = 1:3
            for itype2 = 1:3
                stat(itype1,itype2).corr(:,iint) = tempStat(itype1,itype2).corr(:);
                run(itype1,itype2).corr(:,iint) = tempRun(itype1,itype2).corr(:);
            end
        end
  
end

totalStatType = stat;
totalRunType = run;

clear tempStat tempRun stat run

%%%%%%% SIGNAL $$$$$$$

for itype1=1:3
    for itype2=1:3
        stat(itype1,itype2).corr = [];
        run(itype1,itype2).corr = [];
    end
end

for iint = 1:nTimeBins

    for itype1=1:3
        for itype2=1:3
            tempStat(itype1,itype2).corr = [];
            tempRun(itype1,itype2).corr = [];
        end
    end

    for isession = sesh2include

        for itype1 = 1:3
            for itype2 = 1:3

                t1idx = [session(isession).units.cellType]==itype1;
                t2idx = [session(isession).units.cellType]==itype2;

                temp = session(isession).stat.signalCorrArray(:,:,iint);
%                 temp = setUpperTri2NaN(temp);

                thisTemp = temp(t1idx,t2idx);

                tempStat(itype1,itype2).corr = cat(1,tempStat(itype1,itype2).corr, thisTemp(:));

                temp = session(isession).run.signalCorrArray(:,:,iint);
%                 temp = setUpperTri2NaN(temp);

                thisTemp = temp(t1idx,t2idx);

                tempRun(itype1,itype2).corr = cat(1,tempRun(itype1,itype2).corr, thisTemp(:));


            end
        end
    end

        for itype1 = 1:3
            for itype2 = 1:3
                stat(itype1,itype2).corr(:,iint) = tempStat(itype1,itype2).corr(:);
                run(itype1,itype2).corr(:,iint) = tempRun(itype1,itype2).corr(:);
            end
        end
  
end

signalStatType = stat;
signalRunType = run;


clear tempStat tempRun stat run

%%%%%%% NOISE $$$$$$$

for itype1=1:3
    for itype2=1:3
        stat(itype1,itype2).corr = [];
        run(itype1,itype2).corr = [];
    end
end

for iint = 1:nTimeBins

    for itype1=1:3
        for itype2=1:3
            tempStat(itype1,itype2).corr = [];
            tempRun(itype1,itype2).corr = [];
        end
    end

    for isession = sesh2include

        for itype1 = 1:3
            for itype2 = 1:3

                t1idx = [session(isession).units.cellType]==itype1;
                t2idx = [session(isession).units.cellType]==itype2;

                temp = session(isession).stat.noiseCorrArray(:,:,iint);
%                 temp = setUpperTri2NaN(temp);

                thisTemp = temp(t1idx,t2idx);

                tempStat(itype1,itype2).corr = cat(1,tempStat(itype1,itype2).corr, thisTemp(:));

                temp = session(isession).run.noiseCorrArray(:,:,iint);
%                 temp = setUpperTri2NaN(temp);

                thisTemp = temp(t1idx,t2idx);

                tempRun(itype1,itype2).corr = cat(1,tempRun(itype1,itype2).corr, thisTemp(:));


            end
        end
    end

        for itype1 = 1:3
            for itype2 = 1:3
                stat(itype1,itype2).corr(:,iint) = tempStat(itype1,itype2).corr(:);
                run(itype1,itype2).corr(:,iint) = tempRun(itype1,itype2).corr(:);
            end
        end
  
end

noiseStatType = stat;
noiseRunType = run;

%% plot means
figure
for itype1=1:3
    for itype2=1:3
    subplot(3,3,itype1*3-3+itype2)
    hold on
    shadedErrorBar(-100:10:1700, mean(totalStatType(itype1,itype2).corr,1,'omitnan'),nansem(totalStatType(itype1,itype2).corr,1),'lineProps','k')
    shadedErrorBar(-100:10:1700, mean(totalRunType(itype1,itype2).corr,1,'omitnan'),nansem(totalRunType(itype1,itype2).corr,1),'lineProps','r')
    ylim([0 0.22])
    title(num2str([itype1 itype2]))
    xlim([-200 1800]);
    ax = gca; ax.XTick = -200:200:1800;  ax.XTickLabel=[];
    defaultAxesProperties(gca,false)
    xlim()
    end
end

figure
for itype1=1:3
    for itype2=1:3
    subplot(3,3,itype1*3-3+itype2)
    hold on
    shadedErrorBar(-100:10:1700, mean(signalStatType(itype1,itype2).corr,1,'omitnan'),nansem(signalStatType(itype1,itype2).corr,1),'lineProps','k')
    shadedErrorBar(-100:10:1700, mean(signalRunType(itype1,itype2).corr,1,'omitnan'),nansem(signalRunType(itype1,itype2).corr,1),'lineProps','r')
    ylim([-0.01 0.06])
    xlim([-200 1800]);
     ax = gca; ax.XTick = -200:200:1800; ax.XTickLabel=[];
    defaultAxesProperties(gca,false)
    title(num2str([itype1 itype2]))
    end
end

figure
for itype1=1:3
    for itype2=1:3
    subplot(3,3,itype1*3-3+itype2)
    hold on
    shadedErrorBar(-100:10:1700, mean(noiseStatType(itype1,itype2).corr,1,'omitnan'),nansem(noiseStatType(itype1,itype2).corr,1),'lineProps','k')
    shadedErrorBar(-100:10:1700, mean(noiseRunType(itype1,itype2).corr,1,'omitnan'),nansem(noiseRunType(itype1,itype2).corr,1),'lineProps','r')
    ylim([0 0.2])
    xlim([-200 1800]);
     ax = gca; ax.XTick = -200:200:1800;  ax.XTickLabel=[];
    defaultAxesProperties(gca,false)
    title(num2str([itype1 itype2]))
    end
end


%% plot std
figure
for itype1=1:3
    for itype2=1:3
    subplot(3,3,itype1*3-3+itype2)
    hold on
    plot(-100:10:1700, nanstd(totalStatType(itype1,itype2).corr),'k')
    plot(-100:10:1700, nanstd(totalRunType(itype1,itype2).corr),'r')
      ylim([0.12 0.27])
    title(num2str([itype1 itype2]))
    xlim([-200 1800])
    ax = gca; ax.XTick = -200:200:1800;  ax.XTickLabel=[];
    defaultAxesProperties(gca,false)
    xlim()
    end
end

figure
for itype1=1:3
    for itype2=1:3
    subplot(3,3,itype1*3-3+itype2)
    hold on
    plot(-100:10:1700, nanstd(signalStatType(itype1,itype2).corr),'k')
    plot(-100:10:1700, nanstd(signalRunType(itype1,itype2).corr),'r')
     ylim([0 0.2])
    xlim([-200 1800]);
     ax = gca; ax.XTick = -200:200:1800; ax.XTickLabel=[];
    defaultAxesProperties(gca,false)
    title(num2str([itype1 itype2]))
    end
end

figure
for itype1=1:3
    for itype2=1:3
    subplot(3,3,itype1*3-3+itype2)
    hold on
    plot(-100:10:1700, nanstd(noiseStatType(itype1,itype2).corr),'k')
    plot(-100:10:1700, nanstd(noiseRunType(itype1,itype2).corr),'r')
    ylim([0.09 0.25])
    xlim([-200 1800])
     ax = gca; ax.XTick = -200:200:1800;  ax.XTickLabel=[];
    defaultAxesProperties(gca,false)
    title(num2str([itype1 itype2]))
    end
end



%% plot absmean
figure
for itype1=1:3
    for itype2=1:3
    subplot(3,3,itype1*3-3+itype2)
    hold on
    shadedErrorBar(-100:10:1700, mean(abs(totalStatType(itype1,itype2).corr),1,'omitnan'),nansem(abs(totalStatType(itype1,itype2).corr),1),'lineProps','k')
    shadedErrorBar(-100:10:1700, mean(abs(totalRunType(itype1,itype2).corr),1,'omitnan'),nansem(abs(totalRunType(itype1,itype2).corr),1),'lineProps','r')
    ylim([0.1 0.25])
    title(num2str([itype1 itype2]))
    xlim([-200 1800])
    ax = gca; ax.XTick = -200:200:1800;  ax.XTickLabel=[];
    defaultAxesProperties(gca,false)
    xlim()
    end
end

figure
for itype1=1:3
    for itype2=1:3
    subplot(3,3,itype1*3-3+itype2)
    hold on
    shadedErrorBar(-100:10:1700, mean(abs(signalStatType(itype1,itype2).corr),1,'omitnan'),nansem(abs(signalStatType(itype1,itype2).corr),1),'lineProps','k')
    shadedErrorBar(-100:10:1700, mean(abs(signalRunType(itype1,itype2).corr),1,'omitnan'),nansem(abs(signalRunType(itype1,itype2).corr),1),'lineProps','r')
     ylim([0 0.15])
    xlim([-200 1800]);
     ax = gca; ax.XTick = -200:200:1800; ax.XTickLabel=[];
    defaultAxesProperties(gca,false)
    title(num2str([itype1 itype2]))
    end
end

figure
for itype1=1:3
    for itype2=1:3
    subplot(3,3,itype1*3-3+itype2)
    hold on
    shadedErrorBar(-100:10:1700, mean(abs(noiseStatType(itype1,itype2).corr),1,'omitnan'),nansem(abs(noiseStatType(itype1,itype2).corr),1),'lineProps','k')
    shadedErrorBar(-100:10:1700, mean(abs(noiseRunType(itype1,itype2).corr),1,'omitnan'),nansem(abs(noiseRunType(itype1,itype2).corr),1),'lineProps','r')
     ylim([0.05 0.2])
    xlim([-200 1800])
     ax = gca; ax.XTick = -200:200:1800;  ax.XTickLabel=[];
    defaultAxesProperties(gca,false)
    title(num2str([itype1 itype2]))
    end
end

%% regression between signal and noise

for itype1 = 1:3
    for itype2 = 1:3

    statrho=[]; runrho=[];
    for iint = 1:181
        statSigTemp = signalStatType(itype1,itype2).corr(:,iint); 
        statNoiseTemp = noiseStatType(itype1,itype2).corr(:,iint);
        statSigTemp(isnan(statSigTemp)) = []; statNoiseTemp(isnan(statNoiseTemp))=[];
        [statrho(iint,:), statS(iint)] = polyfit(statSigTemp, statNoiseTemp,1);

        runSigTemp = signalRunType(itype1,itype2).corr(:,iint); 
        runNoiseTemp = noiseRunType(itype1,itype2).corr(:,iint); 
        runSigTemp(isnan(runSigTemp)) = []; runNoiseTemp(isnan(runNoiseTemp))=[];
        [runrho(iint,:), runS(iint)] = polyfit(runSigTemp, runNoiseTemp,1);
    end
    

    statrhoSlope = statrho(:,1)';
    statrhoIntercept = statrho(:,2)';

    runrhoSlope = runrho(:,1)';
    runrhoIntercept = runrho(:,2)';

    subplot(3,3,itype1*3-3+itype2), hold on
    plot([-200 1800], [0 0],'k')
    plot(-100:10:1700, statrhoSlope,'k')
    plot(-100:10:1700, runrhoSlope,'r')
    xlim([-200 1800])
    ylim([-0.5 1])
    end
end


%% distributions of noise corrs at different timepoints
for itime = [1, 34, 96, 136]

    statDist = [];
    runDist = [];

    for isession = 1:5
        statDist = vertcat(statDist, nanmean(session(isession).stat.interval(itime).array,2));
        runDist = vertcat(runDist, nanmean(session(isession).run.interval(itime).array,2));
    end
    figure
    subplot(211)
    hold off
    histogram(statDist, 'normalization', 'probability', 'binedges', -0.5:0.05:0.5, 'FaceColor','k')
    hold on
    histogram(runDist, 'normalization', 'probability', 'binedges', -0.5:0.05:0.5, 'FaceColor','r')
    plot(nanmean(statDist), 0.38, 'kv')
    plot(nanmean(runDist), 0.38, 'rv')
    ylim([0, 0.4])
    title(itime)
    defaultAxesProperties(gca, true)
    subplot(212)
    scatter_kde(statDist, runDist, 'filled', 'MarkerSize', 5), colormap(viridis)
    xlim([-1 1])
    ylim([-1 1])
    hold on
    plot([-1 1], [-1 1], 'k')
    axis equal
    xlim([-0.5 0.5])
    ylim([-0.5 0.5])
    defaultAxesProperties(gca, false)
    print(gcf, ['noiseCorrTime_', num2str(itime)], '-vector', '-dsvg')

    % plot([nanmean(statDist), nanmean(statDist)],...
    %     [nanmean(runDist)-nansem(runDist)*1.96, nanmean(runDist)+nansem(runDist)*1.96], 'c')
    % plot([nanmean(statDist)-nansem(statDist)*1.96, nanmean(statDist)+nansem(statDist)*1.96],...
    %     [nanmean(runDist), nanmean(runDist)], 'c')

end


%% example





%% plot example units
% 30 and 21 from session 1
%figure
%ax=gca;
figure
cols=inferno(6);
for iunit = 43:numel(session(1).units)

    data = smoothdata(cell2mat(cellfun(@(x) mean(x,2), session(1).units(iunit).allSpikes(:,1), 'UniformOutput', false)').*100,1,'gaussian', 17.5);
    for icond = 1:6
        plot(data(:,icond),'Color', cols(icond,:))
        hold on

    end
    title(iunit)
    pause
    hold off
end
%% plot 2 units PSTHs +  correlations
% session 1: [21 30], [30 31], [18 21] [36 42] [44 21]
% good HIGH signal: [35 40], [35 102]
% noise spike early: [68 123], [44 158]
% big signal and noise: [81 173]

% session 2: nise spike at dip [88 50], [5 64], [84 52]x, [79, 76]
isession = 2;

%
 iunit1 = round(rand*numel(session(isession).units));
 iunit2 = round(rand*numel(session(isession).units));
% 
%  session2:
%    iunit1 = 84; 53, 83, 35, 53*
%   iunit2 = 52; 50, 26, 52, 57*
iunit1 = 53; iunit2=57;
cols=inferno(6);




figure

for istate = 1:2

subplot(5,2,0+istate), hold off
data = smoothdata(cell2mat(cellfun(@(x) mean(x,2), session(isession).units(iunit1).allSpikes(:,istate), 'UniformOutput', false)').*100,1,'gaussian', 17.5);
for icond = 1:6
    plot(-200:10:1790, data(:,icond),'Color', cols(icond,:))
    hold on
end
title(iunit1)
xlim([-200, 1800])
ax= gca; ax.XTick = -200:200:1800; % ax.YLim = [0 66]; ax.YTick = ax.YLim;
defaultAxesProperties(gca, true)

subplot(5,2,2+istate), hold off
data = smoothdata(cell2mat(cellfun(@(x) mean(x,2), session(isession).units(iunit2).allSpikes(:,istate), 'UniformOutput', false)').*100,1,'gaussian', 17.5);
for icond = 1:6
    plot(-200:10:1790, data(:,icond),'Color', cols(icond,:))
    hold on
end
title(iunit2)
xlim([-200, 1800])
ax= gca; ax.XTick = -200:200:1800; %ax.YLim = [0 55]; ax.YTick = ax.YLim;
defaultAxesProperties(gca, true)

if istate ==1
    fname = 'stat';
elseif istate==2
    fname = 'run';
end
subplot(5,2,4+istate)
plot(-100:10:1700, squeeze(session(isession).(fname).totalCorrArray(iunit1,iunit2,:)),'k')
title('total r_{sc}')
xlim([-200, 1800])
ax= gca; ax.XTick = -200:200:1800; ax.YLim = [-0.25 0.8]; 
defaultAxesProperties(gca, true)


subplot(5,2,6+istate)
plot(-100:10:1700, squeeze(session(isession).(fname).signalCorrArray(iunit1,iunit2,:)),'k')
title('signal r_{sc}')
xlim([-200, 1800])
ax= gca; ax.XTick = -200:200:1800; ax.YLim = [-0.25 0.4]; 
defaultAxesProperties(gca, true)

subplot(5,2,8+istate)
plot(-100:10:1700, squeeze(session(isession).(fname).noiseCorrArray(iunit1,iunit2,:)),'k')
title('noise r_{sc}')
xlim([-200, 1800])
ax= gca; ax.XTick = -200:200:1800; ax.YLim = [-0.25 0.8];
defaultAxesProperties(gca, true)

% pause
end
%%

figure

for istate = 1:2
if istate ==1
    fname = 'stat';
elseif istate==2
    fname = 'run';
end
subplot(3,2,0+istate)
plot(-100:10:1700, squeeze(session(isession).(fname).totalCorrArray(iunit1,iunit2,:)),'k')
title('total r_{sc}')
xlim([-200, 1800])
ax= gca; ax.XTick = -200:200:1800; ax.YLim = [-0.2 0.8]; ax.YTick=[-0.2:0.2:0.8];
defaultAxesProperties(gca, true)


subplot(3,2,2+istate)
plot(-100:10:1700, squeeze(session(isession).(fname).signalCorrArray(iunit1,iunit2,:)),'k')
title('signal r_{sc}')
xlim([-200, 1800])
ax= gca; ax.XTick = -200:200:1800; ax.YLim = [-0.2 0.4]; ax.YTick=[-0.2:0.2:0.4];
defaultAxesProperties(gca, true)

subplot(3,2,4+istate)
plot(-100:10:1700, squeeze(session(isession).(fname).noiseCorrArray(iunit1,iunit2,:)),'k')
title('noise r_{sc}')
xlim([-200, 1800])
ax= gca; ax.XTick = -200:200:1800; ax.YLim = [-0.2 0.8]; ax.YTick=[-0.2:0.2:0.8];
defaultAxesProperties(gca, true)
end


%% sort and plot correlation matrices at different time points
figure
isession = 2;
bv = -100:10:1700;
time_idx = [21 31 41 91];

% choose which matrix to sort with
dm = session(isession).run.signalCorrArray(:,:,time_idx(4));
dm(isnan(dm)) = 0;
% idx2remove = find(all(isnan(dm),2));
% dm(idx2remove,:) = [];
% dm(:,idx2remove) = [];

full_dm = 1-dm;
min_dm=min(full_dm(:));
dm2 = full_dm-min_dm;
dm3 = dm2/max(dm2(:));


Z_sub = linkage(dm3, 'average');
% leafOrder = optimalleaforder(Z_sub,dm3');
[denHandle, leafnode_idx, outpermNodeOrder] = dendrogram(Z_sub, size(dm3,1), 'orientation', 'right'); 
unique(cluster(Z_sub,"cutoff",1.15))
newOrder = [];


for ileaf = 1:numel(outpermNodeOrder)
    idx = find(leafnode_idx==outpermNodeOrder(ileaf));
    newOrder = vertcat(newOrder, idx(:));
end


for itime = 1:numel(time_idx)

    dm = session(isession).run.signalCorrArray(:,:,time_idx(itime));
    dm(isnan(dm)) = 0;

%     dm(idx2remove,:) = [];
%     dm(:,idx2remove) = [];

dm_reordered=[];
for ipsth1 = 1:numel(newOrder)
    for ipsth2 = 1:numel(newOrder)
        dm_reordered(ipsth1,ipsth2) = dm(newOrder(ipsth1),newOrder(ipsth2));
    end
end

M=dm_reordered;


M = M.*-(eye(height(M))-1);
dm_reordered = M;


subplot(1,4,itime)

imagesc(dm_reordered)
title(bv(time_idx(itime)))
caxis([-1 1])
colormap(redblue)
axis xy
ax = gca;
box off
ax = gca; ax.XTick=[]; ax.YTick=[];

defaultAxesProperties(gca, false)
end
%%

time_idx = [1 24 96 136];
for itime = 1:numel(time_idx)

    dm = session(isession).run.totalCorrArray(:,:,time_idx(itime));
    dm(isnan(dm)) = 0;

%     dm(idx2remove,:) = [];
%     dm(:,idx2remove) = [];

dm_reordered=[];
for ipsth1 = 1:numel(newOrder)
    for ipsth2 = 1:numel(newOrder)
        dm_reordered(ipsth1,ipsth2) = dm(newOrder(ipsth1),newOrder(ipsth2));
    end
end

    subplot(2,4,itime+4)
imagesc(dm_reordered)
caxis([-1 1])
colormap(redblue)
axis xy
ax = gca; ax.XTick=[]; ax.YTick=[];

end

%% example timewise corr plots

time_idx = [21, 41];
isession = 2;

% stat
dm1 = session(isession).stat.signalCorrArray(:,:,time_idx(1));
dm1 = setUpperTri2NaN(dm1);

dm2 = session(isession).stat.signalCorrArray(:,:,time_idx(2));
dm2 = setUpperTri2NaN(dm2);
X = removeNANrows([dm1(:), dm2(:)]);
rstat = corr(X(:,1), X(:,2));

figure
subplot(121)
plot(X(:,1), X(:,2), 'k.')
title(rstat)
axis equal
xlim([-1 1]), ylim([-1 1])

% run
dm1 = session(isession).run.signalCorrArray(:,:,time_idx(1));
dm1 = setUpperTri2NaN(dm1);

dm2 = session(isession).run.signalCorrArray(:,:,time_idx(2));
dm2 = setUpperTri2NaN(dm2);
X = removeNANrows([dm1(:), dm2(:)]);
rrun = corr(X(:,1), X(:,2));

subplot(122)
plot(X(:,1), X(:,2), 'r.')
title(rrun)
axis equal
xlim([-1 1]), ylim([-1 1])



%% reorder
% choose which matrix to sort with
dm = session(1).run.signalCorrArray(:,:,time_idx(end));
dm(isnan(dm)) = 0;

full_dm = 1-dm;
min_dm=min(full_dm(:));
dm2 = full_dm-min_dm;
dm3 = dm2/max(dm2(:));


Z_sub = linkage(dm3, 'average');
% leafOrder = optimalleaforder(Z_sub,dm3');
[denHandle, leafnode_idx, outpermNodeOrder] = dendrogram(Z_sub, size(dm3,1), 'orientation', 'right'); 

newOrder = [];


for ileaf = 1:numel(outpermNodeOrder)
    idx = find(leafnode_idx==outpermNodeOrder(ileaf));
    newOrder = vertcat(newOrder, idx(:));
end


%% 3d matrix of correlations


% consider sorting by t=1

  %newOrder = randperm(140);

M = shiftdim(session(1).run.signalCorrArray,3);
% M(140,:,:) = [];
% M(:,140,:) = [];
M(isnan(M))=0;


dm_reordered=zeros(size(M));

for itime = 1:size(M,3)

for ipsth1 = 1:numel(newOrder)
    for ipsth2 = 1:numel(newOrder)
        dm_reordered(ipsth1,ipsth2,itime) = M(newOrder(ipsth1),newOrder(ipsth2),itime);
    end
end
end

M=dm_reordered;

M1 = M(:,:,1);
M1 = M1.*-(eye(height(M1))-1);
M(:,:,1)=M1;


xslice = size(M,1);                             % define the cross sections to view
yslice = size(M,2);
zslice = time_idx;

x = 1:size(M,1); y = 1:size(M,2); z = 1:size(M,3);
figure
h = slice(x, y, z, M, xslice, yslice, zslice);    % display the slices
  set(h,'EdgeColor','none');
% set(h,'EdgeAlpha',0.1)
% ylim([-3 3])
view(-34,24)

% cb = colorbar;                                 
colormap(redblue)
caxis([-.5 .5])
ax=gca; ax.CameraUpVector = [0 1 0];
xlabel('x')
ylabel('y')
zlabel('z')

set(gca, 'DataAspectRatio', [diff(get(gca, 'XLim')) diff(get(gca, 'XLim')) diff(get(gca, 'ZLim')/2)])
defaultAxesProperties(gca, false)
grid off
ax.ZTick = 1:20:181;
% 
direction = [0 1 1];
rotate(h,direction,5);

direction = [1 1 0];
rotate(h,direction,5)
rotate(h,direction,5)



















%% slope calculation v2

for isession = 1:5
    statRho = [];
    statCosTheta = [];
    for iint = 1:181
        temp = session(isession).stat.signalCorrArray(:,:,iint);
        temp = setUpperTri2NaN(temp);
        tempStatSignal = temp;

        temp = session(isession).stat.noiseCorrArray(:,:,iint);
        temp = setUpperTri2NaN(temp);
        tempStatNoise = temp;

        allVals = [tempStatSignal(:), tempStatNoise(:)];
        allVals = removeNANrows(allVals);


        [statRho(iint,:), statS(iint)] = polyfit(allVals(:,1),allVals(:,2),1);
        statCosTheta(iint) = calcCosineTheta(allVals(:,1),allVals(:,2));
    end

    tempS(isession).statRho = statRho(:,1);
    tempS(isession).statCosTheta = statCosTheta;
    
end



for isession = 1:5
    runRho = [];
    runCosTheta = [];
    for iint = 1:181
        temp = session(isession).run.signalCorrArray(:,:,iint);
        temp = setUpperTri2NaN(temp);
        tempRunSignal = temp;

        temp = session(isession).run.noiseCorrArray(:,:,iint);
        temp = setUpperTri2NaN(temp);
        tempRunNoise = temp;

        allVals = [tempRunSignal(:), tempRunNoise(:)];
        allVals = removeNANrows(allVals);


        [runRho(iint,:), runS(iint)] = polyfit(allVals(:,1),allVals(:,2),1);
        runCosTheta(iint) = calcCosineTheta(allVals(:,1),allVals(:,2));

    end

    tempS(isession).runRho = runRho(:,1);
    tempS(isession).runCosTheta = runCosTheta;
end

%%

idx2use = 15:120;
bv = -100:10:1700;

allStat = cat(2,tempS.statRho)';
allRun = cat(2,tempS.runRho)';


figure, hold on
shadedErrorBar(bv(idx2use), mean(allStat(:,idx2use),1), sem(allStat(:,idx2use),1),'lineProps','k')
shadedErrorBar(bv(idx2use), mean(allRun(:,idx2use),1), sem(allRun(:,idx2use),1),'lineProps','r')
xlim([-200 1800])
ax = gca; ax.XTick = -200:200:1800;
plot([-200 1800], [0 0], 'k:')
ylim([-.2 1.3])
ax.XTick = 0:200:1800;
xlim([0 1200])
defaultAxesProperties(gca,true)


%% rm-anova on slope

statVals = allStat(:,idx2use)';
runVals  = allRun(:,idx2use)';


slopeVec = cat(1,statVals(:), runVals(:));
stateVec = categorical(cat(1,repelem(1, numel(statVals(:)))', repelem(2,numel(runVals(:)))'));
timeVec = categorical(cat(1,repmat(1:size(statVals,1), 1,size(statVals,2))',repmat(1:size(statVals,1), 1,size(runVals,2))'));
subjVec = categorical(cat(1,repelem(1:size(statVals,2),1, size(statVals,1))', repelem(1:size(runVals,2),1, size(statVals,1))'));


[p,tbl,stats,terms] = anovan(slopeVec,{timeVec,stateVec,subjVec},'model',2,'random',3,'varnames',{'Time','State','Subj'});






%% angle between signal and noise

idx2use = 15:120;
bv = -100:10:1700;

allStat = acosd(cat(1,tempS.statCosTheta));
allRun = acosd(cat(1,tempS.runCosTheta));

figure, hold on
shadedErrorBar(bv(idx2use), mean(allStat(:,idx2use),1), sem(allStat(:,idx2use),1),'lineProps','k')
shadedErrorBar(bv(idx2use), mean(allRun(:,idx2use),1), sem(allRun(:,idx2use),1),'lineProps','r')
xlim([-200 1800])
ax = gca; ax.XTick = -200:200:1800;
plot([-200 1800], [90 90], 'k:')
% ylim([-.2 1.3])
ax.XTick = 0:200:1800;
xlim([0 1200])
defaultAxesProperties(gca,true)



