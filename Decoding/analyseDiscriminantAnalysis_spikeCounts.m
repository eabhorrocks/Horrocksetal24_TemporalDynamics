%% analyse various discriminant analyses

sessionTags = {'M22027', '20220517';...
    'M22029', '20220607';...
    'M22031', '20220524';
    'M22032', '20220621';
    'M22033', '20220706'};

dataDir = 'C:\Users\edward.horrocks\Documents\Code\V1Dynamics\Data\basic_111022';
% dataDir = 'E:\V1Data\Data\basic_111022';


for isession = 1:size(sessionTags,1)
    isession
    fname = [sessionTags{isession,1},'_', sessionTags{isession,2},'_', 'LDAsingle_normal.mat']; %_cumulativeDiscrimAnalysis

    load(fullfile(dataDir,fname))


    session1(isession).stat = session.stat;
    session1(isession).run = session.run;
    session1(isession).units = session.goodUnits;

    % correct incorrectly named struct field
    session1(isession).stat = renameStructField(session1(isession).stat,'shuffleGammaOne', 'shuffleGammaZero');
    session1(isession).run = renameStructField(session1(isession).run,'shuffleGammaOne', 'shuffleGammaZero');


end


session = session1;


%% plot gamma opt decoding, stat vs run, w/ difference plot

statNormOpt=[];
runNormOpt=[];

for isession = 1:5
    statNormOpt=cat(1,statNormOpt,mean(cat(1,session(isession).stat.GammaOpt.perm.meanPerf),1));
    runNormOpt=cat(1,runNormOpt,mean(cat(1,session(isession).run.GammaOpt.perm.meanPerf),1));
end


figure, hold on
shadedErrorBar(-150:10:1750, mean(statNormOpt,1),sem(statNormOpt,1), 'lineProps', 'k')
shadedErrorBar(-150:10:1750, mean(runNormOpt,1),sem(runNormOpt,1), 'lineProps', 'r')
title('Gamma optimised')


diffArray = runNormOpt-statNormOpt;
figure, hold on
shadedErrorBar(-150:10:1750, mean(diffArray,1),sem(diffArray,1), 'lineProps', 'm')
plot([-200 1800], [0 0], 'k')

%% plot individual sessions, normal spike counts, different gamma, LDA

for isession = 1:5

figure, 
tiledlayout('flow'),
nexttile
hold on
plot(mean(cat(1,session(isession).stat.GammaOpt.perm.meanPerf),1));
plot(mean(cat(1,session(isession).stat.GammaOne.perm.meanPerf),1));
plot(mean(cat(1,session(isession).stat.GammaZero.perm.meanPerf),1));
legend({'opt','one','zero'})
ylim([0 1])
nexttile, hold on
plot(mean(cat(1,session(isession).run.GammaOpt.perm.meanPerf),1));
plot(mean(cat(1,session(isession).run.GammaOne.perm.meanPerf),1));
plot(mean(cat(1,session(isession).run.GammaZero.perm.meanPerf),1));
legend({'opt','one','zero'})
ylim([0 1])
end

%% plot averages over sessions, normal counts, average gamma, LDA

allStatGammaOpt=[];
allStatGammaOne=[];
allStatGammaZero=[];
allRunGammaOpt=[];
allRunGammaOne=[];
allRunGammaZero=[];

for isession = 1:5

    allStatGammaOpt=cat(1,allStatGammaOpt,mean(cat(1,session(isession).stat.GammaOpt.perm.meanPerf),1));
    allStatGammaOne=cat(1,allStatGammaOne,mean(cat(1,session(isession).stat.GammaOne.perm.meanPerf),1));
    allStatGammaZero=cat(1,allStatGammaZero,mean(cat(1,session(isession).stat.GammaZero.perm.meanPerf),1));

    allRunGammaOpt=cat(1,allRunGammaOpt,mean(cat(1,session(isession).run.GammaOpt.perm.meanPerf),1));
    allRunGammaOne=cat(1,allRunGammaOne,mean(cat(1,session(isession).run.GammaOne.perm.meanPerf),1));
    allRunGammaZero=cat(1,allRunGammaZero,mean(cat(1,session(isession).run.GammaZero.perm.meanPerf),1));
end

figure, tiledlayout('flow')
nexttile, hold on
plot(mean(allStatGammaOpt,1), 'k')
plot(mean(allStatGammaOne,1), 'k--')
plot(mean(allStatGammaZero,1),'k:')
ylim([0 1])
nexttile, hold on
plot(mean(allRunGammaOpt,1),'r')
plot(mean(allRunGammaOne,1),'r--')
plot(mean(allRunGammaZero,1),'r:')
ylim([0 1])


%% compare LDA and QDA

for isession = 1:5

figure, 
tiledlayout('flow'),
nexttile
hold on
plot(mean(cat(1,session(isession).stat.GammaOpt.perm.meanPerf),1));
plot(mean(cat(1,session(isession).stat.quadratic.perm.meanPerf),1));
plot(mean(cat(1,session(isession).stat.diagquadratic.perm.meanPerf),1));
legend({'LDA','QDA','QDAdiag'})
ylim([0 1])

nexttile, hold on
plot(mean(cat(1,session(isession).run.GammaOpt.perm.meanPerf),1));
plot(mean(cat(1,session(isession).run.quadratic.perm.meanPerf),1));
plot(mean(cat(1,session(isession).run.diagquadratic.perm.meanPerf),1));
legend({'LDA','QDA','QDAdiag'})
ylim([0 1])
end



%% compare normal vs shuffled spike counts, LDA, gamma optimised


for isession = 1:5
figure, 
tiledlayout('flow'),
nexttile
hold on
plot(mean(cat(1,session(isession).stat.GammaOpt.perm.meanPerf),1),'k');
plot(mean(cat(1,session(isession).stat.shuffleGammaOpt.perm.meanPerf),1),'k:');
ylim([0 1])

nexttile, hold on
plot(mean(cat(1,session(isession).run.GammaOpt.perm.meanPerf),1),'r');
plot(mean(cat(1,session(isession).run.shuffleGammaOpt.perm.meanPerf),1),'r:');
ylim([0 1])
end

%% mean normal vs shuffled spike counts, gamma optimised

statNormOpt=[];
statShufOpt=[];
runNormOpt=[];
runShufOpt=[];

for isession = 1:5
    statNormOpt=cat(1,statNormOpt,mean(cat(1,session(isession).stat.GammaOpt.perm.meanPerf),1));
    statShufOpt=cat(1,statShufOpt,mean(cat(1,session(isession).stat.shuffleGammaOpt.perm.meanPerf),1));
    runNormOpt=cat(1,runNormOpt,mean(cat(1,session(isession).run.GammaOpt.perm.meanPerf),1));
    runShufOpt=cat(1,runShufOpt,mean(cat(1,session(isession).run.shuffleGammaOpt.perm.meanPerf),1));
end


figure, hold on
shadedErrorBar(-100:10:1700, mean(statNormOpt,1),sem(statNormOpt,1), 'lineProps', 'k')
shadedErrorBar(-100:10:1700, mean(statShufOpt,1),sem(statShufOpt,1), 'lineProps', 'k:')
shadedErrorBar(-100:10:1700, mean(runNormOpt,1),sem(runNormOpt,1), 'lineProps', 'r')
shadedErrorBar(-100:10:1700, mean(runShufOpt,1),sem(runShufOpt,1), 'lineProps', 'r:')

figure, hold on
plot(-100:10:1700, mean(statNormOpt,1),'k')
plot(-100:10:1700, mean(statShufOpt,1), 'k:')
plot(-100:10:1700, mean(runNormOpt,1), 'r')
plot(-100:10:1700, mean(runShufOpt,1), 'r:')



%% mean normal vs shuffled spike counts, gamma=0


statNormZero=[];
statShufZero=[];
runNormZero=[];
runShufZero=[];

for isession = 1:5
    statNormZero=cat(1,statNormZero,mean(cat(1,session(isession).stat.GammaZero.perm.meanPerf),1));
    statShufZero=cat(1,statShufZero,mean(cat(1,session(isession).stat.shuffleGammaZero.perm.meanPerf),1));
    runNormZero=cat(1,runNormZero,mean(cat(1,session(isession).run.GammaZero.perm.meanPerf),1));
    runShufZero=cat(1,runShufZero,mean(cat(1,session(isession).run.shuffleGammaZero.perm.meanPerf),1));
end


figure, hold on
shadedErrorBar(-100:10:1700, mean(statNormZero,1),sem(statNormZero,1), 'lineProps', 'k')
shadedErrorBar(-100:10:1700, mean(statShufZero,1),sem(statShufZero,1), 'lineProps', 'k:')
shadedErrorBar(-100:10:1700, mean(runNormZero,1),sem(runNormZero,1), 'lineProps', 'r')
shadedErrorBar(-100:10:1700, mean(runShufZero,1),sem(runShufZero,1), 'lineProps', 'r:')

title('Gamma = 0')



%% mean normal vs shuffled spike counts, gamma optimised

statNormOpt=[];
statShufOpt=[];
runNormOpt=[];
runShufOpt=[];

for isession = 1:5
    statNormOpt=cat(1,statNormOpt,mean(cat(1,session(isession).stat.GammaOpt.perm.meanPerf),1));
    statShufOpt=cat(1,statShufOpt,mean(cat(1,session(isession).stat.shuffleGammaOpt.perm.meanPerf),1));
    runNormOpt=cat(1,runNormOpt,mean(cat(1,session(isession).run.GammaOpt.perm.meanPerf),1));
    runShufOpt=cat(1,runShufOpt,mean(cat(1,session(isession).run.shuffleGammaOpt.perm.meanPerf),1));
end


figure, hold on
shadedErrorBar(-100:10:1700, mean(statNormOpt,1),sem(statNormOpt,1), 'lineProps', 'k')
shadedErrorBar(-100:10:1700, mean(statShufOpt,1),sem(statShufOpt,1), 'lineProps', 'k:')
shadedErrorBar(-100:10:1700, mean(runNormOpt,1),sem(runNormOpt,1), 'lineProps', 'r')
shadedErrorBar(-100:10:1700, mean(runShufOpt,1),sem(runShufOpt,1), 'lineProps', 'r:')
title('Gamma optimised')

%% mean normal vs shuffled spike counts, gamma=0


statNormZero=[];
statShufZero=[];
runNormZero=[];
runShufZero=[];

for isession = 1:5
    statNormZero=cat(1,statNormZero,mean(cat(1,session(isession).stat.GammaZero.perm.meanPerf),1));
    statShufZero=cat(1,statShufZero,mean(cat(1,session(isession).stat.shuffleGammaZero.perm.meanPerf),1));
    runNormZero=cat(1,runNormZero,mean(cat(1,session(isession).run.GammaZero.perm.meanPerf),1));
    runShufZero=cat(1,runShufZero,mean(cat(1,session(isession).run.shuffleGammaZero.perm.meanPerf),1));
end


figure, hold on
shadedErrorBar(-100:10:1700, mean(statNormZero,1),sem(statNormZero,1), 'lineProps', 'k')
shadedErrorBar(-100:10:1700, mean(statShufZero,1),sem(statShufZero,1), 'lineProps', 'k:')
shadedErrorBar(-100:10:1700, mean(runNormZero,1),sem(runNormZero,1), 'lineProps', 'r')
shadedErrorBar(-100:10:1700, mean(runShufZero,1),sem(runShufZero,1), 'lineProps', 'r:')

title('Gamma = 0')


%% difference by session

statDiff = (statShufZero)-(statNormZero);
runDiff = (runShufZero)-(runNormZero);

for isesh = 1:5
subplot(1,5,isesh), hold on
plot(-100:10:1700, statShufZero(isesh,:),'k:')
plot(-100:10:1700, statNormZero(isesh,:),'k')

plot(-100:10:1700, runShufZero(isesh,:),'r:')
plot(-100:10:1700, runNormZero(isesh,:),'r')
end

%% cumulative decoding

statGammaOpt = [];
statGammaOne = [];
statShuffleGammaOpt = [];

runGammaOpt = [];
runGammaOne = [];
runShuffleGammaOpt = [];

for isession = 1:5
    statGammaOpt = cat(1,statGammaOpt, mean(cat(1,session(isession).stat.GammaOpt.perm.meanPerf),1));
    statGammaOne = cat(1,statGammaOne, mean(cat(1,session(isession).stat.GammaOne.perm.meanPerf),1));
    statShuffleGammaOpt = cat(1,statShuffleGammaOpt, mean(cat(1,session(isession).stat.shuffleGammaOpt.perm.meanPerf),1));
    
    runGammaOpt = cat(1,runGammaOpt, mean(cat(1,session(isession).run.GammaOpt.perm.meanPerf),1));
    runGammaOne = cat(1,runGammaOne, mean(cat(1,session(isession).run.GammaOne.perm.meanPerf),1));
    runShuffleGammaOpt = cat(1,runShuffleGammaOpt, mean(cat(1,session(isession).run.shuffleGammaOpt.perm.meanPerf),1));
end


figure, hold on
shadedErrorBar(-190:10:1800, mean(statGammaOpt,1), sem(statGammaOpt,1), 'LineProps', 'k')
shadedErrorBar(-190:10:1800, mean(runGammaOpt,1), sem(runGammaOpt,1), 'LineProps', 'r')


diffArray = runGammaOpt-statGammaOpt;
figure, hold on
shadedErrorBar(-190:10:1800, mean(diffArray,1),sem(diffArray,1), 'lineProps', 'm')
plot([-200 1800], [0 0], 'k')


%% plot stat vs run
figure, 
subplot(131), hold on
shadedErrorBar(-190:10:1800, mean(statGammaOpt,1), sem(statGammaOpt,1), 'LineProps', 'k')
shadedErrorBar(-190:10:1800, mean(runGammaOpt,1), sem(runGammaOpt,1), 'LineProps', 'r')
title('Gamma Opt')

subplot(132), hold on
shadedErrorBar(-190:10:1800, mean(statGammaOne,1), sem(statGammaOne,1), 'LineProps', 'k')
shadedErrorBar(-190:10:1800, mean(runGammaOne,1), sem(runGammaOne,1), 'LineProps', 'r')
title('Gamma One')

subplot(133), hold on
shadedErrorBar(-190:10:1800, mean(statShuffleGammaOpt,1), sem(statShuffleGammaOpt,1), 'LineProps', 'k')
shadedErrorBar(-190:10:1800, mean(runShuffleGammaOpt,1), sem(runShuffleGammaOpt,1), 'LineProps', 'r')
title('Shuffle, Gamma Opt')


%% compare each decoding case within locomotor state

figure
subplot(121), hold on
shadedErrorBar(-190:10:1800, mean(statGammaOpt,1), sem(statGammaOpt,1), 'LineProps', 'k')
shadedErrorBar(-190:10:1800, mean(statGammaOne,1), sem(statGammaOne,1), 'LineProps', 'k:')
shadedErrorBar(-190:10:1800, mean(statShuffleGammaOpt,1), sem(statShuffleGammaOpt,1), 'LineProps', 'k--')

subplot(122), hold on
shadedErrorBar(-190:10:1800, mean(runGammaOpt,1), sem(runGammaOpt,1), 'LineProps', 'r')
shadedErrorBar(-190:10:1800, mean(runGammaOne,1), sem(runGammaOne,1), 'LineProps', 'r:')
shadedErrorBar(-190:10:1800, mean(runShuffleGammaOpt,1), sem(runShuffleGammaOpt,1), 'LineProps', 'r--')



%% just gamma =1

statNormOpt=[];
runNormOpt=[];

for isession = 1:5
    statNormOpt=cat(1,statNormOpt,mean(cat(1,session(isession).stat.GammaOne.perm.meanPerf),1));
    runNormOpt=cat(1,runNormOpt,mean(cat(1,session(isession).run.GammaOne.perm.meanPerf),1));
end


figure, hold on
shadedErrorBar(-150:10:1750, mean(statNormOpt,1),sem(statNormOpt,1), 'lineProps', 'k')
shadedErrorBar(-150:10:1750, mean(runNormOpt,1),sem(runNormOpt,1), 'lineProps', 'r')
title('Gamma=1')
yline(1/6)



diffArray = runNormOpt-statNormOpt;
figure, hold on
shadedErrorBar(-150:10:1750, mean(diffArray,1),sem(diffArray,1), 'lineProps', 'm')
yline(0)




%% rm anova

dec_stat = statGammaOne';
dec_run = runGammaOne';
decVec = cat(1,dec_stat(:), dec_run(:));
stateVec = categorical(cat(1,repelem(1, numel(dec_stat(:)))', repelem(2,numel(dec_run(:)))'));
timeVec = cat(1,repmat(1:191, 1,size(dec_stat,2))',repmat(1:191, 1,size(dec_run,2))');
subjVec = cat(1,repelem(1:size(dec_stat,2),1, 191)', repelem(1:size(dec_run,2),1, 191)');


[p,tbl,stats,terms] = anovan(decVec,{timeVec,stateVec,subjVec},'model',2,'random',3,'varnames',{'Time','State','Subj'});