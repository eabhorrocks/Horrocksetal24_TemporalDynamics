%% analyse FA fit on sliding windows

sessionTags = {'M22027', '20220517';...
    'M22029', '20220607';...
    'M22031', '20220524';
    'M22032', '20220621';
    'M22033', '20220706'};

dataDir = 'C:\Users\edward.horrocks\Documents\Code\V1Dynamics\Data\basic_111022';

s = struct;

for isession = 1:5
    tic
    fname = [sessionTags{isession,1},'_', sessionTags{isession,2},'_FAoverTime.mat']; 

    load(fullfile(dataDir,fname))

    s(isession).session = session;
end

%% plot population metrics

allStatSim=[];
allRunSim=[];
allStatSV=[];
allRunSV=[];
allStatDims=[];
allRunDims=[];

allStatSim_shuf=[];
allRunSim_shuf=[];
allStatSV_shuf=[];
allRunSV_shuf=[];
allStatDims_shuf=[];
allRunDims_shuf=[];

ii=0;
for isession=1:5
    ii=ii+1;
    allStatSim(ii,:) = cellfun(@(x) x(1), s(isession).session.statLoadingSim);
    allRunSim(ii,:) = cellfun(@(x) x(1), s(isession).session.runLoadingSim);

    allStatSV(ii,:) = s(isession).session.statSV;
    allRunSV(ii,:) = s(isession).session.runSV;
    
    allStatDims(ii,:) = s(isession).session.statDims;
    allRunDims(ii,:) = s(isession).session.runDims;

    allStatSim_shuf(ii,:) = cellfun(@(x) x(1), s(isession).session.statShufLoadingSim);
    allRunSim_shuf(ii,:) = cellfun(@(x) x(1), s(isession).session.runShufLoadingSim);

    allStatSV_shuf(ii,:) = s(isession).session.statShufSV;
    allRunSV_shuf(ii,:) = s(isession).session.runShufSV;
    
    allStatDims_shuf(ii,:) = s(isession).session.statShufDims;
    allRunDims_shuf(ii,:) = s(isession).session.runShufDims;
end


figure
subplot(3,2,1)
% tiledlayout('flow')
% nexttile
hold on
shadedErrorBar(-100:10:1700, mean(allStatSV,1), sem(allStatSV,1),'lineProps','k')
shadedErrorBar(-100:10:1700, mean(allRunSV,1), sem(allRunSV,1),'lineProps','r')
% ylabel('Proportion of variance shared')
ax = gca; ax.XTick = -200:200:1800; xlim([-200 1800]), ylim([0 40])
defaultAxesProperties(gca, true)

subplot(3,2,3)
hold on
shadedErrorBar(-100:10:1700, mean(allStatDims,1), sem(allStatDims,1),'lineProps','k')
shadedErrorBar(-100:10:1700, mean(allRunDims,1), sem(allRunDims,1),'lineProps','r')
% ylabel('Dimensionality of shared variance')
ax = gca; ax.XTick = -200:200:1800; xlim([-200 1800]), ylim([1 6]);  ax.YTick = 1:6;
defaultAxesProperties(gca, true)

subplot(3,2,5)
hold on
shadedErrorBar(-100:10:1700, mean(allStatSim,1), sem(allStatSim,1),'lineProps','k')
shadedErrorBar(-100:10:1700, mean(allRunSim,1), sem(allRunSim,1),'lineProps','r')
% ylabel('Loading similarity of first shared dim')
ax = gca; ax.XTick = -200:200:1800; xlim([-200 1800]), ylim([0, 0.6])
defaultAxesProperties(gca, true)


% plot population metrics of shuffled

% figure
% tiledlayout('flow')
subplot(3,2,2)
hold on
shadedErrorBar(-100:10:1700, mean(allStatSV_shuf,1), sem(allStatSV_shuf,1),'lineProps','k')
shadedErrorBar(-100:10:1700, mean(allRunSV_shuf,1), sem(allRunSV_shuf,1),'lineProps','r')
% ylabel('Proportion of variance shared')
ax = gca; ax.XTick = -200:200:1800; xlim([-200 1800]), ylim([0 40])
defaultAxesProperties(gca, true)

subplot(3,2,4)
hold on
shadedErrorBar(-100:10:1700, mean(allStatDims_shuf,1), sem(allStatDims_shuf,1),'lineProps','k')
shadedErrorBar(-100:10:1700, mean(allRunDims_shuf,1), sem(allRunDims_shuf,1),'lineProps','r')
% ylabel('Dimensionality of shared variance')
ax = gca; ax.XTick = -200:200:1800; xlim([-200 1800]),  ylim([1 6]); ax.YTick = 1:6;
defaultAxesProperties(gca, true)

subplot(3,2,6)
hold on
shadedErrorBar(-100:10:1700, mean(allStatSim_shuf,1), sem(allStatSim_shuf,1),'lineProps','k')
shadedErrorBar(-100:10:1700, mean(allRunSim_shuf,1), sem(allRunSim_shuf,1),'lineProps','r')
% ylabel('Loading similarity of first shared dim')
ax = gca; ax.XTick = -200:200:1800; xlim([-200 1800]),ylim([0, 0.6])
defaultAxesProperties(gca, true)

%% compare metrics within state w/ and w/o noise corrs

%%% stat %%%
figure
subplot(311)
hold on
shadedErrorBar(-100:10:1700, mean(allStatSV,1), sem(allStatSV,1),'lineProps','k')
shadedErrorBar(-100:10:1700, mean(allStatSV_shuf,1), sem(allStatSV_shuf,1),'lineProps','k:')
title('Proportion of variance shared')
ax = gca; ax.XTick = -200:200:1800; xlim([-200 1800])
defaultAxesProperties(gca, false)

subplot(312)
hold on
shadedErrorBar(-100:10:1700, mean(allStatDims,1), sem(allStatDims,1),'lineProps','k')
shadedErrorBar(-100:10:1700, mean(allStatDims_shuf,1), sem(allStatDims_shuf,1),'lineProps','k:')
title('Dimensionality of shared variance')
ax = gca; ax.XTick = -200:200:1800; xlim([-200 1800])
defaultAxesProperties(gca, false)

subplot(313)
hold on
shadedErrorBar(-100:10:1700, mean(allStatSim,1), sem(allStatSim,1),'lineProps','k')
shadedErrorBar(-100:10:1700, mean(allStatSim_shuf,1), sem(allStatSim_shuf,1),'lineProps','k:')
title('Loading similarity of first shared dim')
ax = gca; ax.XTick = -200:200:1800; xlim([-200 1800])
defaultAxesProperties(gca, false)


%%% run %%%


figure
subplot(311)
hold on
shadedErrorBar(-100:10:1700, mean(allRunSV,1), sem(allRunSV,1),'lineProps','r')
shadedErrorBar(-100:10:1700, mean(allRunSV_shuf,1), sem(allRunSV_shuf,1),'lineProps','r:')
title('Proportion of variance shared')
ax = gca; ax.XTick = -200:200:1800; xlim([-200 1800])
defaultAxesProperties(gca, false)

subplot(312)
hold on
shadedErrorBar(-100:10:1700, mean(allRunDims,1), sem(allRunDims,1),'lineProps','r')
shadedErrorBar(-100:10:1700, mean(allRunDims_shuf,1), sem(allRunDims_shuf,1),'lineProps','r:')
title('Dimensionality of shared variance')
ax = gca; ax.XTick = -200:200:1800; xlim([-200 1800])
defaultAxesProperties(gca, false)

subplot(313)
hold on
shadedErrorBar(-100:10:1700, mean(allRunSim,1), sem(allRunSim,1),'lineProps','r')
shadedErrorBar(-100:10:1700, mean(allRunSim_shuf,1), sem(allRunSim_shuf,1),'lineProps','r:')
title('Loading similarity of first shared dim')
ax = gca; ax.XTick = -200:200:1800; xlim([-200 1800])
defaultAxesProperties(gca, false)



%% RM-ANOVAs

bv = -100:10:1700;
time = find(bv==-100):find(bv==1700);

statVals = allStatSim(:,time)';
runVals = allRunSim(:,time)';

statVals2 = allStatSim_shuf(:,time)';
runVals2 = allRunSim_shuf(:,time)';


vals = cat(1,statVals(:), runVals(:));
stateVec = categorical((cat(1,repelem(1, numel(statVals(:)))', repelem(2,numel(runVals(:)))')));
timeVec = categorical(cat(1,repmat(1:numel(time), 1,size(statVals,2))',repmat(1:numel(time), 1,size(runVals,2))'));
subjVec = categorical(cat(1,repelem(1:size(statVals,2),1, numel(time))', repelem(1:size(runVals,2),1, numel(time))'));

% effect  of state on shared variance (intact)
[p,tbl,stats,terms] = anovan(vals,{timeVec,stateVec,subjVec},'model',2,'random',3,'varnames',{'Time','State','Subj'});

vals = cat(1,statVals2(:), runVals2(:));
% effect of state on shared variance (disrupted)
[p,tbl,stats,terms] = anovan(vals,{timeVec,stateVec,subjVec},'model',2,'random',3,'varnames',{'Time','State','Subj'});



vals = cat(1,statVals(:), runVals(:), statVals2(:), runVals2(:));
corrVec = repelem(1:2,1,numel(subjVec));
stateVec2 = cat(1,stateVec,stateVec);
timeVec2 = cat(1,timeVec,timeVec);
subjVec2 = cat(1,subjVec,subjVec);

[p,tbl,stats,terms] = anovan(vals,{timeVec2,stateVec2,subjVec2, corrVec},'model','interaction','random',3,'varnames',{'Time','State','Subj','Corr'});


%% interaction plot

figure

propSVrun = allRunSV_shuf./allRunSV;
propSVstat = allStatSV_shuf./allStatSV;
subplot(311), hold on
plot([-200 1800], [0.5 0.5], 'k:')
shadedErrorBar(-100:10:1700, mean(propSVstat,1), sem(propSVstat,1),'lineProps','k')
shadedErrorBar(-100:10:1700, mean(propSVrun,1), sem(propSVrun,1),'lineProps','r')
ax = gca; ax.XTick = -200:200:1800; xlim([-200 1800])
ylim([0 1])
defaultAxesProperties(gca, true)



propDimsrun = (allRunDims_shuf-1)./(allRunDims-1);
propDimsstat = (allStatDims_shuf-1)./(allStatDims-1);
propDimsstat(isinf(propDimsstat))=nan;
subplot(312), hold on
plot([-200 1800], [0.5 0.5], 'k:')
shadedErrorBar(-100:10:1700, nanmean(propDimsstat,1), nansem(propDimsstat,1),'lineProps','k')
shadedErrorBar(-100:10:1700, nanmean(propDimsrun,1), nansem(propDimsrun,1),'lineProps','r')
ax = gca; ax.XTick = -200:200:1800; xlim([-200 1800])
ylim([0 1])
defaultAxesProperties(gca, true)


absRunDiff = abs(allRunSim_shuf - allRunSim);
absStatDiff = abs(allStatSim_shuf - allStatSim);

propSimrun = allRunSim_shuf./allRunSim;
propSimstat = allStatSim_shuf./allStatSim;
% deal with outliers
propSimrun(propSimrun>1)=1;
propSimstat(propSimstat>1)=1;
subplot(313), hold on
plot([-200 1800], [0.5 0.5], 'k:')
shadedErrorBar(-100:10:1700, nanmean(propSimstat,1), nansem(propSimstat,1),'lineProps','k')
shadedErrorBar(-100:10:1700, nanmean(propSimrun,1), nansem(propSimrun,1),'lineProps','r')
ax = gca; ax.XTick = -200:200:1800; xlim([-200 1800])
ylim([0 1])
defaultAxesProperties(gca, true)



