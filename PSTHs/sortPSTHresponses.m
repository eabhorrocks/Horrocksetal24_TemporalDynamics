%% sort psth responses


sessionTags = {'M22027', '20220517';...
    'M22029', '20220607';...
    'M22031', '20220524';
    'M22032', '20220621';
    'M22033', '20220706'};

dataDir = 'C:\Users\edward.horrocks\Documents\Code\V1Dynamics\Data\basic_111022';


for isession = 1:5%size(sessionTags,1)
    fname = [sessionTags{isession,1},'_', sessionTags{isession,2},'_psth3rds.mat'];

    load(fullfile(dataDir,fname))

    session(isession).units = units;

end

%%

allUnits = cat(1,session([1 2 4 5]).units);
goodUnits = allUnits([allUnits.isi_viol]<=0.1...
    & [allUnits.amplitude_cutoff]<=0.1 & [allUnits.amplitude]>=50 & [allUnits.firing_rate]>=0);

%%

for isession = 1:5
    sesh_nUnits(isession) = numel(session(isession).units);
    sesh_nGood(isession) = numel(session(isession).units([session(isession).units.isi_viol]<=0.1...
    & [session(isession).units.amplitude_cutoff]<=0.1 & [session(isession).units.amplitude]>=50 &...
    [session(isession).units.firing_rate]>=0));
end

%% get tau rise time for cell-type classificaiton

for iunit = 1:numel(goodUnits)
    [ccg, t] = CCG(goodUnits(iunit).spike_times, ones(size(goodUnits(iunit).spike_times)),...
    'binSize', 0.0005, 'duration', 0.1,'norm', 'rate');
    goodUnits(iunit).acg = ccg;
    fit_params_out = fit_ACG(ccg,false);

    goodUnits(iunit).tau_rise = fit_params_out.acg_tau_rise;
end

%%
narrow_idx = find([goodUnits.duration]<=0.45);
wide_idx = find([goodUnits.duration]>0.45 & [goodUnits.tau_rise]>6);
pyr_idx = find(~ismember(1:numel(goodUnits), [narrow_idx,wide_idx]));

[goodUnits.cellType] = deal(nan);
[goodUnits(pyr_idx).cellType] = deal(1);
[goodUnits(narrow_idx).cellType] = deal(2);
[goodUnits(wide_idx).cellType] = deal(3);

%%
stat.speedVec = [];
stat.stateVec = [];
stat.psthVec = [];
stat.zpsthVec = [];
stat.reliVec = [];
stat.respVec = [];
stat.excVec =  [];
stat.suppVec  = [];
stat.onsetLatVec = [];
stat.offsetLatVec = [];
stat.peakLatVec = [];
stat.peak90LatVec = [];
stat.peakAbsZVec = [];
stat.meanFRVec = [];
stat.rangeFRVec = [];
stat.susIndexVec = [];
stat.baselineFRVec = [];
stat.waveformdur = [];
stat.cellType=[];

run.speedVec = [];
run.stateVec = [];
run.psthVec = [];
run.zpsthVec = [];
run.reliVec = [];
run.respVec = [];
run.excVec =  [];
run.suppVec  = [];
run.onsetLatVec = [];
run.offsetLatVec = [];
run.peakLatVec = [];
run.peak90LatVec = [];
run.peakAbsZVec = [];
run.meanFRVec = [];
run.rangeFRVec = [];
run.susIndexVec = [];
run.baselineFRVec = [];
run.waveformdur = [];
run.cellType=[];


statSmall.speedVec = [];
statSmall.statSmallVec = [];
statSmall.psthVec = [];
statSmall.zpsthVec = [];
statSmall.reliVec = [];
statSmall.respVec = [];
statSmall.excVec =  [];
statSmall.suppVec  = [];
statSmall.onsetLatVec = [];
statSmall.offsetLatVec = [];
statSmall.peakLatVec = [];
statSmall.peak90LatVec = [];
statSmall.peakAbsZVec = [];
statSmall.meanFRVec = [];
statSmall.rangeFRVec = [];
statSmall.susIndexVec = [];
statSmall.baselineFRVec = [];
statSmall.waveformdur = [];


statBig.speedVec = [];
statBig.statBigeVec = [];
statBig.psthVec = [];
statBig.zpsthVec = [];
statBig.reliVec = [];
statBig.respVec = [];
statBig.excVec =  [];
statBig.suppVec  = [];
statBig.onsetLatVec = [];
statBig.offsetLatVec = [];
statBig.peakLatVec = [];
statBig.peak90LatVec = [];
statBig.peakAbsZVec = [];
statBig.meanFRVec = [];
statBig.rangeFRVec = [];
statBig.susIndexVec = [];
statBig.baselineFRVec = [];
statBig.waveformdur = [];
tic
for iunit= 1:numel(goodUnits)

    for ispeed = 1:6

        % psths while stationary
        stat.speedVec = cat(1, stat.speedVec, ispeed);
        stat.stateVec = cat(1, stat.stateVec, 1);
        stat.psthVec = cat(1, stat.psthVec, goodUnits(iunit).stat_psth(ispeed).psth);
        stat.zpsthVec = cat(1, stat.zpsthVec, goodUnits(iunit).stat_psth(ispeed).zpsth);


        stat.reliVec = cat(1, stat.reliVec, goodUnits(iunit).stat_psth(ispeed).reliability_z);
        stat.excVec = cat(1, stat.excVec, goodUnits(iunit).stat_psth(ispeed).excited);
        stat.suppVec = cat(1, stat.suppVec, goodUnits(iunit).stat_psth(ispeed).suppressed);

        stat.onsetLatVec = cat(1, stat.onsetLatVec, goodUnits(iunit).stat_psth(ispeed).onsetLatency);
        stat.offsetLatVec = cat(1, stat.offsetLatVec, goodUnits(iunit).stat_psth(ispeed).offsetLatency);
        stat.peakLatVec = cat(1, stat.peakLatVec, goodUnits(iunit).stat_psth(ispeed).peakLatency);
        stat.peak90LatVec = cat(1, stat.peak90LatVec, goodUnits(iunit).stat_psth(ispeed).peak90Latency);

        stat.peakAbsZVec = cat(1, stat.peakAbsZVec, goodUnits(iunit).stat_psth(ispeed).peakAbsZ);
        stat.meanFRVec = cat(1, stat.meanFRVec, goodUnits(iunit).stat_psth(ispeed).meanFR);
        stat.rangeFRVec = cat(1, stat.rangeFRVec, goodUnits(iunit).stat_psth(ispeed).rangeFR);
        stat.susIndexVec = cat(1, stat.susIndexVec, goodUnits(iunit).stat_psth(ispeed).susIdx);
        stat.baselineFRVec = cat(1, stat.baselineFRVec, goodUnits(iunit).stat_blank.mu);
        stat.waveformdur = cat(1, stat.waveformdur, goodUnits(iunit).duration);
%         stat.cellType = cat(1, stat.cellType, goodUnits(iunit).cellType);


        % psths while running
        run.speedVec = cat(1, run.speedVec, ispeed);
        run.stateVev = cat(1, run.stateVec, 2);
        run.psthVec = cat(1, run.psthVec, goodUnits(iunit).run_psth(ispeed).psth);
        run.zpsthVec = cat(1, run.zpsthVec, goodUnits(iunit).run_psth(ispeed).zpsth);


        run.reliVec = cat(1, run.reliVec, goodUnits(iunit).run_psth(ispeed).reliability_z);
        run.excVec = cat(1, run.excVec, goodUnits(iunit).run_psth(ispeed).excited);
        run.suppVec = cat(1, run.suppVec, goodUnits(iunit).run_psth(ispeed).suppressed);

        run.onsetLatVec = cat(1, run.onsetLatVec, goodUnits(iunit).run_psth(ispeed).onsetLatency);
        run.offsetLatVec = cat(1, run.offsetLatVec, goodUnits(iunit).run_psth(ispeed).offsetLatency);
        run.peakLatVec = cat(1, run.peakLatVec, goodUnits(iunit).run_psth(ispeed).peakLatency);
        run.peak90LatVec = cat(1, run.peak90LatVec, goodUnits(iunit).run_psth(ispeed).peak90Latency);

        run.peakAbsZVec = cat(1, run.peakAbsZVec, goodUnits(iunit).run_psth(ispeed).peakAbsZ);
        run.meanFRVec = cat(1, run.meanFRVec, goodUnits(iunit).run_psth(ispeed).meanFR);
        run.rangeFRVec = cat(1, run.rangeFRVec, goodUnits(iunit).run_psth(ispeed).rangeFR);
        run.susIndexVec = cat(1, run.susIndexVec, goodUnits(iunit).run_psth(ispeed).susIdx);
        run.baselineFRVec = cat(1, run.baselineFRVec, goodUnits(iunit).run_blank.mu);
        run.waveformdur = cat(1, run.waveformdur, goodUnits(iunit).duration);
%         run.cellType = cat(1, run.cellType, goodUnits(iunit).cellType);



        % psths while stationary w/ small pupil
        statSmall.speedVec = cat(1, statSmall.speedVec, ispeed);
        statSmall.statSmallVec = cat(1, statSmall.statSmallVec, 1);
        statSmall.psthVec = cat(1, statSmall.psthVec, goodUnits(iunit).statSmall_psth(ispeed).psth);
        statSmall.zpsthVec = cat(1, statSmall.zpsthVec, goodUnits(iunit).statSmall_psth(ispeed).zpsth);


        statSmall.reliVec = cat(1, statSmall.reliVec, goodUnits(iunit).statSmall_psth(ispeed).reliability_z);
        statSmall.excVec = cat(1, statSmall.excVec, goodUnits(iunit).statSmall_psth(ispeed).excited);
        statSmall.suppVec = cat(1, statSmall.suppVec, goodUnits(iunit).statSmall_psth(ispeed).suppressed);

        statSmall.onsetLatVec = cat(1, statSmall.onsetLatVec, goodUnits(iunit).statSmall_psth(ispeed).onsetLatency);
        statSmall.offsetLatVec = cat(1, statSmall.offsetLatVec, goodUnits(iunit).statSmall_psth(ispeed).offsetLatency);
        statSmall.peakLatVec = cat(1, statSmall.peakLatVec, goodUnits(iunit).statSmall_psth(ispeed).peakLatency);
        statSmall.peak90LatVec = cat(1, statSmall.peak90LatVec, goodUnits(iunit).statSmall_psth(ispeed).peak90Latency);

        statSmall.peakAbsZVec = cat(1, statSmall.peakAbsZVec, goodUnits(iunit).statSmall_psth(ispeed).peakAbsZ);
        statSmall.meanFRVec = cat(1, statSmall.meanFRVec, goodUnits(iunit).statSmall_psth(ispeed).meanFR);
        statSmall.rangeFRVec = cat(1, statSmall.rangeFRVec, goodUnits(iunit).statSmall_psth(ispeed).rangeFR);
        statSmall.susIndexVec = cat(1, statSmall.susIndexVec, goodUnits(iunit).statSmall_psth(ispeed).susIdx);
        statSmall.baselineFRVec = cat(1, statSmall.baselineFRVec, goodUnits(iunit).statSmall_blank.mu);
        statSmall.waveformdur = cat(1, statSmall.waveformdur, goodUnits(iunit).duration);


        % psths while stationary w/ big pupil
        statBig.speedVec = cat(1, statBig.speedVec, ispeed);
        statBig.statBigeVec = cat(1, statBig.statBigeVec, 1);
        statBig.psthVec = cat(1, statBig.psthVec, goodUnits(iunit).statBig_psth(ispeed).psth);
        statBig.zpsthVec = cat(1, statBig.zpsthVec, goodUnits(iunit).statBig_psth(ispeed).zpsth);

        statBig.reliVec = cat(1, statBig.reliVec, goodUnits(iunit).statBig_psth(ispeed).reliability_z);
        statBig.excVec = cat(1, statBig.excVec, goodUnits(iunit).statBig_psth(ispeed).excited);
        statBig.suppVec = cat(1, statBig.suppVec, goodUnits(iunit).statBig_psth(ispeed).suppressed);

        statBig.onsetLatVec = cat(1, statBig.onsetLatVec, goodUnits(iunit).statBig_psth(ispeed).onsetLatency);
        statBig.offsetLatVec = cat(1, statBig.offsetLatVec, goodUnits(iunit).statBig_psth(ispeed).offsetLatency);
        statBig.peakLatVec = cat(1, statBig.peakLatVec, goodUnits(iunit).statBig_psth(ispeed).peakLatency);
        statBig.peak90LatVec = cat(1, statBig.peak90LatVec, goodUnits(iunit).statBig_psth(ispeed).peak90Latency);

        statBig.peakAbsZVec = cat(1, statBig.peakAbsZVec, goodUnits(iunit).statBig_psth(ispeed).peakAbsZ);
        statBig.meanFRVec = cat(1, statBig.meanFRVec, goodUnits(iunit).statBig_psth(ispeed).meanFR);
        statBig.rangeFRVec = cat(1, statBig.rangeFRVec, goodUnits(iunit).statBig_psth(ispeed).rangeFR);
        statBig.susIndexVec = cat(1, statBig.susIndexVec, goodUnits(iunit).statBig_psth(ispeed).susIdx);
        statBig.baselineFRVec = cat(1, statBig.baselineFRVec, goodUnits(iunit).statBig_blank.mu);
        statBig.waveformdur = cat(1, statBig.waveformdur, goodUnits(iunit).duration);


    end
end
toc
%% get reliably responsive

rangeThresh = 3;
reliThrsh =-1.645;
zThresh = 3.29;

statGoodVals =stat.reliVec<=reliThrsh & stat.rangeFRVec>=rangeThresh & stat.peakAbsZVec>=zThresh;
runGoodVals = run.reliVec<=reliThrsh & run.rangeFRVec>=rangeThresh & run.peakAbsZVec>=zThresh;
bothGoodVals = statGoodVals & runGoodVals;

statSmallGoodVals =statSmall.reliVec<=reliThrsh & statSmall.rangeFRVec>=rangeThresh & statSmall.peakAbsZVec>=zThresh;
statBigGoodVals =statBig.reliVec<=reliThrsh & statBig.rangeFRVec>=rangeThresh & statBig.peakAbsZVec>=zThresh;

%% mcnemar test for good responses

                     % run good
%                  +         -
%             --------------------
%         +   |   1,1   |   1,2   |
% stat good   |-------------------           
%         -   |   2,1   |   2,2   |
%             --------------------
%                                       
%
%   x=[101 59; 121 33];

x(1,1) = sum(statGoodVals & runGoodVals);
x(1,2) = sum(statGoodVals & ~runGoodVals);
x(2,1) = sum(~statGoodVals & runGoodVals);
x(2,2) = sum(~statGoodVals & ~runGoodVals);

x
 p = mcnemar(x)

 

%% plot average PSTH responses

figure, hold on
shadedErrorBar(1:200, mean(stat.psthVec(statGoodVals,:),1), sem(stat.psthVec(statGoodVals,:),1))
shadedErrorBar(1:200, mean(run.psthVec(runGoodVals,:),1), sem(run.psthVec(runGoodVals,:),1), 'lineProps', {'Color', 'r'})

figure
plot(mean(run.psthVec(runGoodVals,:),1) - mean(stat.psthVec(statGoodVals,:),1))

figure, hold on
shadedErrorBar(1:200, mean(statSmall.psthVec(statSmallGoodVals,:),1), sem(statSmall.psthVec(statSmallGoodVals,:),1))
shadedErrorBar(1:200, mean(statBig.psthVec(statBigGoodVals,:),1), sem(statBig.psthVec(statBigGoodVals,:),1), 'lineProps', {'Color', 'r'})

figure
plot(mean(statBig.psthVec(statBigGoodVals,:),1) - mean(statSmall.psthVec(statSmallGoodVals,:),1))


figure
plot(-195:10:1795,mean(stat.psthVec(statGoodVals,:),1),'k')
 defaultAxesProperties(gca, false);
 ax = gca; ax.XTick = -200:200:1800;
 ylim([0 20])

 figure
plot(-195:10:1795,mean(run.psthVec(runGoodVals,:),1), 'r'), hold on
% plot(1:200,mean(stat.psthVec(statGoodVals,:),1), 'k')
 defaultAxesProperties(gca, false);
 ax = gca; ax.XTick = -200:200:1800;
 ylim([0 20])



 %% plot small, big and running responses together

 figure, hold on
shadedErrorBar(-195:10:1795, mean(statSmall.psthVec(statSmallGoodVals,:),1), sem(statSmall.psthVec(statSmallGoodVals,:),1),'lineProps', {'Color', 'g'})
shadedErrorBar(-195:10:1795, mean(statBig.psthVec(statBigGoodVals,:),1), sem(statBig.psthVec(statBigGoodVals,:),1), 'lineProps', {'Color', 'c'})
shadedErrorBar(-195:10:1795, mean(run.psthVec(runGoodVals,:),1), sem(run.psthVec(runGoodVals,:),1), 'lineProps', {'Color', 'r'})
 ax = gca; ax.XTick = -200:200:1800;
 ylim([0 20])
  defaultAxesProperties(gca, true);

 %% plot average responses for each speed
speedcols = inferno(6);

 figure, hold on
 for ispeed = 1:6
     shadedErrorBar(1:200,mean(stat.psthVec(statGoodVals & stat.speedVec==ispeed,:),1),...
         sem(stat.psthVec(statGoodVals & stat.speedVec==ispeed,:),1),...
         'lineProps', {'Color', speedcols(ispeed,:)})
 end

 figure, hold on
 for ispeed = 1:6
     shadedErrorBar(1:200,mean(run.psthVec(runGoodVals & run.speedVec==ispeed,:),1),...
         sem(run.psthVec(runGoodVals & run.speedVec==ispeed,:),1),...
         'lineProps', {'Color', speedcols(ispeed,:)})
 end


%% do hierarchical clustering on stat responses

allPSTH_toUse = normalize(stat.psthVec(bothGoodVals,:),2,'range');
maxStretch = 10; % bins

nPSTH = size(allPSTH_toUse,1);
dm = nan([nPSTH]);

% nest loop version
tic
for ipsth1 = 1:nPSTH
    psth_temp = allPSTH_toUse(ipsth1,:);
    parfor ipsth2 = ipsth1:nPSTH
        dm_temp(:,ipsth2) = dtw(psth_temp,allPSTH_toUse(ipsth2,:), maxStretch);
    end
    dm(ipsth1,:) = dm_temp;
    dm_temp = [];
end
toc
full_dm = triu(dm)+triu(dm,1)';

% get dendrogram
Z_sub = linkage(full_dm, 'average');
% get optimal leaf ordering
tic
leafOrder = optimalleaforder(Z_sub,full_dm);
toc
[denHandle, leafnode_idx, outpermNodeOrder] = dendrogram(Z_sub, nPSTH, 'reorder', leafOrder, 'orientation', 'right');


%% generate new allPSTH array...

% outpermNodeOrder is the order of the leaf nodes...
% leafnode_idx is the label for individual responses to each leaf node

newPSTHorder = [];
speedVecReordered = [];

tic
for ileaf = 1:numel(outpermNodeOrder)
    idx = find(leafnode_idx==outpermNodeOrder(ileaf));
    newPSTHorder = vertcat(newPSTHorder, idx);
end
toc

allStat = normalize(stat.psthVec(bothGoodVals,:),2,'range');
allRun = normalize(run.psthVec(bothGoodVals,:),2,'range');



allStatReordered = allStat(newPSTHorder,:);
allRunReordered = allRun(newPSTHorder,:);

speedVec = stat.speedVec(bothGoodVals);
speedVecReordered = speedVec(newPSTHorder);
waveformVec = stat.waveformdur(bothGoodVals);
waveformVecReordered = waveformVec(newPSTHorder);
cellTypeVec = stat.cellType(bothGoodVals);
cellTypeVecReordered = cellTypeVec(newPSTHorder);


figure
ax1 = subplot(1,3,1);
imagesc(allStatReordered)
defaultAxesProperties(gca, true)
ax = gca; ax.XTick = 1:20:200;
ax1.Colormap = crameri('bilbao');

ax2 = subplot(1,3,2)
imagesc(allRunReordered)
defaultAxesProperties(gca, true)
ax = gca; ax.XTick = 1:20:200;
ax2.Colormap = crameri('bilbao')

%figure
ax3 = subplot(1,3,3)
imagesc(allRunReordered-allStatReordered)
caxis([-1 1]),
% differences in the normalized PSTHs indicat that running doesn't simply
% scale responses. There is a change in shape.
%diffPSTH = flipud(allRunReordered)-flipud(allStatReordered);
ax3.Colormap = redblue;
defaultAxesProperties(gca, true)
ax = gca; ax.XTick = 1:20:200;

%% differences by speed

for ispeed=1:6
    figure
    subplot(1,2,1)
    imagesc(flipud(allStatReordered(speedVecReordered==ispeed,:)))
    subplot(1,2,2)
    imagesc(flipud(allRunReordered(speedVecReordered==ispeed,:)))
end
%
% for ispeed=1:6
%     figure
%     imagesc(allRunReordered(speedVecReordered==ispeed,:) - allStatReordered(speedVecReordered==ispeed,:))
%     caxis([-1 1]), colormap(redblue)
% end

%% plot average psths by speed
for ispeed=1:6
    figure, hold on
    shadedErrorBar(1:200, mean(stat.psthVec(statGoodVals & stat.speedVec==ispeed,:),1), sem(stat.psthVec(statGoodVals & stat.speedVec==ispeed,:),1))
    shadedErrorBar(1:200, mean(run.psthVec(runGoodVals & run.speedVec==ispeed,:),1), sem(run.psthVec(runGoodVals & run.speedVec==ispeed,:),1), 'lineProps', {'Color', 'r'})
    title(num2str(ispeed))
    ylim([0 22])
end

%% plot mean responses by cell type

cellTypes = {'pyramidal', 'narrow', 'wide'};
cellCols = {'r', 'b', 'c'};
figure
for itype = 1:3
    subplot(1,3,itype), hold on
      shadedErrorBar(-195:10:1795, mean(stat.psthVec(statGoodVals & stat.cellType==itype,:),1), sem(stat.psthVec(statGoodVals & stat.cellType==itype,:),1))
      shadedErrorBar(-195:10:1795, mean(run.psthVec(runGoodVals & run.cellType==itype,:),1), sem(run.psthVec(runGoodVals & run.cellType==itype,:),1), 'lineProps', {'Color', 'r'})
      title(cellTypes{itype})
      ax=gca; ax.XTick = -200:200:1800;
      xlim([-200 1800]), ax.YLim(1) = 0;
      defaultAxesProperties(gca,false)
end






%% clusters of all PSTHs, irrespective of state

% allPSTH = vertcat(stat.psthVec(statGoodVals,:), run.psthVec(runGoodVals,:));
% 
% allPSTH_toUse = normalize(allPSTH,2,'range');
% maxStretch = 10; % bins
% 
% nPSTH = size(allPSTH_toUse,1);
% dm = nan([nPSTH]);
% 
% % nest loop version
% tic
% for ipsth1 = 1:nPSTH
%     psth_temp = allPSTH_toUse(ipsth1,:);
%     parfor ipsth2 = ipsth1:nPSTH
%         dm_temp(:,ipsth2) = dtw(psth_temp,allPSTH_toUse(ipsth2,:), maxStretch);
%     end
%     dm(ipsth1,:) = dm_temp;
%     dm_temp = [];
% end
% toc
% full_dm = triu(dm)+triu(dm,1)';
% 
% % get dendrogram
% Z_sub = linkage(full_dm, 'average');
% % get optimal leaf ordering
% tic
% leafOrder = optimalleaforder(Z_sub,full_dm);
% toc
% [denHandle, leafnode_idx, outpermNodeOrder] = dendrogram(Z_sub, nPSTH, 'reorder', leafOrder, 'orientation', 'right');


%% plot clusters

subclust_idx = cluster(Z_sub,'cutoff', 1.15); % 1.15 for average works ~best

nSubClust = max(subclust_idx)

figure
h = histogram(subclust_idx,nSubClust);
idx = find(h.Values>20);
nHighCountClusters = numel(idx) % number of high count clusters
medCount = median(h.Values) % median cluster count


%% plot each cluster with original members and medoid

for i = 1:numel(idx)
    figure
    idx2 = find(subclust_idx==idx(i));
    sequences = {};
    for iseq =1:numel(idx2)
        sequences{iseq} = allPSTH_toUse(idx2(iseq),:);
    end
    hold off
    average = DBA(sequences);
    medidx = medoidIndex(sequences, maxStretch);
    average = sequences{medidx};
    plot(allPSTH_toUse(subclust_idx==idx(i),:)'), hold on
    plot(average, 'k', 'LineWidth', 2)
    
    plot(mean(allPSTH_toUse(subclust_idx==idx(i),:)), 'k--', 'LineWidth', 2)
    ylim([0 1])
    pause

end

%% plot examples

statresps = allStatReordered;
runresps = allRunReordered;

figure
subplot(121)
imagesc(statresps)
subplot(122)
imagesc(runresps)
colormap(crameri('bilbao'))

% 969, 971

 figure
for idx = 44%[613, 1143, 1133]
figure 
% hold off
plot(statresps(idx,:), 'k')
 hold on
plot(runresps(idx,:), 'r')
% title(num2str(idx));
end

% statOG = stat.psthVec(bothGoodVals,:);
% runOG = run.psthVec(bothGoodVals,:);
% 
% oldidx = newPSTHorder(idx);
% plot(statOG(oldidx,:), 'k')
% hold on
% plot(runOG(oldidx,:), 'r')
% title(num2str(idx))

% end


%% looking for a specific example here..

statOG = stat.psthVec(bothGoodVals,:);
runOG = run.psthVec(bothGoodVals,:);

% ex1 = 1559; ex2 = 1453
% alt: 14, 33, 201, 224, 245, 266, 303, 305, 436, 450, 606, 607, 613, 614,
% 624, 626, 631, 632, 633, 634, 636, 637, 638, 666, 684, *694, 795, 796,
% 800, 807, 844, 845, 857, 872, 935, *939, *983 (checked to 1k)
%idx = find(max(statOG,[],2)>10);
%idx = [694, 939, 983]
idx = [1559, 983]% x1453, 983]
figure
for i =1:numel(idx)
    figure
    hold off,
    plot(statOG(idx(i),:),'k'), hold on
    plot(runOG(idx(i),:),'r'),
    title(idx(i))
    pause
end
% reference on flipud plot is size(statOG,1)-idx+1 ?

% get location on sorted responses (they are flipped ud!)
nResp = size(allStatReordered,1);
idx1 = find(newPSTHorder==1559);
idx2 = find(newPSTHorder==983);
loc1 = idx1/nResp
loc2 = idx2/nResp
% statGoodPSTH = stat.psthVec(bothGoodVals,:);
% runGoodPSTH = run.psthVec(bothGoodVals,:);
% 
% for iresp= 1500:size(statGoodPSTH,1)
%     hold off
%     plot(statGoodPSTH(iresp,:), 'k'),
%     hold on
%     plot(runGoodPSTH(iresp,:), 'r')
%     title(num2str(iresp))
%     pause
% end
% 
% 
%% plot examples based on location in sorted resps



% sorted order: 1405, 1425, 1427+, 1437, 1449, 1478*, 1516, *1535, 1555*,
% 1560, 1611, 1614*, 1621*, 1632, 1719*
statOG = stat.psthVec(bothGoodVals,:);
runOG = run.psthVec(bothGoodVals,:);

allStatReordered = statOG(newPSTHorder,:);
allRunReordered = runOG(newPSTHorder,:);

figure
idx = [1478, 709]
for  i = 1:numel(idx)
    figure
    plot(allStatReordered(idx(i),:),'k'), hold on
    plot(allRunReordered(idx(i),:),'r')
    title(idx(i));
end

nResp = size(allStatReordered,1);
loc(1) = idx(1)/nResp
loc(2) = idx(2)/nResp

%%
figure
subplot(121)
imagesc(allStat(newPSTHorder,:))
subplot(122)
imagesc(allRun(newPSTHorder,:))
colormap(crameri('bilbao'))
