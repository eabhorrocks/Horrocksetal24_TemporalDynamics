%% analyse factor analysis (June 2023 version)

sessionTags = {'M22027', '20220517';...
    'M22029', '20220607';...
    'M22031', '20220524';
    'M22032', '20220621';
    'M22033', '20220706'};

dataDir = 'C:\Users\edward.horrocks\Documents\Code\V1Dynamics\Data\basic_111022';
dataDir = 'E:\V1Data\Data\v1_fromC24';
s = struct;

for isession = 1:5
    tic
    fname = [sessionTags{isession,1},'_', sessionTags{isession,2},'_FA_normal.mat'];

    load(fullfile(dataDir,fname))

    s(isession).session = session;
end


%% plot explained variance
 figure, hold on
for isession = 1:5
totalSV(isession) = s(isession).session.s.SV;
plot(cumsum(s(isession).session.s.propSharedVariance),'LineStyle','-','Marker','.')
%plot(s(isession).session.s.propSharedVariance,'LineStyle','-','Marker','.')

end
legend({'1,','2','3','4','5'})
%%

figure, hold on
for isession = 1:5
t(isession).sv = cumsum(s(isession).session.s.propSharedVariance);
t(isession).nsv = s(isession).session.s.propSharedVariance;

end

allSV = padcat(t.nsv,2);
allSV=allSV(:,1:5);
% allSV = allSV.*totalSV;

figure
errorbar(1:size(allSV,1),nanmean(allSV,2), nansem(allSV,2))

figure
shadedErrorBar(1:size(allSV,1),nanmean(allSV,2), nansem(allSV,2))

%
figure
errorbar(1:size(allSV,1),nanmean(allSV,2), nansem(allSV,2),'lineStyle','-','Color','k','Marker','.')

% 
% figure
% shadedErrorBar(1:size(allSV,1),nanmean(allSV,2), nansem(allSV,2),'lineProps','k.-')
% errorbar(1:size(allSV,1),nanmean(allSV,2), nansem(allSV,2),'lineStyle','-','Color','k','Marker','.')

%% get distance metrics with no weights

bin_N = 1; % bin every N elements
binSize = 10; % bin size in ms
binVector = -200:binSize*bin_N:(1800-binSize*bin_N);
timeBinVector = round(binVector,2); % time associated with each index (1st edge of bin)
nWins = 10; % number of consecutive windows required to reach steady state

for isession = 1:5

    % set parameters
    %weights = s(isession).session.s.propSharedVariance(1:s(isession).session.s.qOpt)'; % weights are frac. of shared variance explained by each component
    weights = repelem(1, 1, s(isession).session.s.qOpt);
    cond = s(isession).session.s.cond;

    % onset trajectories
    initialSSIdx = 1:find(timeBinVector==0)-1;
    finalSSIdx = find(timeBinVector==500):find(timeBinVector==1000)-1;



    for icond = 1:numel(cond)
        neuralTrajectory = cond(icond).meanTrajectory(:,1:s(isession).session.s.qOpt);
        output_onset(icond) = calcPopulationGeoMetrics(neuralTrajectory, weights, timeBinVector, initialSSIdx, finalSSIdx, nWins);
    end

    % offset trajectories
    initialSSIdx = find(timeBinVector==500):find(timeBinVector==1000)-1;
    finalSSIdx = find(timeBinVector==1500):find(timeBinVector==1800-binSize*bin_N);

    for icond = 1:numel(cond)
        neuralTrajectory = cond(icond).meanTrajectory(:,1:s(isession).session.s.qOpt);
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

   % statPerf  =cat(1, statPerf, s(isession).session.stat.meanPerf);
   % runPerf  =cat(1, runPerf, s(isession).session.run.meanPerf);

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

%% plot distance metrics as bar charts

% onset
figure,
subplot(121), hold on
bar(1:2, [mean(statDD(:)), mean(runDD(:))])
errorbar(1:2, mean([statDD(:),runDD(:)],1), sem([statDD(:),runDD(:)],1))
ylabel('Direct Distance')
defaultAxesProperties(gca, false)

subplot(122), hold on
bar(1:2, [mean(statCD(:)), mean(runCD(:))])
errorbar(1:2, mean([statCD(:),runCD(:)],1), sem([statCD(:),runCD(:)],1))
ylabel('Distance Travelled')
defaultAxesProperties(gca, false)

% offset
figure,
subplot(121), hold on
bar(1:2, [mean(statDDoff(:)), mean(runDDoff(:))])
errorbar(1:2, mean([statDDoff(:),runDDoff(:)],1), sem([statDDoff(:),runDDoff(:)],1))
ylabel('Direct Distance')
defaultAxesProperties(gca, false)

subplot(122), hold on
bar(1:2, [mean(statCDoff(:)), mean(runCDoff(:))])
errorbar(1:2, mean([statCDoff(:),runCDoff(:)],1), sem([statCDoff(:),runCDoff(:)],1))
ylabel('Distance Travelled')
defaultAxesProperties(gca, false)

%% box plot versions

statVals = statDD;
runVals = runDD;

cols_cellArray = {[0 0 0], [1 0 1]};
allVals = [statVals(:), runVals(:)];
xvals = [1, 2];

figure, hold on
for ival = 1:numel(statVals)
plot([1 2], [statVals(ival), runVals(ival)],'Color',[.7 .7 .7],'LineStyle','-',...
    'Marker', 'o');
end
boxplot(allVals,'Colors','kr')
ylabel('onsetDD')
ylim([0 16])
defaultAxesProperties(gca, true)


statVals = statDDoff;
runVals = runDDoff;

cols_cellArray = {[0 0 0], [1 0 1]};
allVals = [statVals(:), runVals(:)];
xvals = [1, 2];

figure, hold on
for ival = 1:numel(statVals)
plot([1 2], [statVals(ival), runVals(ival)],'Color',[.7 .7 .7],'LineStyle','-',...
    'Marker', 'o');
end
boxplot(allVals,'Colors','kr')
ylabel('offsetDD')
ylim([0 16])
defaultAxesProperties(gca, true)

%%% disttance travelled %%%
statVals = statCD;
runVals = runCD;

cols_cellArray = {[0 0 0], [1 0 1]};
allVals = [statVals(:), runVals(:)];
xvals = [1, 2];

figure, hold on
for ival = 1:numel(statVals)
plot([1 2], [statVals(ival), runVals(ival)],'Color',[.7 .7 .7],'LineStyle','-',...
    'Marker', 'o');
end
boxplot(allVals,'Colors','kr')
ylabel('dist travelled')
ylim([0 45])
defaultAxesProperties(gca, true)


statVals = statCDoff;
runVals = runCDoff;

cols_cellArray = {[0 0 0], [1 0 1]};
allVals = [statVals(:), runVals(:)];
xvals = [1, 2];

figure, hold on
for ival = 1:numel(statVals)
plot([1 2], [statVals(ival), runVals(ival)],'Color',[.7 .7 .7],'LineStyle','-',...
    'Marker', 'o');
end
boxplot(allVals,'Colors','kr')
ylabel('dist travlled offsett')
ylim([0 45])
defaultAxesProperties(gca, true)
%% distribution plot versions
statVals = statDD;
runVals = runDD;

cols_cellArray = {[0 0 0], [1 0 1]};
allVals = {statVals, runVals};
xvals = [1, 2];

% violin plot
figure
distributionPlot(allVals,'globalNorm',3,'histOpt',1,'divFactor',3,...
    'xValues', xvals, 'addSpread', 0, 'distWidth', 0.8, 'Color', cols_cellArray,...
    'showMM',6)
ylabel('DD onset')
ylim([0 16])
defaultAxesProperties(gca, true)
statVals = statCD;
runVals = runCD;

allVals = {statVals, runVals};
xvals = [1, 2];

% violin plot
figure
distributionPlot(allVals,'globalNorm',3,'histOpt',1,'divFactor',3,...
    'xValues', xvals, 'addSpread', 0, 'distWidth', 0.8, 'Color', cols_cellArray,...
    'showMM',6)
ylabel('CD onset')
ylim([0 45])
defaultAxesProperties(gca, true)


statVals = statDDoff;
runVals = runDDoff;

allVals = {statVals, runVals};
xvals = [1, 2];

% violin plot
figure
distributionPlot(allVals,'globalNorm',3,'histOpt',1,'divFactor',3,...
    'xValues', xvals, 'addSpread', 0, 'distWidth', 0.8, 'Color', cols_cellArray,...
    'showMM',6)
ylabel('DD offset')
ylim([0 16])
defaultAxesProperties(gca, true)
statVals = statCDoff;
runVals = runCDoff;

allVals = {statVals, runVals};
xvals = [1, 2];

% violin plot
figure
distributionPlot(allVals,'globalNorm',3,'histOpt',1,'divFactor',3,...
    'xValues', xvals, 'addSpread', 0, 'distWidth', 0.8, 'Color', cols_cellArray,...
    'showMM',6)
ylabel('CD offset')
ylim([0 45])
defaultAxesProperties(gca, true)


%% stats for distance measures

statVals = statDDoff';
runVals = runDDoff';
valsVec = cat(1,statVals(:),runVals(:));

sessionVec = repelem(1:5,1,6)';
speedVec = repmat(1:6,1,5)';

stateVec = categorical(cat(1,repelem(1,numel(statVals),1), repelem(2,numel(runVals),1)));
speedVec = categorical(cat(1,speedVec,speedVec));
seshVec = categorical(cat(1,sessionVec, sessionVec));

tbl = table(valsVec,speedVec,stateVec,seshVec,...
    'VariableNames',{'vals','speed','state','sesh'});

f = 'vals ~ state + (1|speed) + (1|sesh)';

lme = fitlme(tbl, f, 'DummyVarCoding', 'reference')


[mean(statVals(:)), sem(statVals(:)); mean(runVals(:)), sem(runVals(:))]




%% speed and acceleration of mean trajectories


allStat=[];
allRun=[];

allStatAcc = [];
allRunAcc = [];

allStatJerk = [];
allRunJerk = [];


for isesh = 1:5
    qOpt = s(isesh).session.s.qOpt;
    weights = repelem(1, 1, qOpt);
    nFactors = size(s(1).session.s.loadings,2);

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

ylim([0 12])
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

%% plot by speed

figure
for ispeed = 1:6
    allStat = [];
    allRun = [];
    for isesh = 1:5
        allStat = cat(1,allStat, s(isesh).session.s.cond(ispeed).speed);
        allRun = cat(1,allRun, s(isesh).session.s.cond(ispeed+6).speed);
    end

    subplot(221), hold on
    plot(-190:10:1795, mean(allStat,1),'Color', speedcols(ispeed,:))
    subplot(223), hold on
    plot(-190:10:1795, mean(allRun,1),'Color', speedcols(ispeed,:))
end

subplot(221)
xlim([-200 1800])
ylim([0 12])
defaultAxesProperties(gca, true)
subplot(223)
xlim([-200 1800])
ylim([0 12])
defaultAxesProperties(gca, true)
%


for ispeed = 1:6
    allStat = [];
    allRun = [];
    for isesh = 1:5
        allStat = cat(1,allStat, s(isesh).session.s.cond(ispeed).speedAcc);
        allRun = cat(1,allRun, s(isesh).session.s.cond(ispeed+6).speedAcc);
    end

    subplot(222), hold on
    plot(-180:10:1795, mean(allStat*100,1),'Color', speedcols(ispeed,:))
    subplot(224), hold on
    plot(-180:10:1795, mean(allRun*100,1),'Color', speedcols(ispeed,:))
end

subplot(222)
xlim([-200 1800])
ylim([-150 200])
defaultAxesProperties(gca, true)
subplot(224)
xlim([-200 1800])
ylim([-150 200])
defaultAxesProperties(gca, true)


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


valsVec = offsetMaxSpdVec;

tbl = table(valsVec,speedVec,stateVec,seshVec,...
    'VariableNames',{'vals','speed','state','sesh'});

f = 'vals ~ state + (1|speed) + (1|sesh)';

lme = fitlme(tbl, f, 'DummyVarCoding', 'reference')


[mean(valsVec(1:30)), sem(valsVec(1:30)); mean(valsVec(31:60)), sem(valsVec(31:60))]



%% power spectrum of accelerations

freqBandThresh = 6;

for isesh = 1:5
    for icond = 1:12
        [s(isesh).session.s.cond(icond).accPower,f] = pspectrum(s(isesh).session.s.cond(icond).speedAcc,100,'FrequencyLimits',[0 20],'FrequencyResolution',2);
        s(isesh).session.s.cond(icond).relPower = s(isesh).session.s.cond(icond).accPower/mean(s(isesh).session.s.cond(icond).accPower);
        s(isesh).session.s.cond(icond).lowPower = mean(s(isesh).session.s.cond(icond).relPower(1:find(f<freqBandThresh,1,'last')));
        s(isesh).session.s.cond(icond).highPower = mean(s(isesh).session.s.cond(icond).relPower(find(f>freqBandThresh,1,'first'):end));

    end
end

allStat=[];
allRun=[];

allStatLow=[];
allRunLow=[];

allStatHigh=[];
allRunHigh=[];

for isesh = 1:5
    for icond = 1:6

        allStat = cat(2,allStat,s(isesh).session.s.cond(icond).relPower);
        allRun = cat(2,allRun,s(isesh).session.s.cond(icond+6).relPower);

        allStatLow = cat(1,allStatLow,s(isesh).session.s.cond(icond).lowPower);
        allRunLow = cat(1,allRunLow,s(isesh).session.s.cond(icond+6).lowPower);

        allStatHigh = cat(1,allStatHigh,s(isesh).session.s.cond(icond).highPower);
        allRunHigh = cat(1,allRunHigh,s(isesh).session.s.cond(icond+6).highPower);

    end
end

figure, hold on
shadedErrorBar(f, mean(allStat,2), sem(allStat,2));
shadedErrorBar(f, mean(allRun,2), sem(allRun,2),'lineProps','r');

figure
subplot(121), hold on
bar([1, 2], [mean(allStatLow),mean(allRunLow)]);
errorbar(1:2,[mean(allStatLow),mean(allRunLow)], [sem(allStatLow),sem(allRunLow)], 'lineStyle','none')
subplot(122), hold on
bar([1, 2], [mean(allStatHigh),mean(allRunHigh)])
errorbar(1:2,[mean(allStatHigh),mean(allRunHigh)], [sem(allStatHigh),sem(allRunHigh)],'lineStyle','none')



%%  plot power analysis  by stim speed

figure,
allStat=[];
allRun=[];

allStatLow=[];
allRunLow=[];

allStatHigh=[];
allRunHigh=[];

for icond = 1:6
    allStat=[];
    allRun=[];

    allStatLow=[];
    allRunLow=[];

    allStatHigh=[];
    allRunHigh=[];
    for isesh = 1:5


        allStat = cat(2,allStat,s(isesh).session.s.cond(icond).relPower);
        allRun = cat(2,allRun,s(isesh).session.s.cond(icond+6).relPower);

        allStatLow = cat(1,allStatLow,s(isesh).session.s.cond(icond).lowPower);
        allRunLow = cat(1,allRunLow,s(isesh).session.s.cond(icond+6).lowPower);

        allStatHigh = cat(1,allStatHigh,s(isesh).session.s.cond(icond).highPower);
        allRunHigh = cat(1,allRunHigh,s(isesh).session.s.cond(icond+6).highPower);

    end

    subplot(1,2,1), hold on
    plot(f, mean(allStat,2), 'Color', speedcols(icond,:))

    subplot(1,2,2), hold on
    plot(f, mean(allRun,2), 'Color', speedcols(icond,:))

end


%% analyse  'dynamic range' of trajectories


for isession = 1:5
    stat = s(isession).session.s.cond(1:6);
    run = s(isession).session.s.cond(7:12);
    weights = s(isession).session.s.propSharedVariance';
    weights=ones(size(weights));
    clear statWeight runWeight

    for iint = 1:200
        for ispeed1 =1:6
            for ispeed2 = 1:6
                statWeight(iint,ispeed1,ispeed2) = sqrt(sum(weights.*(stat(ispeed1).meanTrajectory(iint,:)-stat(ispeed2).meanTrajectory(iint,:)).^2));
                runWeight(iint,ispeed1,ispeed2) = sqrt(sum(weights.*(run(ispeed1).meanTrajectory(iint,:)-run(ispeed2).meanTrajectory(iint,:)).^2));
            end
        end

        statWeight(iint,:,:) = setUpperTri2NaN(squeeze(statWeight(iint,:,:)));
        runWeight(iint,:,:) = setUpperTri2NaN(squeeze(runWeight(iint,:,:)));

    end

    s(isession).session.stat.pairRangeDists = nanmean(statWeight,[2,3]);
    s(isession).session.run.pairRangeDists = nanmean(runWeight,[2,3]);

end

allStat=[];
allRun=[];
for isession =1:5
    allStat=cat(2,allStat,s(isession).session.stat.pairRangeDists);
    allRun=cat(2,allRun,s(isession).session.run.pairRangeDists);
end

figure, hold on
shadedErrorBar(-195:10:1795, mean(allStat,2), sem(allStat,2))
shadedErrorBar(-195:10:1795, mean(allRun,2), sem(allRun,2), 'lineProps', 'r')
xlim([-200 1800])
ax=gca; ax.XTick = -200:200:1800;
defaultAxesProperties(gca, true)
xlabel('Time'), ylabel('Mean distance between speed-mean trajectories')


%% stats on dynamic range

meanStat = mean(allStat(21:121,:),1);
meanRun = mean(allRun(21:121,:),1);

[h,p] = adtest(meanStat-meanRun)

[h,p] = ttest(meanStat(:),meanRun(:))



%% average decoding performance

% plot average decoding performance
bmp = -175:10:1775; bmp(end) = [];

allStat = [];
allRun = [];

for isession = 1:5
    allStat = cat(1,allStat, s(isession).session.stat.meanPerf);
    allRun = cat(1,allRun, s(isession).session.run.meanPerf);
end

figure,  hold on
plot([-0.2 1.8], [1/6 1/6], 'k:')
shadedErrorBar(bmp, mean(allStat,1), sem(allStat,1))
shadedErrorBar(bmp, mean(allRun,1), sem(allRun,1),'lineProps', 'r')
ax = gca; ax.XTick = -200:200:1800;
xlim([-200 1800])
defaultAxesProperties(gca,true)

% figure, hold on
% shadedErrorBar(bmp, mean(diffArray,1), sem(diffArray,1), 'lineProps','m')
% defaultAxesProperties(gca,false)
% ax = gca; ax.XTick = -200:200:1800;
% xlim([-200 1800])
% plot([-200 1800], [0 0], 'k:')

%%

dec_stat = allStat';
dec_run = allRun';
decVec = cat(1,dec_stat(:), dec_run(:));
stateVec = categorical(cat(1,repelem(1, numel(dec_stat(:)))', repelem(2,numel(dec_run(:)))'));
timeVec = cat(1,repmat(1:190, 1,size(dec_stat,2))',repmat(1:190, 1,size(dec_run,2))');
subjVec = cat(1,repelem(1:size(dec_stat,2),1, 190)', repelem(1:size(dec_run,2),1, 190)');


[p,tbl,stats,terms] = anovan(decVec,{timeVec,stateVec,subjVec},'model',2,'random',3,'varnames',{'Time','State','Subj'});


%% cross-condition tangling


for isesh = 1:5
    qOpt =s(isesh).session.s.qOpt;
    %qOpt = find(cumsum(s(isesh).session.s.propSharedVariance)>=0.85,1,'first');%6;% s(isesh).session.s.qOpt;
    for ispeed1 = 1:6
        for ispeed2 = 1:6

            idim = 1:qOpt;
            itime = 1:200;

            weights = s(isesh).session.s.propSharedVariance(1:qOpt);
            deltaT = 10;
            epsVal = 0.1;


            statTraj1 = s(isesh).session.cond(ispeed1).meanTrajectory(itime,idim);
            statTraj2 = s(isesh).session.cond(ispeed2).meanTrajectory(itime,idim);
            Qstat=[];
            for it = 2:numel(itime)
                for it2 = 2:numel(itime)

                    xt = statTraj1(itime(it),:);
                    xt2 = statTraj2(itime(it2),:);

                    xderiv_t = (statTraj1(itime(it),:) - statTraj1(itime(it-1),:))/deltaT;
                    xderiv_t2 = (statTraj2(itime(it2),:) - statTraj2(itime(it2-1),:))/deltaT;

                    Qstat(it,it2) = (norm(xderiv_t-xderiv_t2)^2) / (norm(xt-xt2)^2 + epsVal);

                end
            end


            runTraj1 = s(isesh).session.cond(ispeed1+6).meanTrajectory(itime,idim);
            runTraj2 = s(isesh).session.cond(ispeed2+6).meanTrajectory(itime,idim);
            Qrun=[];
            for it = 2:numel(itime)
                for it2 = 2:numel(itime)

                    xt = runTraj1(itime(it),:);
                    xt2 = runTraj2(itime(it2),:);
                    %
                    xderiv_t = (runTraj1(itime(it),:) - runTraj1(itime(it-1),:))/deltaT;
                    xderiv_t2 = (runTraj2(itime(it2),:) - runTraj2(itime(it2-1),:))/deltaT;


                    Qrun(it,it2) = (norm(xderiv_t-xderiv_t2)^2) / (norm(xt-xt2)^2 + epsVal);

                end
            end

            sesh(isesh).crossSpeed(ispeed1,ispeed2).Qstat = Qstat;
            sesh(isesh).crossSpeed(ispeed1,ispeed2).Qrun = Qrun;
        end



    end
end

%%
binSize = 10;
windowSize = 200/binSize;

for isesh = 1:5
    for ispeed1 = 1:6
        for ispeed2 = 1:6
            this_Qstat = sesh(isesh).crossSpeed(ispeed1,ispeed2).Qstat;
            this_Qrun = sesh(isesh).crossSpeed(ispeed1,ispeed2).Qrun;

            for itime = 1:180
                tqs = this_Qstat(itime:itime+windowSize, itime:itime+windowSize);
                tqr = this_Qrun(itime:itime+windowSize, itime:itime+windowSize);


                sesh(isesh).crossSpeed(ispeed1,ispeed2).statLocalTangling(itime) = ...
                    prctile(tqs(1+windowSize/2,:),90,'all');

                sesh(isesh).crossSpeed(ispeed1,ispeed2).runLocalTangling(itime) = ...
                    prctile(tqr(1+windowSize/2,:),90,'all');
            end
        end
    end
end


%% plot all combinations

figure
tiledlayout(6,6)
for ispeed1 = 1:6
    for ispeed2 = 1:6
        nexttile, hold on
        allStat=[];
        allRun=[];
        for isesh = 1:5
            allStat = cat(1,allStat,sesh(isesh).crossSpeed(ispeed1,ispeed2).statLocalTangling);
            allRun = cat(1,allRun,sesh(isesh).crossSpeed(ispeed1,ispeed2).runLocalTangling);
        end

        shadedErrorBar(-95:10:1695, mean(allStat,1), sem(allStat,1),'lineProps','k')
        shadedErrorBar(-95:10:1695, mean(allRun,1), sem(allRun,1),'lineProps','r')
        xlim([-200 1800])
        ax=gca; ax.XTick = -200:200:1800;
        defaultAxesProperties(gca, true)
        ylabel('Local Tangling')
        xlabel('Time (ms)')
        title([ispeed1, ispeed2])

    end
end

%% plot average of cross-combinations

allStat=[];
allRun=[];
for ispeed1 = 1:6
    for ispeed2 = 1:6

        if ispeed1 ~= ispeed2

            for isesh = 1:5
                allStat = cat(1,allStat,sesh(isesh).crossSpeed(ispeed1,ispeed2).statLocalTangling);
                allRun = cat(1,allRun,sesh(isesh).crossSpeed(ispeed1,ispeed2).runLocalTangling);
            end
        end

    end
end


figure, hold on
shadedErrorBar(-95:10:1695, mean(allStat,1), sem(allStat,1),'lineProps','k')
shadedErrorBar(-95:10:1695, mean(allRun,1), sem(allRun,1),'lineProps','r')
xlim([-200 1800])
ax=gca; ax.XTick = -200:200:1800;
defaultAxesProperties(gca, true)
ylabel('Between-condition local tangling')
xlabel('Time (ms)')


%% stats for cross-tangling

statVals = allStat';
runVals = allRun';

subjVec = categorical(repelem(1:5,1,30*180)); subjVec = cat(1,subjVec(:),subjVec(:));
timeVec = categorical(repmat(1:180,1,300)');
allVals = cat(1,statVals(:),runVals(:));
stateVec = categorical(repelem(1:2,1,numel(statVals))');

[p,tbl,stats,terms] = anovan(allVals,{timeVec,stateVec,subjVec},'model',2,'random',3,'varnames',{'Time','State','Subj'});




%% example of high and low tanglig

isession_stat = 2;
ispeed_stat=1; %3
isession_run = 2;
ispeed_run=1;
time_idx = 41:80;

cols = repelem(0,256,3);

statTraj = s(isession_stat).session.cond(ispeed_stat).meanTrajectory(:,1:3);
runTraj = s(isession_run).session.cond(isession_run+6).meanTrajectory(:,1:3);
deltaT = 10;
epsVal = 0.1;
Qstat=[];
Qrun=[];

itime = 1:200;
% get tangling of trajectories
for it = 2:numel(itime)
    for it2 = 2:numel(itime)

        xt = statTraj(itime(it),:);
        xt2 = statTraj(itime(it2),:);

        xderiv_t = (statTraj(itime(it),:) - statTraj(itime(it-1),:))/deltaT;
        xderiv_t2 = (statTraj(itime(it2),:) - statTraj(itime(it2-1),:))/deltaT;

        Qstat(it,it2) = (norm(xderiv_t-xderiv_t2)^2) / (norm(xt-xt2)^2 + epsVal);


        xt = runTraj(itime(it),:);
        xt2 = runTraj(itime(it2),:);

        xderiv_t = (runTraj(itime(it),:) - runTraj(itime(it-1),:))/deltaT;
        xderiv_t2 = (runTraj(itime(it2),:) - runTraj(itime(it2-1),:))/deltaT;

        Qrun(it,it2) = (norm(xderiv_t-xderiv_t2)^2) / (norm(xt-xt2)^2 + epsVal);



    end
end

for itime = 1:180
    tqs = this_Qstat(itime:itime+windowSize, itime:itime+windowSize);
    tqr = this_Qrun(itime:itime+windowSize, itime:itime+windowSize);


    tanglingVecStat(itime) = ...
        prctile(tqs(1+windowSize/2,:),90,'all');

    tanglingVecRun(itime) = ...
        prctile(tqr(1+windowSize/2,:),90,'all');

%         tanglingVecStat(itime) = ...
%         max(tqs(1+windowSize/2,:));
% 
%     tanglingVecRun(itime) = ...
%         max(tqr(1+windowSize/2,:));
end

tanglingVecStatOrig = tanglingVecStat(time_idx-10); % save real values
tanglingVecRunOrig = tanglingVecRun(time_idx-10);

tanglingVecStat = tanglingVecStat(time_idx-10);
tanglingVecRun = tanglingVecRun(time_idx-10);



[~,binEdges] = discretize([tanglingVecStat,tanglingVecRun],256);

tanglingVecStat = discretize(tanglingVecStat,binEdges);
tanglingVecRun = discretize(tanglingVecRun,binEdges);

figure,
plot(tanglingVecStat), hold on
plot(tanglingVecRun)

figure,
ax(1) = subplot(121);, hold on
% plot start and finsh point
plot3(statTraj(time_idx(1),1),...
    statTraj(time_idx(1),2),...
    statTraj(time_idx(1),3),'^', 'Color',cols(tanglingVecStat(1),:))

plot3(statTraj(time_idx(end),1),...
    statTraj(time_idx(end),2),...
    statTraj(time_idx(end),3),'o','Color',cols(tanglingVecStat(end),:))

for itime = 1:numel(time_idx)-1

    plot3(s(isession_stat).session.cond(ispeed_stat).meanTrajectory(time_idx(itime:itime+1),1),...
        s(isession_stat).session.cond(ispeed_stat).meanTrajectory(time_idx(itime:itime+1),2),...
        s(isession_stat).session.cond(ispeed_stat).meanTrajectory(time_idx(itime:itime+1),3),...
        '.','LineStyle','-','Color', cols(tanglingVecStat(itime),:))

end
grid off
view(-170, 35)
title(mean(tanglingVecStatOrig))

% run
ax(2) = subplot(122); hold on
% plot start and finsh point
plot3(s(isession_run).session.cond(ispeed_run+6).meanTrajectory(time_idx(1),1),...
    s(isession_run).session.cond(ispeed_run+6).meanTrajectory(time_idx(1),2),...
    s(isession_run).session.cond(ispeed_run+6).meanTrajectory(time_idx(1),3),'^', 'Color',cols(tanglingVecRun(1),:))

plot3(s(isession_run).session.cond(ispeed_run+6).meanTrajectory(time_idx(end),1),...
    s(isession_run).session.cond(ispeed_run+6).meanTrajectory(time_idx(end),2),...
    s(isession_run).session.cond(ispeed_run+6).meanTrajectory(time_idx(end),3),'o','Color',cols(tanglingVecRun(end),:))

for itime = 1:numel(time_idx)-1

    plot3(s(isession_run).session.cond(ispeed_run+6).meanTrajectory(time_idx(itime:itime+1),1),...
        s(isession_run).session.cond(ispeed_run+6).meanTrajectory(time_idx(itime:itime+1),2),...
        s(isession_run).session.cond(ispeed_run+6).meanTrajectory(time_idx(itime:itime+1),3),...
        '.','LineStyle','-','Color', cols(tanglingVecRun(itime),:))


end
title(mean(tanglingVecRunOrig))
grid off
view(-170, 35)


