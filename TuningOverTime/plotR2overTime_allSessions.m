%% load data

sessionTags = {'M22027', '20220517';
    'M22029', '20220607';
    'M22031', '20220524';
    'M22032', '20220621';
    'M22033', '20220706'};

for isession=1:size(sessionTags,1)

    fname = [sessionTags{isession,1},'_', sessionTags{isession,2},'_r2overTime_dec22.mat'];

    load(fname)

    for iunit = 1:numel(tempUnits)
        tempUnits(iunit).session = isession;
    end

    session(isession).units = tempUnits;
end


%% get good units only
t2 = vertcat(session.units);
tempUnits = t2([t2.isi_viol]<=0.1 & [t2.amplitude_cutoff]<=0.1 & [t2.amplitude]>=50);

%% get cell type
for iunit = 1:numel(tempUnits)
    [ccg, t] = CCG(tempUnits(iunit).spike_times, ones(size(tempUnits(iunit).spike_times)),...
        'binSize', 0.0005, 'duration', 0.1,'norm', 'rate');
    tempUnits(iunit).acg = ccg;
    fit_params_out = fit_ACG(ccg,false);

    tempUnits(iunit).tau_rise = fit_params_out.acg_tau_rise;
end


% Criteria:
% narrow waveform (trough-to-peak ≤ 450 μs),
% wide waveform (trough-to-peak > 450 μs and τrise > 6 ms), putative interneurons,
% and the rest, as pyramidal cells


narrow_idx = find([tempUnits.duration]<=0.45);
wide_idx = find([tempUnits.duration]>0.45 & [tempUnits.tau_rise]>6);
pyr_idx = find(~ismember(1:numel(tempUnits), [narrow_idx,wide_idx]));

narrowUnits = tempUnits(narrow_idx);
wideUnits = tempUnits(wide_idx);
pyrUnits = tempUnits(pyr_idx);


%% mcnemar test on p(tuned)

                     % run good
%                  +         -
%             --------------------
%         +   |   1,1   |   1,2   |
% stat good   |-------------------           
%         -   |   2,1   |   2,2   |
%             --------------------
%                                       
%



statGoodVals = ~isnan(vertcat(tempUnits.stat_tunedStart));
runGoodVals = ~isnan(vertcat(tempUnits.run_tunedStart));


x(1,1) = sum(statGoodVals & runGoodVals);
x(1,2) = sum(statGoodVals & ~runGoodVals);
x(2,1) = sum(~statGoodVals & runGoodVals);
x(2,2) = sum(~statGoodVals & ~runGoodVals);

x
 p = mcnemar(x)


%% plot and stats for tuning start times

allTunedStart = [vertcat(tempUnits.stat_tunedStart), vertcat(tempUnits.run_tunedStart)];
allTunedStart(allTunedStart<0)=0;
allTunedStart_valid = allTunedStart(all(~isnan(allTunedStart),2),:);

sessionVec = vertcat(tempUnits.session);
sessionVec = sessionVec(all(~isnan(allTunedStart),2));

[statCI_lower, statCI_upper] = medianCI(allTunedStart_valid(:,1), 1.96);
[runCI_lower, runCI_upper] = medianCI(allTunedStart_valid(:,2), 1.96);

figure, hold on
bar([1 2], median(allTunedStart_valid))
plot([1 1], [statCI_lower, statCI_upper], 'k')
plot([2 2], [runCI_lower, runCI_upper], 'k')

figure, hold on
scatter_kde(allTunedStart_valid(:,1), allTunedStart_valid(:,2), 'filled', 'MarkerSize', 10)
xlim([0 1.8]); ylim([0 1.8])
ax = gca; ax.XTick = 0:0.2:1.8; ax.YTick = 0:0.2:1.8;
plot([-0.2 1.8], [-0.2 1.8], 'k-')
axis equal
defaultAxesProperties(gca, true)
colormap(viridis)

valsVec = cat(1,allTunedStart_valid(:,1), allTunedStart_valid(:,2));
unitVec = categorical(cat(2,1:size(allTunedStart_valid,1), 1:size(allTunedStart_valid,1))');
stateVec = categorical(cat(1,repelem(1,size(allTunedStart_valid,1),1), repelem(2,size(allTunedStart_valid,1),1)));
seshVec = categorical(cat(1,sessionVec, sessionVec));

tbl = table(valsVec,unitVec,stateVec,seshVec,...
    'VariableNames',{'vals','unit','state','sesh'});

    f = 'vals ~ state + (1|unit) + (1|sesh)';
      
    lme = fitlme(tbl, f, 'DummyVarCoding', 'reference')

%% plot and stats for tuning finish times

allTunedFinish = [vertcat(tempUnits.stat_tunedFinish), vertcat(tempUnits.run_tunedFinish)];
allTunedFinish_valid = allTunedFinish(all(~isnan(allTunedFinish),2),:);

figure, hold on
scatter_kde(allTunedFinish_valid(:,1), allTunedFinish_valid(:,2), 'filled', 'MarkerSize', 10)
xlim([0 1.8]); ylim([0 1.8])
ax = gca; ax.XTick = 0:0.2:1.8; ax.YTick = 0:0.2:1.8;
plot([-0.2 1.8], [-0.2 1.8], 'k-')
axis equal
defaultAxesProperties(gca, true)
colormap(viridis)

figure, hold on
[statCI_lower, statCI_upper] = medianCI(allTunedFinish_valid(:,1), 1.96);
[runCI_lower, runCI_upper] = medianCI(allTunedFinish_valid(:,2), 1.96);
bar([1 2], median(allTunedFinish_valid))
plot([1 1], [statCI_lower, statCI_upper], 'k')
plot([2 2], [runCI_lower, runCI_upper], 'k')


valsVec = cat(1,allTunedFinish_valid(:,1), allTunedFinish_valid(:,2));
unitVec = categorical(cat(2,1:size(allTunedStart_valid,1), 1:size(allTunedStart_valid,1))');
stateVec = categorical(cat(1,repelem(1,size(allTunedStart_valid,1),1), repelem(2,size(allTunedStart_valid,1),1)));
seshVec = categorical(cat(1,sessionVec, sessionVec));

tbl = table(valsVec,unitVec,stateVec,seshVec,...
    'VariableNames',{'vals','unit','state','sesh'});

    f = 'vals ~ state + (1|unit) + (1|sesh)';
      
    lme = fitlme(tbl, f, 'DummyVarCoding', 'reference')

%% plot and stats for tuning duration times

allTunedDuration = [vertcat(tempUnits.stat_tunedDuration), vertcat(tempUnits.run_tunedDuration)];
allTunedDuration_valid = allTunedDuration(all(~isnan(allTunedDuration),2),:);

figure, hold on
[statCI_lower, statCI_upper] = medianCI(allTunedDuration_valid(:,1), 1.96);
[runCI_lower, runCI_upper] = medianCI(allTunedDuration_valid(:,2), 1.96);
bar([1 2], median(allTunedDuration_valid))
plot([1 1], [statCI_lower, statCI_upper], 'k')
plot([2 2], [runCI_lower, runCI_upper], 'k')


figure, hold on
scatter_kde(allTunedDuration_valid(:,1)./100, allTunedDuration_valid(:,2)./100, 'filled', 'MarkerSize', 10)
xlim([0 1.8]); ylim([0 1.8])
ax = gca; ax.XTick = 0:0.2:1.8; ax.YTick = 0:0.2:1.8;
plot([-0.2 1.8], [-0.2 1.8], 'k-')
axis equal
defaultAxesProperties(gca, true)
colormap(viridis)

valsVec = cat(1,allTunedDuration_valid(:,1), allTunedDuration_valid(:,2));
unitVec = categorical(cat(2,1:size(allTunedStart_valid,1), 1:size(allTunedStart_valid,1))');
stateVec = categorical(cat(1,repelem(1,size(allTunedStart_valid,1),1), repelem(2,size(allTunedStart_valid,1),1)));
seshVec = categorical(cat(1,sessionVec, sessionVec));

tbl = table(valsVec,unitVec,stateVec,seshVec,...
    'VariableNames',{'vals','unit','state','sesh'});

    f = 'vals ~ state + (1|unit) + (1|sesh)';
      
    lme = fitlme(tbl, f, 'DummyVarCoding', 'reference')


%% proportion tuned over time

nUnits = numel(tempUnits);

figure
subplot(131), hold on
plot(1:181, sum(vertcat(tempUnits.stat_tuneFlagWeak))./nUnits, 'k')
plot(1:181, sum(vertcat(tempUnits.run_tuneFlagWeak))./nUnits, 'r')
ax = gca; ax.XTick = 1:20:181;
ylim([0 0.8])
defaultAxesProperties(gca, true)

subplot(132), hold on
plot(1:181, sum(vertcat(tempUnits.stat_tuneFlagMed))./nUnits, 'k')
plot(1:181, sum(vertcat(tempUnits.run_tuneFlagMed))./nUnits, 'r')
ylim([0 0.8])
ax = gca; ax.XTick = 1:20:181;
defaultAxesProperties(gca, true)


subplot(133), hold on
plot(1:181, sum(vertcat(tempUnits.stat_tuneFlagStrong))./nUnits, 'k')
plot(1:181, sum(vertcat(tempUnits.run_tuneFlagStrong))./nUnits, 'r')
ylim([0 0.8])
ax = gca; ax.XTick = 1:20:181;
defaultAxesProperties(gca, true)

%% histogram versions of plots (like biorxiv manuscript).

allTunedStart = [vertcat(tempUnits.stat_tunedStart), vertcat(tempUnits.run_tunedStart)];
allTunedStart_valid = allTunedStart(all(~isnan(allTunedStart),2),:);

figure, subplot(211), hold on
allTunedStart_valid = allTunedStart(all(~isnan(allTunedStart),2),:);

statVals = allTunedStart_valid(:,1);
runVals = allTunedStart_valid(:,2);

subplot(211)
h_stat = histogram([statVals], 'BinEdges', [-0.1:0.05:1.75],'normalization', 'probability', 'FaceColor', [0 0 0]);
h_run = histogram([runVals], 'BinEdges', [-0.1:0.05:1.75],'normalization', 'probability', 'FaceColor', [.7 .7 .7]);


statQuan = quantile([statVals], 3);
runQuan = quantile([runVals], 3);

plot([statQuan(1), statQuan(3)], [0.35 0.35], 'Color', [0 0 0])
plot(statQuan(2), 0.35, 'o', 'Color', [0 0 0])
plot([runQuan(1), runQuan(3)], [0.38 0.38], '-', 'Color', [.7 .7 .7])
plot(runQuan(2), 0.38, 'o','Color', [.7 .7 .7])
ax = gca; ax.XTick = -0.2:0.2:1.8;
xlim([-0.2 1.8]); ylim([0 0.4])
title('Tuning Start')
defaultAxesProperties(gca, true)

figure, hold on
plot([-0.2, 1.8], [0 0], 'k-')
sh = stairs([-0.1:0.05:1.75], [h_run.Values-h_stat.Values, h_run.Values(end)-h_stat.Values(end)]);
xlim([-0.2 1.8]); ax = gca; ax.XTick = -0.2:0.2:1.8;
ylim([-0.07 0.07]), ax.YTick = [-0.07,0, 0.07];
defaultAxesProperties(gca, true)

for iseg = 1:numel(sh.XData)-1
    if sh.YData(iseg)>=0
        col2use = 'r';
    else
        col2use = 'k';
    end
    fill([sh.XData(iseg), sh.XData(iseg), sh.XData(iseg+1), sh.XData(iseg+1)],...
        [0, sh.YData(iseg), sh.YData(iseg), 0], col2use)
end


%% tuning finish
figure, subplot(211), hold on

statVals = allTunedFinish_valid(:,1);
runVals = allTunedFinish_valid(:,2);
statVals(statVals<0)=0;
runVals(runVals<0)=0;

subplot(211)
h_stat = histogram([statVals], 'BinEdges', [-0.1:0.05:1.75],'normalization', 'probability', 'FaceColor', [0 0 0]);
h_run = histogram([runVals], 'BinEdges', [-0.1:0.05:1.75],'normalization', 'probability', 'FaceColor', [.7 .7 .7]);


statQuan = quantile([statVals], 3);
runQuan = quantile([runVals], 3);

plot([statQuan(1), statQuan(3)], [0.2 0.2], 'Color', [0 0 0])
plot(statQuan(2), 0.2, 'o', 'Color', [0 0 0])
plot([runQuan(1), runQuan(3)], [0.22 0.22], '-', 'Color', [.7 .7 .7])
plot(runQuan(2), 0.22, 'o','Color', [.7 .7 .7])
ax = gca; ax.XTick = -0.2:0.2:1.8;
xlim([-0.2 1.8]); ylim([0 0.25])
title('Tuning Finish')
defaultAxesProperties(gca, true)

figure, hold on
plot([-0.2, 1.8], [0 0], 'k-')
sh = stairs([-0.1:0.05:1.75], [h_run.Values-h_stat.Values, h_run.Values(end)-h_stat.Values(end)]);
xlim([-0.2 1.8]); ax = gca; ax.XTick = -0.2:0.2:1.8;
ylim([-0.05 0.05]), ax.YTick = [-0.05,0, 0.05];
defaultAxesProperties(gca, true)

for iseg = 1:numel(sh.XData)-1
    if sh.YData(iseg)>=0
        col2use = 'r';
    else
        col2use = 'k';
    end
    fill([sh.XData(iseg), sh.XData(iseg), sh.XData(iseg+1), sh.XData(iseg+1)],...
        [0, sh.YData(iseg), sh.YData(iseg), 0], col2use)
end


%% tuning duration
figure, subplot(211), hold on

statVals = allTunedDuration_valid(:,1)./100;
runVals = allTunedDuration_valid(:,2)./100;

subplot(211)
h_stat = histogram([statVals], 'BinEdges', [-0.1:0.05:1.75],'normalization', 'probability', 'FaceColor', [0 0 0]);
h_run = histogram([runVals], 'BinEdges', [-0.1:0.05:1.75],'normalization', 'probability', 'FaceColor', [.7 .7 .7]);


statQuan = quantile([statVals], 3);
runQuan = quantile([runVals], 3);

plot([statQuan(1), statQuan(3)], [0.12 0.12], 'Color', [0 0 0])
plot(statQuan(2), 0.12, 'o', 'Color', [0 0 0])
plot([runQuan(1), runQuan(3)], [0.14 0.14], '-', 'Color', [.7 .7 .7])
plot(runQuan(2), 0.14, 'o','Color', [.7 .7 .7])
ax = gca; ax.XTick = -0.2:0.2:1.8;
xlim([-0.2 1.8]); ylim([0 0.15])
title('Tuning Duration')
defaultAxesProperties(gca, true)

figure, hold on
plot([-0.2, 1.8], [0 0], 'k-')
sh = stairs([-0.1:0.05:1.75], [h_run.Values-h_stat.Values, h_run.Values(end)-h_stat.Values(end)]);
xlim([-0.2 1.8]); ax = gca; ax.XTick = -0.2:0.2:1.8;
ylim([-0.07 0.07]), ax.YTick = [-0.07,0, 0.07];
defaultAxesProperties(gca, true)

for iseg = 1:numel(sh.XData)-1
    if sh.YData(iseg)>=0
        col2use = 'r';
    else
        col2use = 'k';
    end
    fill([sh.XData(iseg), sh.XData(iseg), sh.XData(iseg+1), sh.XData(iseg+1)],...
        [0, sh.YData(iseg), sh.YData(iseg), 0], col2use)
end





%% dynamic range



dr_stat = [];
dr_run = [];


for iunit = 1:numel(tempUnits)
    dr_stat = [dr_stat; [tempUnits(iunit).stat_ints.dynamicRange]];
    %pdmean_stat = [pdmean_stat; [tempUnits(iunit).stat_ints.meanPairwiseDiffs]];

    dr_run = [dr_run; [tempUnits(iunit).run_ints.dynamicRange]];
    %pdmean_run = [pdmean_run; [tempUnits(iunit).run_ints.meanPairwiseDiffs]];

end

% valid_idx = all(~isnan(ffmean_stat),2) & all(~isnan(ffmean_run),2)

figure, hold on
shadedErrorBar(10:190, nanmean(dr_stat,1), nansem(dr_stat,1), 'lineProps', {'Color', 'k'})
shadedErrorBar(10:190, nanmean(dr_run,1), nansem(dr_run,1), 'lineProps', {'Color', 'r'})


%% dynamic range stats

dr_stat_mean = mean(dr_stat(:,11:111),2);
dr_run_mean = mean(dr_run(:,11:111),2);

sessionVec = cat(1,tempUnits.session);

valsVec = cat(1,dr_stat_mean, dr_run_mean);
unitVec = categorical(cat(2,1:numel(dr_stat_mean), 1:numel(dr_stat_mean))');
stateVec = categorical(cat(1,repelem(1,numel(dr_stat_mean),1), repelem(2,numel(dr_stat_mean),1)));
seshVec = categorical(cat(1,sessionVec, sessionVec));

tbl = table(valsVec,unitVec,stateVec,seshVec,...
    'VariableNames',{'vals','unit','state','sesh'});

    f = 'vals ~ state + (1|unit) + (1|sesh)';
      
    lme = fitlme(tbl, f, 'DummyVarCoding', 'reference')



%% Cell-type tuning metrics

narrowUnits = tempUnits(narrow_idx);
wideUnits = tempUnits(wide_idx);
pyrUnits = tempUnits(pyr_idx);

typeNames = {'pyr', 'narrow', 'wide'};

%% tuning start
for itype = 1:3
    figure
    switch itype
        case 1
            allTunedStart = [vertcat(pyrUnits.stat_tunedStart), vertcat(pyrUnits.run_tunedStart)];
        case 2
            allTunedStart = [vertcat(narrowUnits.stat_tunedStart), vertcat(narrowUnits.run_tunedStart)];
        case 3
            allTunedStart = [vertcat(wideUnits.stat_tunedStart), vertcat(wideUnits.run_tunedStart)];
    end

    allTunedStart(allTunedStart<0)=0;
    allTunedStart_valid = allTunedStart(all(~isnan(allTunedStart),2),:);

    %
    [statCI_lower, statCI_upper] = medianCI(allTunedStart_valid(:,1), 1.96);
    [runCI_lower, runCI_upper] = medianCI(allTunedStart_valid(:,2), 1.96);

    subplot(122), hold on
    bar([1 2], median(allTunedStart_valid))
    plot([1 1], [statCI_lower, statCI_upper], 'k')
    plot([2 2], [runCI_lower, runCI_upper], 'k')
    p = signrank(allTunedStart_valid(:,1), allTunedStart_valid(:,2));
    ylim([0 0.3])
    title(num2str(p,3))

    subplot(121), hold on
    scatter_kde(allTunedStart_valid(:,1), allTunedStart_valid(:,2), 'filled', 'MarkerSize', 10)
    xlim([0 1.8]); ylim([0 1.8])
    ax = gca; ax.XTick = 0:0.2:1.8; ax.YTick = 0:0.2:1.8;
    plot([-0.2 1.8], [-0.2 1.8], 'k-')
    % axis equal
    title(typeNames{itype})
    defaultAxesProperties(gca, true)
    colormap(viridis)
end


%% tuning finish
for itype = 1:3
    figure
    switch itype
        case 1
            allTunedStart = [vertcat(pyrUnits.stat_tunedFinish), vertcat(pyrUnits.run_tunedFinish)];
        case 2
            allTunedStart = [vertcat(narrowUnits.stat_tunedFinish), vertcat(narrowUnits.run_tunedFinish)];
        case 3
            allTunedStart = [vertcat(wideUnits.stat_tunedFinish), vertcat(wideUnits.run_tunedFinish)];
    end

    allTunedStart(allTunedStart<0)=0;
    allTunedStart_valid = allTunedStart(all(~isnan(allTunedStart),2),:);

    %
    [statCI_lower, statCI_upper] = medianCI(allTunedStart_valid(:,1), 1.96);
    [runCI_lower, runCI_upper] = medianCI(allTunedStart_valid(:,2), 1.96);

    subplot(122), hold on
    bar([1 2], median(allTunedStart_valid))
    plot([1 1], [statCI_lower, statCI_upper], 'k')
    plot([2 2], [runCI_lower, runCI_upper], 'k')
    p = signrank(allTunedStart_valid(:,1), allTunedStart_valid(:,2));
    ylim([0 1.8])
    title(num2str(p,3))

    subplot(121), hold on
    scatter_kde(allTunedStart_valid(:,1), allTunedStart_valid(:,2), 'filled', 'MarkerSize', 10)
    xlim([0 1.8]); ylim([0 1.8])
    ax = gca; ax.XTick = 0:0.2:1.8; ax.YTick = 0:0.2:1.8;
    plot([-0.2 1.8], [-0.2 1.8], 'k-')
    % axis equal
    title(typeNames{itype})
    defaultAxesProperties(gca, true)
    colormap(viridis)
end

%% tuning duration
for itype = 1:3
    figure
    switch itype
        case 1
            allTunedStart = [vertcat(pyrUnits.stat_tunedDuration), vertcat(pyrUnits.run_tunedDuration)]./100;
        case 2
            allTunedStart = [vertcat(narrowUnits.stat_tunedDuration), vertcat(narrowUnits.run_tunedDuration)]./100;
        case 3
            allTunedStart = [vertcat(wideUnits.stat_tunedDuration), vertcat(wideUnits.run_tunedDuration)]./100;
    end

    allTunedStart(allTunedStart<0)=0;
    allTunedStart_valid = allTunedStart(all(~isnan(allTunedStart),2),:);

    %
    [statCI_lower, statCI_upper] = medianCI(allTunedStart_valid(:,1), 1.96);
    [runCI_lower, runCI_upper] = medianCI(allTunedStart_valid(:,2), 1.96);

    subplot(122), hold on
    bar([1 2], median(allTunedStart_valid))
    plot([1 1], [statCI_lower, statCI_upper], 'k')
    plot([2 2], [runCI_lower, runCI_upper], 'k')
    p = signrank(allTunedStart_valid(:,1), allTunedStart_valid(:,2));
    ylim([0 1.8])
    title(num2str(p,3))

    subplot(121), hold on
    scatter_kde(allTunedStart_valid(:,1), allTunedStart_valid(:,2), 'filled', 'MarkerSize', 10)
    xlim([0 1.8]); ylim([0 1.8])
    ax = gca; ax.XTick = 0:0.2:1.8; ax.YTick = 0:0.2:1.8;
    plot([-0.2 1.8], [-0.2 1.8], 'k-')
    % axis equal
        title(typeNames{itype})

    defaultAxesProperties(gca, true)
    colormap(viridis)
end


%% tuned over time

for itype = 1:3

    switch itype
        case 1
            theseUnits = pyrUnits;
        case 2
            theseUnits = narrowUnits;
        case 3
            theseUnits = wideUnits;
    end

    nUnits = numel(theseUnits);

    figure
%     subplot(131), 
    hold on
    plot(-100:10:1700, sum(vertcat(theseUnits.stat_tuneFlagWeak))./nUnits, 'k')
    plot(-100:10:1700, sum(vertcat(theseUnits.run_tuneFlagWeak))./nUnits, 'r')
    ax = gca; ax.XTick = -200:200:1800;
        title(typeNames{itype})
% 
    ylim([0 0.5])
    defaultAxesProperties(gca, true)
% 
%     subplot(132), hold on
%     plot(1:181, sum(vertcat(theseUnits.stat_tuneFlagMed))./nUnits, 'k')
%     plot(1:181, sum(vertcat(theseUnits.run_tuneFlagMed))./nUnits, 'r')
%     ylim([0 0.5])
%         title(typeNames{itype})
% 
%     ax = gca; ax.XTick = 1:20:181;
%     defaultAxesProperties(gca, true)
% 
% 
%     subplot(133), hold on
%     plot(1:181, sum(vertcat(theseUnits.stat_tuneFlagStrong))./nUnits, 'k')
%     plot(1:181, sum(vertcat(theseUnits.run_tuneFlagStrong))./nUnits, 'r')
%     ylim([0 0.5])
%         title(typeNames{itype})
% 
%     ax = gca; ax.XTick = 1:20:181;
%     defaultAxesProperties(gca, true)
end


%%

for itype = 1:3

    switch itype
        case 1
            theseUnits = pyrUnits;
        case 2
            theseUnits = narrowUnits;
        case 3
            theseUnits = wideUnits;
    end

    theseUnits = theseUnits([theseUnits.stat_nBouts]>=0 | [theseUnits.run_nBouts]>=0);




    ffmean_stat = [];
    dr_stat = [];
    pdmean_stat = [];

    ffmean_run = [];
    dr_run = [];
    pdmean_run = [];


    for iunit = 1:numel(theseUnits)
        ffmean_stat = [ffmean_stat; [theseUnits(iunit).stat_ints.meanFanoFactor]];
        dr_stat = [dr_stat; [theseUnits(iunit).stat_ints.dynamicRange]];
        %pdmean_stat = [pdmean_stat; [tempUnits(iunit).stat_ints.meanPairwiseDiffs]];

        ffmean_run = [ffmean_run; [theseUnits(iunit).run_ints.meanFanoFactor]];
        dr_run = [dr_run; [theseUnits(iunit).run_ints.dynamicRange]];
        %pdmean_run = [pdmean_run; [tempUnits(iunit).run_ints.meanPairwiseDiffs]];

    end

     valid_idx = all(~isnan(ffmean_stat),2) & all(~isnan(ffmean_run),2)

    figure,
    subplot(211), hold on
    shadedErrorBar(-100:10:1700, nanmean(dr_stat(valid_idx,:),1), nansem(dr_stat(valid_idx,:),1), 'lineProps', {'Color', 'k'})
    shadedErrorBar(-100:10:1700, nanmean(dr_run(valid_idx,:),1), nansem(dr_run(valid_idx,:),1), 'lineProps', {'Color', 'r'})
    ylim([0 20])
     ax = gca; ax.XTick = -200:200:1800;
        title(typeNames{itype})
    defaultAxesProperties(gca, true)

    subplot(212), hold on
    shadedErrorBar(-100:10:1700, nanmean(ffmean_stat,1), nansem(ffmean_stat,1), 'lineProps', {'Color', 'k'})
    shadedErrorBar(-100:10:1700, nanmean(ffmean_run,1), nansem(ffmean_run,1), 'lineProps', {'Color', 'r'})
    ylim([0.9, 2])
         ax = gca; ax.XTick = -200:200:1800;
        title(typeNames{itype})
    defaultAxesProperties(gca, true)

end

%% plot example tuning strength neurons


%
%     good examples (m22029) [tempUnits idx]
%     239, excellent (onset + offset) fr range: [0 26]
%     300 is ok, [5 75ish]
%     226 for onset [0 28]
%     199 pretty good for both

idx = 239 ;%= find([tempUnits.stat_tunedStart]>0.15 & [tempUnits.run_tunedStart]<0.1);
%     229 is good for onset

binWidth = 0.01;
options.smoothType = 'gaussian';
options.smoothWidth = 175;
options.plot = true;
options.cols = inferno(6);
options.preTime = 200;
options.postTime = 800;
options.plotFRrange = [];


%
for iunit = 1:numel(idx)
    figure,
    subplot(311), hold on
    stimes = tempUnits(idx(iunit)).spike_times'.*1000;
    [pout, info] = spikeTimesPSTH(stimes,stat_intervalsCell,binWidth,options);
    title(idx(iunit))
    subplot(312)
    [pout, info] = spikeTimesPSTH(stimes,run_intervalsCell,binWidth,options);
    title(idx(iunit))

    subplot(313), hold on
    plot(0:10:1800, [tempUnits(idx(iunit)).stat_ints.R2], 'k')
    plot(0:10:1800, [tempUnits(idx(iunit)).run_ints.R2], 'r')
    plot([-200 1800], [0.1 0.1], 'k')
    title(num2str(idx(iunit)))
    xlim([-200 1800])
    ylim([-1 1])
end




%% plot tuning curve with tuning strength values for exmaple neuron
% subj=M22029, tempUnits, #239

binVector = -200:10:1790;
r2binVector = -100:10:1700;

thisUnit = tempUnits(239);


figure,


time_idx = find(binVector==200):find(binVector==400);
title_idx = find(r2binVector==300);
theseSpikes = cellfun(@(x) mean(x(time_idx,:)), thisUnit.allSpikes, 'UniformOutput', false);
statSpikes = cell2mat(theseSpikes(:,1));
runSpikes = cell2mat(theseSpikes(:,2));

nTrials = size(statSpikes,2);

subplot(221), hold on
errorbar(1:6, mean(statSpikes,2)*100, sem(statSpikes,2)*100, 'k')
ylim([0 10])
xlim([0.5 6.5])
title(thisUnit.stat_ints(title_idx).R2)
defaultAxesProperties(gca,true)


subplot(223)
errorbar(1:6, mean(runSpikes,2)*100, sem(runSpikes,2)*100, 'r')
title(thisUnit.run_ints(title_idx).R2)
xlim([0.5 6.5])
ylim([0 25])
defaultAxesProperties(gca,true)


time_idx = find(binVector==800):find(binVector==990);
title_idx = find(r2binVector==900);

theseSpikes = cellfun(@(x) mean(x(time_idx,:)), thisUnit.allSpikes, 'UniformOutput', false);
statSpikes = cell2mat(theseSpikes(:,1));
runSpikes = cell2mat(theseSpikes(:,2));

nTrials = size(statSpikes,2);

subplot(222), hold on
errorbar(1:6, mean(statSpikes,2)*100, sem(statSpikes,2)*100, 'k')
xlim([0.5 6.5])
title(thisUnit.stat_ints(title_idx).R2)
ylim([0 10])
defaultAxesProperties(gca,true)

subplot(224)
errorbar(1:6, mean(runSpikes,2)*100, sem(runSpikes,2)*100, 'r')
xlim([0.5 6.5])
title(thisUnit.run_ints(title_idx).R2)
ylim([0 25])
defaultAxesProperties(gca,true)