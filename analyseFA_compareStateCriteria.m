%% compare factor analysis results using diff state criteria

fileNames = {'FA_normal';...
            'FA_strict';...
            'FA_cp'};

plotTitles = {'original'; 'strict'; 'changepoints'};

sessionTags = {'M22027', '20220517';...
    'M22029', '20220607';...
    'M22031', '20220524';
    'M22032', '20220621';
    'M22033', '20220706'};

dataDir = 'C:\Users\edward.horrocks\Documents\Code\V1Dynamics\Data\basic_111022';

%%


for ifile = 1:numel(fileNames)
    ifile
    clear s
    s = struct;
    for isession = 1:5
        tic
        fname = [sessionTags{isession,1},'_', sessionTags{isession,2},'_', fileNames{ifile}, '.mat'];

        load(fullfile(dataDir,fname))

        s(isession).session = session;
    end


    %% get distance metrics with no weights

    bin_N = 1; % bin every N elements
    binSize = 10; % bin size in ms
    binVector = -200:binSize*bin_N:(1800-binSize*bin_N);
    timeBinVector = round(binVector,2); % time associated with each index (1st edge of bin)
    nWins = 10; % number of consecutive windows required to reach steady state

    for isession = 1:5

        clear output_offset output_onset
        % set parameters
        weights = s(isession).session.s.propSharedVariance(1:s(isession).session.s.qOpt)'; % weights are frac. of shared variance explained by each component
        %weights = repelem(1, 1, s(isession).session.s.qOpt);
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

        statPerf  =cat(1, statPerf, s(isession).session.stat.meanPerf);
        runPerf  =cat(1, runPerf, s(isession).session.run.meanPerf);

    end


    crit(ifile).statDR = cell2mat(statDR);
    crit(ifile).runDR  = cell2mat(runDR);

    crit(ifile).statDD = cell2mat(statDD);
    crit(ifile).runDD = cell2mat(runDD);

    crit(ifile).statCD = cell2mat(statCD);
    crit(ifile).runCD = cell2mat(runCD);

    crit(ifile).statDRoff = cell2mat(statDRoff);
    crit(ifile).runDRoff  = cell2mat(runDRoff);

    crit(ifile).statDDoff = cell2mat(statDDoff);
    crit(ifile).runDDoff = cell2mat(runDDoff);

    crit(ifile).statCDoff = cell2mat(statCDoff);
    crit(ifile).runCDoff = cell2mat(runCDoff);

    crit(ifile).statPerf  = statPerf;
    crit(ifile).runPerf  =runPerf;



end



%% plot results
speedcols = inferno(6);

for ifile = 1:3
    subplot(1,3,ifile), hold on
    plot([0 5], [0 5], 'k')
    for ispeed = 1:6
        plot(log2(crit(ifile).statDR(:,ispeed)), log2(crit(ifile).runDR(:,ispeed)), 'o', 'MarkerFaceColor', speedcols(ispeed,:), 'MarkerEdgeColor', 'none')
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
    title(plotTitles{ifile})
    defaultAxesProperties(gca, true)
end



%% plot delta DR
figure
for ifile = 1:3
    crit(ifile).deltaDR = crit(ifile).statDR./crit(ifile).runDR;
    crit(ifile).deltaDRoff = crit(ifile).statDRoff./crit(ifile).runDRoff;
end

figure,
vals = {log2(crit(1).deltaDR(:)), log2(crit(2).deltaDR(:)), log2(crit(3).deltaDR(:))};
distributionPlot(vals,'histOpt',1,'globalNorm',3,'color',[.7 .7 .7],'showMM',6)
hold on
plot([0 4], [0 0],'k:')
ylabel('delta distance ratio (Stat / Run)')
xlabel('State criteria')
ax = gca; ax.XTick = 1:3; ax.XTickLabel = plotTitles;
ax.YTick = [0 1 2 3 4]; ax.YTickLabel = 2.^(ax.YTick);
defaultAxesProperties(gca, false)


%% plot decoding performance

for ifile = 1:3
    crit(ifile).deltaPerf = crit(ifile).runPerf-crit(ifile).statPerf;
end

    bmp = -175:10:1775; bmp(end) = [];

figure
for ifile = 1:3
    subplot(1,3,ifile), hold on

allStat = crit(ifile).statPerf;
allRun = crit(ifile).runPerf;

plot([-0.2 1.8], [1/6 1/6], 'k:')
shadedErrorBar(bmp, mean(allStat,1), sem(allStat,1))
shadedErrorBar(bmp, mean(allRun,1), sem(allRun,1),'lineProps', 'r')
ax = gca; ax.XTick = -200:200:1800;
xlim([-200 1800])
defaultAxesProperties(gca,false)

title(plotTitles{ifile})

end


figure, hold on
plot([-0.2 1.8], [0 0], 'k:')

for ifile =1:3
    shadedErrorBar(bmp, mean(crit(ifile).deltaPerf,1), sem(crit(ifile).deltaPerf,1))
end
ax = gca; ax.XTick = -200:200:1800;
xlim([-200 1800])
defaultAxesProperties(gca,false)


ylabel('Delta LDA perf (run - stat')