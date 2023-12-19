%% Temporal correlation of total, signal and noise correlatons

%% load dataa
sessionTags = {'M22027', '20220517';...
    'M22029', '20220607';...
    'M22031', '20220524';
    'M22032', '20220621';
    'M22033', '20220706'};

dataDir = 'C:\Users\edward.horrocks\Documents\Code\V1Dynamics\Data\basic_111022';


for isession = 1:size(sessionTags,1)
    isession
    fname = [sessionTags{isession,1},'_', sessionTags{isession,2},'_FAcorrs_dev.mat']; %signalNoiseCorrs _corrsWithLDA

    load(fullfile(dataDir,fname))


    session1(isession).stat = session.stat;
    session1(isession).run = session.run;
    session1(isession).units = session.gU;
    
end


session = session1;



    %% Temporal correlation of total, signal and noise correlatons

statCorrArray = nan(181,181,5);
runCorrArray = nan(181,181,5);
tic
for isession = 1:5
    isession
    parfor itime1 = 1:181
        for itime2 = 1:181
            vals1 = session(isession).stat.totalCorrArray(:,:,itime1); vals1 = vals1(:);
            vals2 = session(isession).stat.totalCorrArray(:,:,itime2); vals2 = vals2(:);

            X = removeNANrows([vals1, vals2]);

            statCorrArray(itime1,itime2,isession) = ...
                corr(X(:,1), X(:,2), 'type', 'Pearson');

            vals1 = session(isession).run.totalCorrArray(:,:,itime1); vals1 = vals1(:);
            vals2 = session(isession).run.totalCorrArray(:,:,itime2); vals2 = vals2(:);

            X = removeNANrows([vals1, vals2]);

            runCorrArray(itime1,itime2,isession) = ...
                corr(X(:,1), X(:,2), 'type', 'Pearson');
        end
    end
end
toc

% plot correlation matrix

figure
ax1 = subplot(131)
imagesc(mean(statCorrArray, 3)), colorbar, caxis([0.4 1]), axis xy
% set(gca, 'XTick', [1 11 34 96 111 136], 'YTick', [1 11 34 96 111 136])
ax1.Colormap = turbo;

ax2 = subplot(132)
imagesc(mean(runCorrArray, 3)), colorbar, caxis([0.4 1]), axis xy
% set(gca, 'XTick', [1 11 34 96 111 136], 'YTick', [1 11 34 96 111 136])
ax2.Colormap = turbo;

ax3 = subplot(133)
imagesc(mean(runCorrArray, 3)-mean(statCorrArray, 3)), colorbar, caxis([0 1]), axis xy
% set(gca, 'XTick', [1 11 34 96 111 136], 'YTick', [1 11 34 96 111 136])
ax3.Colormap = crameri('vik');
caxis([-0.2 0.2])

totalStatCorrArray = statCorrArray;
totalRunCorrArray = runCorrArray;



statCorrArray = nan(181,181,5);
runCorrArray = nan(181,181,5);
tic
for isession = 1:5
    isession
    parfor itime1 = 1:181
        for itime2 = 1:181
            vals1 = session(isession).stat.signalCorrArray(:,:,itime1); vals1 = vals1(:);
            vals2 = session(isession).stat.signalCorrArray(:,:,itime2); vals2 = vals2(:);

            X = removeNANrows([vals1, vals2]);

            statCorrArray(itime1,itime2,isession) = ...
                corr(X(:,1), X(:,2), 'type', 'Pearson');

            vals1 = session(isession).run.signalCorrArray(:,:,itime1); vals1 = vals1(:);
            vals2 = session(isession).run.signalCorrArray(:,:,itime2); vals2 = vals2(:);

            X = removeNANrows([vals1, vals2]);

            runCorrArray(itime1,itime2,isession) = ...
                corr(X(:,1), X(:,2), 'type', 'Pearson');
        end
    end
end
toc

signalStatCorrArray = statCorrArray;
signalRunCorrArray = runCorrArray;


% plot correlation matrix

figure
ax1 = subplot(131)
imagesc(mean(signalStatCorrArray, 3)), colorbar, caxis([0 1]), axis xy
% set(gca, 'XTick', [1 11 34 96 111 136], 'YTick', [1 11 34 96 111 136])
ax1.Colormap = turbo;

ax2 = subplot(132)
imagesc(mean(signalRunCorrArray, 3)), colorbar, caxis([0 1]), axis xy
% set(gca, 'XTick', [1 11 34 96 111 136], 'YTick', [1 11 34 96 111 136])
ax2.Colormap = turbo;

ax3 = subplot(133)
imagesc(mean(signalRunCorrArray, 3)-mean(signalStatCorrArray, 3)), colorbar, caxis([0 1]), axis xy
% set(gca, 'XTick', [1 11 34 96 111 136], 'YTick', [1 11 34 96 111 136])
ax3.Colormap = crameri('vik');
caxis([-0.4 0.4])


statCorrArray = nan(181,181,5);
runCorrArray = nan(181,181,5);
tic
for isession = 1:5
    isession
    parfor itime1 = 1:181
        for itime2 = 1:181
            vals1 = session(isession).stat.noiseCorrArray(:,:,itime1); vals1 = vals1(:);
            vals2 = session(isession).stat.noiseCorrArray(:,:,itime2); vals2 = vals2(:);

            X = removeNANrows([vals1, vals2]);

            statCorrArray(itime1,itime2,isession) = ...
                corr(X(:,1), X(:,2), 'type', 'Pearson');

            vals1 = session(isession).run.noiseCorrArray(:,:,itime1); vals1 = vals1(:);
            vals2 = session(isession).run.noiseCorrArray(:,:,itime2); vals2 = vals2(:);

            X = removeNANrows([vals1, vals2]);

            runCorrArray(itime1,itime2,isession) = ...
                corr(X(:,1), X(:,2), 'type', 'Pearson');
        end
    end
end
toc

noiseStatCorrArray = statCorrArray;
noiseRunCorrArray = runCorrArray;

% plot correlation matrix

figure
ax1 = subplot(131)
imagesc(mean(statCorrArray, 3)), colorbar, caxis([0.4 1]), axis xy
% set(gca, 'XTick', [1 11 34 96 111 136], 'YTick', [1 11 34 96 111 136])
ax1.Colormap = turbo;

ax2 = subplot(132)
imagesc(mean(runCorrArray, 3)), colorbar, caxis([0.4 1]), axis xy
% set(gca, 'XTick', [1 11 34 96 111 136], 'YTick', [1 11 34 96 111 136])
ax2.Colormap = turbo;

ax3 = subplot(133)
imagesc(mean(runCorrArray, 3)-mean(statCorrArray, 3)), colorbar, caxis([0 1]), axis xy
% set(gca, 'XTick', [1 11 34 96 111 136], 'YTick', [1 11 34 96 111 136])
ax3.Colormap = crameri('vik');
caxis([-0.2 0.2])



%% plot by session

figure
for isession = 1:5
    subplot(2,5,isession)
    imagesc(statCorrArray(:,:,isession)), colorbar, caxis([0 1])
    subplot(2,5,isession+5)
    imagesc(runCorrArray(:,:,isession)), colorbar, caxis([0 1])
end
colormap(turbo)


%% plot slices through timepoints

meanStatCorrArray = mean(signalstatCorrArray,[3 4]);
meanRunCorrArray = mean(signalRunCorrArray, [3 4]);

time_idx = [1 24 96 136];

figure
for itime = 1:numel(time_idx)
    subplot(2,2,itime), hold on
    plot(meanStatCorrArray(time_idx(itime),:), 'k')
    plot(meanRunCorrArray(time_idx(itime),:), 'r')
    ax = gca; ax.XTick = 0:20:200; ax.YTick = 0:0.2:1;
    defaultAxesProperties(gca, true)
end


%% plot diagonal elements of temporal correlation matrix (non-overlapping)

% signal
for isession = 1:5
    statSignalCorrDiag(isession,:) = diag(signalStatCorrArray(:,:,isession),20);
    runSignalCorrDiag(isession,:) = diag(signalRunCorrArray(:,:,isession),20);
end

figure
shadedErrorBar(0:0.01:1.6, mean(statSignalCorrDiag,1), sem(statSignalCorrDiag,1))
shadedErrorBar(0:0.01:1.6, mean(runSignalCorrDiag,1), sem(runSignalCorrDiag,1), 'lineProps', 'r')
ylabel('Temporal Correlation (non-overlapping windows)'), xlabel('Intersection of time windows')
ylim([0 1]), xlim([-0.1 1.7])
title('Signal')
defaultAxesProperties(gca, false)

% noise
for isession = 1:5
    statNoiseCorrDiag(isession,:) = diag(noiseStatCorrArray(:,:,isession),20);
    runNoiseCorrDiag(isession,:) = diag(noiseRunCorrArray(:,:,isession),20);
end

figure
shadedErrorBar(0:0.01:1.6, mean(statNoiseCorrDiag,1), sem(statNoiseCorrDiag,1))
shadedErrorBar(0:0.01:1.6, mean(runNoiseCorrDiag,1), sem(runNoiseCorrDiag,1), 'lineProps', 'r')
ylabel('Temporal Correlation (non-overlapping windows)'), xlabel('Intersection of time windows')
ylim([0.4 1]), xlim([-0.1 1.7])
title('Noise')
defaultAxesProperties(gca, false)


%% stats

bv = 0:10:1600;
idx = find(bv==0):find(bv==800);

time = idx;
noise_stat = statNoiseCorrDiag(:,time)';
noise_run = runNoiseCorrDiag(:,time)';
noiseVec = cat(1,noise_stat(:), noise_run(:));
stateVec = categorical(cat(1,repelem(1, numel(noise_stat(:)))', repelem(2,numel(noise_run(:)))'));
timeVec = cat(1,repmat(time, 1,size(noise_stat,2))',repmat(time, 1,size(noise_run,2))');
subjVec = cat(1,repelem(1:size(noise_stat,2),1, numel(time))', repelem(1:size(noise_run,2),1, numel(time))');

[p,tbl,stats,terms] = anovan(noiseVec,{timeVec,stateVec,subjVec},'model',2,'random',3,'varnames',{'Time','State','Subj'});


time = idx;
signal_stat = statSignalCorrDiag(:,time)';
signal_run = runSignalCorrDiag(:,time)';
signalVec = cat(1,signal_stat(:), signal_run(:));
stateVec = categorical(cat(1,repelem(1, numel(signal_stat(:)))', repelem(2,numel(signal_run(:)))'));
timeVec = cat(1,repmat(time, 1,size(signal_stat,2))',repmat(time, 1,size(signal_run,2))');
subjVec = cat(1,repelem(1:size(signal_stat,2),1, numel(time))', repelem(1:size(signal_run,2),1, numel(time))');

[p,tbl,stats,terms] = anovan(signalVec,{timeVec,stateVec,subjVec},'model',2,'random',3,'varnames',{'Time','State','Subj'});



%% How similar is correlation structure between stat and run?

figure


corrVal = nan(5,181);
for isession =1:5
    for iint = 1:181
        statTemp = session(isession).stat.totalCorrArray(:,:,iint);
        runTemp = session(isession).run.totalCorrArray(:,:,iint);

        statTemp = setUpperTri2NaN(statTemp);
        runTemp = setUpperTri2NaN(runTemp);
        X = removeNANrows([statTemp(:), runTemp(:)]);

        corrVal(isession,iint) = corr(X(:,1), X(:,2));
    end
end

subplot(131)
shadedErrorBar(-100:10:1700, mean(corrVal,1), sem(corrVal,1), 'lineProps', 'm')
title('Total')
xlim([-200 1800])
ylim([-0.1 0.7])
ax = gca; ax.XTick = -200:200:1800;
defaultAxesProperties(gca,false)

corrVal = nan(5,181);
for isession =1:5
    for iint = 1:181
        statTemp = session(isession).stat.signalCorrArray(:,:,iint);
        runTemp = session(isession).run.signalCorrArray(:,:,iint);

        statTemp = setUpperTri2NaN(statTemp);
        runTemp = setUpperTri2NaN(runTemp);
        X = removeNANrows([statTemp(:), runTemp(:)]);

        corrVal(isession,iint) = corr(X(:,1), X(:,2));
    end
end

subplot(132)
shadedErrorBar(-100:10:1700, mean(corrVal,1), sem(corrVal,1), 'lineProps', 'm')
title('Signal')
xlim([-200 1800])
ylim([-0.1 0.7])
ax = gca; ax.XTick = -200:200:1800;
defaultAxesProperties(gca,false)


corrVal = nan(5,181);
for isession =1:5
    for iint = 1:181
        statTemp = session(isession).stat.noiseCorrArray(:,:,iint);
        runTemp = session(isession).run.noiseCorrArray(:,:,iint);

        statTemp = setUpperTri2NaN(statTemp);
        runTemp = setUpperTri2NaN(runTemp);
        X = removeNANrows([statTemp(:), runTemp(:)]);

        corrVal(isession,iint) = corr(X(:,1), X(:,2));
    end
end

subplot(133)
shadedErrorBar(-100:10:1700, mean(corrVal,1), sem(corrVal,1), 'lineProps', 'm')
title('Noise')
xlim([-200 1800])
ylim([-0.1 0.7])
ax = gca; ax.XTick = -200:200:1800;
defaultAxesProperties(gca,false)


