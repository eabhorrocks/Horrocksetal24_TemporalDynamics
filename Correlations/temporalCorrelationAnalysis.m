%% Temporal correlation of total, signal and noise correlatons

%% load dataa
sessionTags = {'M22027', '20220517';...
    'M22029', '20220607';...
    'M22031', '20220524';
    'M22032', '20220621';
    'M22033', '20220706'};

dataDir = 'C:\Users\edward.horrocks\Documents\Code\V1Dynamics\Data\basic_111022';
%dataDir = '/mnt/rds01/ibn-vision/USERS/Edd/Code/V1Dynamics2024/data';

for isession = 1:size(sessionTags,1)
    isession
    fname = [sessionTags{isession,1},'_', sessionTags{isession,2},'_corrs_longAxis.mat']; %signalNoiseCorrs _corrsWithLDA

    load(fullfile(dataDir,fname))


    session1(isession).stat = session.stat;
    session1(isession).run = session.run;
    session1(isession).units = session.gU;

end


session = session1;



%% Temporal correlation of total, signal and noise correlatons
nInts = size(session(1).stat.totalCorrArray,3);
statCorrArray = nan(nInts,nInts,5);
runCorrArray = nan(nInts,nInts,5);
tic
for isession = 1:5
    isession
    parfor itime1 = 1:nInts
        for itime2 = 1:nInts
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



statCorrArray = nan(nInts,nInts,5);
runCorrArray = nan(nInts,nInts,5);
tic
for isession = 1:5
    isession
    parfor itime1 = 1:nInts
        for itime2 = 1:nInts
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


statCorrArray = nan(nInts,nInts,5);
runCorrArray = nan(nInts,nInts,5);
tic
for isession = 1:5
    isession
    parfor itime1 = 1:nInts
        for itime2 = 1:nInts
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


%% plot diagonal elements of temporal correlation matrix (non-overlapping)

% signal
for isession = 1:5
    statSignalCorrDiag(isession,:) = diag(signalStatCorrArray(:,:,isession),20);
    runSignalCorrDiag(isession,:) = diag(signalRunCorrArray(:,:,isession),20);
end

figure
shadedErrorBar(-0.1:0.01:1.7, mean(statSignalCorrDiag(:,11:191),1), sem(statSignalCorrDiag(:,11:191),1))
shadedErrorBar(-0.1:0.01:1.7, mean(runSignalCorrDiag(:,11:191),1), sem(runSignalCorrDiag(:,11:191),1), 'lineProps', 'r')
ylabel('Temporal Correlation (non-overlapping windows)'), xlabel('Intersection of time windows')
ylim([0 1]), xlim([-0.2 1.8])
title('Signal')
defaultAxesProperties(gca, false)

% noise
for isession = 1:5
    statNoiseCorrDiag(isession,:) = diag(noiseStatCorrArray(:,:,isession),20);
    runNoiseCorrDiag(isession,:) = diag(noiseRunCorrArray(:,:,isession),20);
end

figure
shadedErrorBar(-0.1:0.01:1.7, mean(statNoiseCorrDiag(:,11:191),1), sem(statNoiseCorrDiag(:,11:191),1))
shadedErrorBar(-0.1:0.01:1.7, mean(runNoiseCorrDiag(:,11:191),1), sem(runNoiseCorrDiag(:,11:191),1), 'lineProps', 'r')
ylabel('Temporal Correlation (non-overlapping windows)'), xlabel('Intersection of time windows')
ylim([0.4 1]), xlim([-0.2 1.8])
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

%[p,tbl,stats,terms] = anovan(noiseVec,{timeVec,stateVec,subjVec},'model',2,'random',3,'varnames',{'Time','State','Subj'});


time = idx;
signal_stat = statSignalCorrDiag(:,time)';
signal_run = runSignalCorrDiag(:,time)';
signalVec = cat(1,signal_stat(:), signal_run(:));
stateVec = categorical(cat(1,repelem(1, numel(signal_stat(:)))', repelem(2,numel(signal_run(:)))'));
timeVec = cat(1,repmat(time, 1,size(signal_stat,2))',repmat(time, 1,size(signal_run,2))');
subjVec = cat(1,repelem(1:size(signal_stat,2),1, numel(time))', repelem(1:size(signal_run,2),1, numel(time))');

%[p,tbl,stats,terms] = anovan(signalVec,{timeVec,stateVec,subjVec},'model',2,'random',3,'varnames',{'Time','State','Subj'});



%% How similar is correlation structure between stat and run?

figure


corrVal = nan(5,nInts);
for isession =1:5
    for iint = 1:nInts
        statTemp = session(isession).stat.totalCorrArray(:,:,iint);
        runTemp = session(isession).run.totalCorrArray(:,:,iint);

        statTemp = setUpperTri2NaN(statTemp);
        runTemp = setUpperTri2NaN(runTemp);
        X = removeNANrows([statTemp(:), runTemp(:)]);

        corrVal(isession,iint) = corr(X(:,1), X(:,2));
    end
end

subplot(131)
shadedErrorBar(-300:10:1900, mean(corrVal,1), sem(corrVal,1), 'lineProps', 'm')
title('Total')
xlim([-200 1800])
ylim([-0.1 0.7])
ax = gca; ax.XTick = -200:200:1800;
defaultAxesProperties(gca,false)

corrVal = nan(5,nInts);
for isession =1:5
    for iint = 1:nInts
        statTemp = session(isession).stat.signalCorrArray(:,:,iint);
        runTemp = session(isession).run.signalCorrArray(:,:,iint);

        statTemp = setUpperTri2NaN(statTemp);
        runTemp = setUpperTri2NaN(runTemp);
        X = removeNANrows([statTemp(:), runTemp(:)]);

        corrVal(isession,iint) = corr(X(:,1), X(:,2));
    end
end

subplot(132)
shadedErrorBar(-300:10:1900, mean(corrVal,1), sem(corrVal,1), 'lineProps', 'm')
title('Signal')
xlim([-200 1800])
ylim([-0.1 0.7])
ax = gca; ax.XTick = -200:200:1800;
defaultAxesProperties(gca,false)


corrVal = nan(5,nInts);
for isession =1:5
    for iint = 1:nInts
        statTemp = session(isession).stat.noiseCorrArray(:,:,iint);
        runTemp = session(isession).run.noiseCorrArray(:,:,iint);

        statTemp = setUpperTri2NaN(statTemp);
        runTemp = setUpperTri2NaN(runTemp);
        X = removeNANrows([statTemp(:), runTemp(:)]);

        corrVal(isession,iint) = corr(X(:,1), X(:,2));
    end
end

subplot(133)
shadedErrorBar(-300:10:1900, mean(corrVal,1), sem(corrVal,1), 'lineProps', 'm')
title('Noise')
xlim([-200 1800])
ylim([-0.1 0.7])
ax = gca; ax.XTick = -200:200:1800;
defaultAxesProperties(gca,false)


%% random pop sampling to match firing rates

nInts = size(session(1).stat.totalCorrArray,3);


tic
for isession = 1:5
    statSigCorrArray=[];
    runSigCorrArray=[];
    statNoiseCorrArray=[];
    runNoiseCorrArray=[];
    allFR = [];


    statNoise_corrArrayFull = session(isession).stat.noiseCorrArray;
    runNoise_corrArrayFull = session(isession).run.noiseCorrArray;

    statSig_corrArrayFull = session(isession).stat.signalCorrArray;
    runSig_corrArrayFull = session(isession).run.signalCorrArray;

    for irep = 1:1000
        irep

        nUnits = 20;

        units2use = randperm(numel(session(isession).units),nUnits);
        allFR(irep,:) = mean(cat(1,session(isession).units(units2use).dmFR));



        statNoise_sampCorrArray = statNoise_corrArrayFull(units2use,units2use,:);
        runNoise_sampCorrArray = runNoise_corrArrayFull(units2use,units2use,:);

        statSig_sampCorrArray = statSig_corrArrayFull(units2use,units2use,:);
        runSig_sampCorrArray = runSig_corrArrayFull(units2use,units2use,:);


        parfor itime1 = 1:nInts
            for itime2 = 1:nInts

                %%%noise%%%
                vals1 = statNoise_sampCorrArray(:,:,itime1); vals1 = vals1(:);
                vals2 = statNoise_sampCorrArray(:,:,itime2); vals2 = vals2(:);

                X = removeNANrows([vals1, vals2]);

                statNoiseCorrArray(itime1,itime2,irep) = ...
                    corr(X(:,1), X(:,2), 'type', 'Pearson');

                vals1 = runNoise_sampCorrArray(:,:,itime1); vals1 = vals1(:);
                vals2 = runNoise_sampCorrArray(:,:,itime2); vals2 = vals2(:);

                X = removeNANrows([vals1, vals2]);

                runNoiseCorrArray(itime1,itime2,irep) = ...
                    corr(X(:,1), X(:,2), 'type', 'Pearson');



                %%%signal%%%
                vals1 = statSig_sampCorrArray(:,:,itime1); vals1 = vals1(:);
                vals2 = statSig_sampCorrArray(:,:,itime2); vals2 = vals2(:);

                X = removeNANrows([vals1, vals2]);

                statSigCorrArray(itime1,itime2,irep) = ...
                    corr(X(:,1), X(:,2), 'type', 'Pearson');

                vals1 = runSig_sampCorrArray(:,:,itime1); vals1 = vals1(:);
                vals2 = runSig_sampCorrArray(:,:,itime2); vals2 = vals2(:);

                X = removeNANrows([vals1, vals2]);

                runSigCorrArray(itime1,itime2,irep) = ...
                    corr(X(:,1), X(:,2), 'type', 'Pearson');

            end
        end

        statNoiseDiag(irep,:) = diag(statNoiseCorrArray(:,:,irep),20);
        runNoiseDiag(irep,:) = diag(runNoiseCorrArray(:,:,irep),20);
        statSigDiag(irep,:) = diag(statSigCorrArray(:,:,irep),20);
        runSigDiag(irep,:) = diag(runSigCorrArray(:,:,irep),20);


    end
    toc

    sesh(isession).statNoiseDiag = statNoiseDiag;
    sesh(isession).runNoiseDiag = runNoiseDiag;
    sesh(isession).statSigDiag = statSigDiag;
    sesh(isession).runSigDiag = runSigDiag;
    sesh(isession).allFR = allFR;

end

%% subsample to match FRs

% function
% inputs: x,y OR [matrix]
% binResolution (OR autocalculate based on histcounts/histogram)
nReps = 10;

for isession = 1:5
    sesh(isession).m_statNoiseDiag = [];
    sesh(isession).m_runNoiseDiag = [];
    sesh(isession).m_statSignalDiag = [];
    sesh(isession).m_runSignalDiag = [];

    x = sesh(isession).allFR(:,1);
    y = sesh(isession).allFR(:,2);

    idx = downsampToMatch(x,y,'resolution',0.5,'nReps',nReps);

    for irep = 1:nReps
        sesh(isession).m_statNoiseDiag = ...
            cat(1,sesh(isession).m_statNoiseDiag,sesh(isession).statNoiseDiag(idx{1,irep},:));

        sesh(isession).m_runNoiseDiag = ...
            cat(1,sesh(isession).m_runNoiseDiag,sesh(isession).runNoiseDiag(idx{2,irep},:));

        sesh(isession).m_statSignalDiag =...
            cat(1,sesh(isession).m_statSignalDiag,sesh(isession).statSigDiag(idx{1,irep},:));

        sesh(isession).m_runSignalDiag =...
            cat(1, sesh(isession).m_runSignalDiag,sesh(isession).runSigDiag(idx{2,irep},:));

    end
end


%% plot means all reps
all_m_statNoiseDiag = cat(1,sesh.m_statNoiseDiag);
all_m_runNoiseDiag = cat(1,sesh.m_runNoiseDiag);

all_m_statSignalDiag = cat(1,sesh.m_statSignalDiag);
all_m_runSignalDiag = cat(1,sesh.m_runSignalDiag);

figure
shadedErrorBar(1:201, mean(all_m_statNoiseDiag), sem(all_m_statNoiseDiag),'lineProps','k')
shadedErrorBar(1:201, mean(all_m_runNoiseDiag), sem(all_m_runNoiseDiag),'lineProps','r')

figure
shadedErrorBar(1:201, mean(all_m_statSignalDiag), sem(all_m_statSignalDiag),'lineProps','k')
shadedErrorBar(1:201, mean(all_m_runSignalDiag), sem(all_m_runSignalDiag),'lineProps','r')

%% plott means of subjects

statNoise=[];
runNoise=[];
statSignal=[];
runSignal=[];

for isesh = 1:5
    statNoise=cat(1,statNoise,mean(sesh(isesh).m_statNoiseDiag));
    runNoise=cat(1,runNoise,mean(sesh(isesh).m_runNoiseDiag));
    statSignal=cat(1,statSignal,mean(sesh(isesh).m_statSignalDiag));
    runSignal=cat(1,runSignal,mean(sesh(isesh).m_runSignalDiag));
end

figure
subplot(121)
shadedErrorBar(1:201, mean(statNoise,1), sem(statNoise,1),'lineProps','k')
shadedErrorBar(1:201, mean(runNoise,1), sem(runNoise,1),'lineProps','r')

subplot(122)
shadedErrorBar(1:201, mean(statSignal,1), sem(statSignal,1),'lineProps','k')
shadedErrorBar(1:201, mean(runSignal,1), sem(runSignal,1),'lineProps','r')


