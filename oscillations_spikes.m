%% oscillations in spike trains?

%% spectrum of psths %%

%% load data
load('PSTHdata_wCellTypes.mat')

%% get reliable responses
rangeThresh = 3;
reliThrsh =-1.645;
zThresh = 3.29;

statGoodVals =stat.reliVec<=reliThrsh & stat.rangeFRVec>=rangeThresh & stat.peakAbsZVec>=zThresh;
runGoodVals = run.reliVec<=reliThrsh & run.rangeFRVec>=rangeThresh & run.peakAbsZVec>=zThresh;
bothGoodVals = statGoodVals & runGoodVals;


allStat = stat.psthVec(statGoodVals,:);
allRun = run.psthVec(runGoodVals,:);


%% matlab
[pstat, fstat] = pspectrum(allStat',100,'FrequencyLimits',[1 10]);
[prun, frun] = pspectrum(allRun',100,'FrequencyLimits',[1 10]);


figure, hold on
shadedErrorBar(fstat,mean(pow2db(pstat/mean(pstat)),2), sem(pow2db(pstat/mean(pstat)),2),'lineProps','k')
shadedErrorBar(frun,mean(pow2db(prun/mean(prun)),2), sem(pow2db(prun/mean(prun)),2),'lineProps','r')


%% chronux

params.Fs=100;
params.fpass = [1 10];
params.trialave = 0;
params.pad = 0;
params.tapers = [3 3];
params.err = [1 0.05];

[Sstat,fstat,errstat] = mtspectrumc( allStat', params );
[Srun,frun,errrun] = mtspectrumc( allRun', params );

figure, hold on
plot(fstat,pow2db(mean(Sstat,2)/mean(mean(Sstat,2))),'k')
plot(frun,pow2db(mean(Srun,2)/mean(mean(Srun,2))),'r')



%% binned spikes %%

%% load data
load('M22029_20220607_r2overTime_dec22.mat')

%%
goodUnits = tempUnits([tempUnits.isi_viol]<=0.1...
    & [tempUnits.amplitude_cutoff]<=0.1 & [tempUnits.amplitude]>=50);

iplot=0;
for istate = 1:2
    for ispeed = 1:6
        iplot=iplot+1;
        nTrials = size(goodUnits(iunit).allSpikes{ispeed,istate},2);
        nTime = size(goodUnits(iunit).allSpikes{ispeed,istate},1);
        allSpikes = nan(numel(goodUnits),nTime,nTrials);
        for iunit = 1:numel(goodUnits)

            allSpikes(iunit,:,:) = goodUnits(iunit).allSpikes{ispeed,istate,:};
        end

        subplot(2,6,iplot)
        meanUnitSpikes = squeeze(mean(allSpikes,1));
        imagesc(meanUnitSpikes'), colorbar,colormap(1-gray)
        ylabel('trial')
        xlabel('time')

    end
end
%
% figure
% meanTrialSpikes = squeeze(mean(allSpikes,3));
% imagesc(meanTrialSpikes), colorbar,colormap(1-gray)
% ylabel('unit')
% xlabel('time')

%% chronux anaalysis of binned points process

sessionTags = {'M22027', '20220517';
    'M22029', '20220607';
    'M22031', '20220524';
    'M22032', '20220621';
    'M22033', '20220706'};

for isession=1:size(sessionTags,1)

    fname = [sessionTags{isession,1},'_', sessionTags{isession,2},'_r2overTime_dec22.mat'];

    load(fname)

goodUnits = tempUnits([tempUnits.isi_viol]<=0.1...
    & [tempUnits.amplitude_cutoff]<=0.1 & [tempUnits.amplitude]>=50);

clear t singleTrials
% figure
iplot=0;
for ispeed = 1:6
    iplot=iplot+1;
    %subplot(2,3,iplot), hold on
    for istate = 1:2
        nTrials = size(goodUnits(1).allSpikes{ispeed,istate},2);
        nTime = size(goodUnits(1).allSpikes{ispeed,istate},1);
        allSpikes = nan(numel(goodUnits),nTime,nTrials);
        for iunit = 1:numel(goodUnits)

            allSpikes(iunit,:,:) = goodUnits(iunit).allSpikes{ispeed,istate,:};
        end

        meanUnitSpikes = squeeze(mean(allSpikes,1));

        params.Fs=100;
        params.fpass = [1 10];
        params.trialave = 0;
        params.pad = 0;
        params.tapers = [3 3];
        params.err = [1 0.05];


        [S,f,R,Serr]=mtspectrumpb(meanUnitSpikes,params);
        meanS = mean(S,2);
        normS = meanS/mean(meanS);

%         plot(f,normS)

        t(ispeed,istate,:) = normS;
        singleTrials(ispeed,istate).pow = S;
        
    end
end


% plot average over all speeds
allStat = pow2db(squeeze(t(:,1,:)));
allRun = pow2db(squeeze(t(:,2,:)));
% allStat = squeeze(t(:,1,:));
% allRun = squeeze(t(:,2,:));

figure, hold on
shadedErrorBar(f,mean(allStat,1), sem(allStat,1),'lineProps','k')
shadedErrorBar(f,mean(allRun,1), sem(allRun,1),'lineProps','r')


s(isession).stat = mean(allStat,1);
s(isession).run = mean(allRun,1);
s(isession).t = t;
s(isession).singleTrials = singleTrials;
s(isession).statTrials = cat(2,singleTrials(:,1).pow);
s(isession).runTrials = cat(2,singleTrials(:,2).pow);
s(isession).allTrials = cat(2,s(isession).statTrials,s(isession).runTrials);
end

%% plot average of subjects

allStat = cat(1,s.stat);
allRun = cat(1,s.run);

figure, hold on
shadedErrorBar(f,mean(allStat,1), sem(allStat,1),'lineProps','k')
shadedErrorBar(f,mean(allRun,1), sem(allRun,1),'lineProps','r')
xlim([1 10]);
ylabel('Normalised power (db)')
xlabel('Frequency (Hz)')

%% by speed

speedcols=inferno(6);

for ispeed = 1:6
    allStat=[];
        allRun=[];

    for isesh = 1:5
        allStat=cat(1,allStat, squeeze(s(isesh).t(ispeed,1,:))');
        allRun=cat(1,allRun, squeeze(s(isesh).t(ispeed,2,:))');
    end
    subplot(2,3,ispeed)
    title(ispeed)
shadedErrorBar(f,mean(allStat,1), sem(allStat,1),'lineProps','k')
shadedErrorBar(f,mean(allRun,1), sem(allRun,1),'lineProps','r')
end


figure
for ispeed = 1:6
    allStat=[];
        allRun=[];

    for isesh = 1:5
        allStat=cat(1,allStat, squeeze(s(isesh).t(ispeed,1,:))');
        allRun=cat(1,allRun, squeeze(s(isesh).t(ispeed,2,:))');
    end

    subplot(121), hold on
    shadedErrorBar(f, mean(allStat), sem(allStat), 'lineProps',{'Color', speedcols(ispeed,:)});
    
    subplot(122), hold on
    shadedErrorBar(f, mean(allRun), sem(allRun), 'lineProps',{'Color', speedcols(ispeed,:)});

end


%% plot single trial powers

for isession = 1:5;
fRange = [2 6];
[~, idx(1)] = min(abs(f-fRange(1))); [~, idx(2)] = min(abs(f-fRange(2)));
idxRange = idx(1):idx(2);

fRangePower = mean(s(isession).statTrials(idxRange,:),1);
figure
histogram(fRangePower,'BinWidth',0.01) 
xlabel([num2str(fRange) ' rel pow.'])
title(num2str(isession))

end

%%
isession = 2;
fRange = [2 6];
[~, idx(1)] = min(abs(f-fRange(1))); [~, idx(2)] = min(abs(f-fRange(2)));
idxRange = idx(1):idx(2);

threshold = 0.055;
allTrials = s(isession).statTrials;
meanPerf = mean(allTrials,1);
relPerf = allTrials./meanPerf;

fRangePower = mean(relPerf(idxRange,:),1);

nTrials = numel(fRangePower);

quants = quantile(fRangePower,4);

%[~, maxindex] = maxk(fRangePower,40);
%[~, minindex] = mink(fRangePower,40);

minindex = find(fRangePower<quants(1));
maxindex = find(fRangePower>quants(4));

figure, 
subplot(211), hold on
plot(f,s(isession).allTrials(:,minindex),'m')
plot(f,s(isession).allTrials(:,maxindex),'c')

subplot(212), hold on
shadedErrorBar(f,mean(s(isession).allTrials(:,minindex),2), sem(s(isession).allTrials(:,minindex),2),'lineProps','m')
shadedErrorBar(f,mean(s(isession).allTrials(:,maxindex),2), sem(s(isession).allTrials(:,maxindex),2),'lineProps','c')

