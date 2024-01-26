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

end

% plot average of subjects

allStat = cat(1,s.stat);
allRun = cat(1,s.run);

figure, hold on
shadedErrorBar(f,mean(allStat,1), sem(allStat,1),'lineProps','k')
shadedErrorBar(f,mean(allRun,1), sem(allRun,1),'lineProps','r')
xlim([1 10]);
ylabel('Normalised power (db)')
xlabel('Frequency (Hz)')
