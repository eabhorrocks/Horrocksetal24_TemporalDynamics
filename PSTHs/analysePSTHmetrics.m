%% Plot and analyse psth metrics


sessionTags = {'M22027', '20220517';...
    'M22029', '20220607';...
    'M22031', '20220524';
    'M22032', '20220621';
    'M22033', '20220706'};

dataDir = 'C:\Users\edward.horrocks\Documents\Code\V1Dynamics\Data\basic_111022';


for isession = 1:size(sessionTags,1)
    fname = [sessionTags{isession,1},'_', sessionTags{isession,2},'_PSTH_strict_noDownSamp.mat']; %PSTH_noEye %psth3rds

    load(fullfile(dataDir,fname))
    [units.session] = deal(isession);
    session(isession).units = units;
    
    

end

%%

allUnits = cat(1,session([1 2 4 5]).units);
goodUnits = allUnits([allUnits.isi_viol]<=0.1...
    & [allUnits.amplitude_cutoff]<=0.1 & [allUnits.amplitude]>=50 & [allUnits.firing_rate]>=0);

%% get cell types

for iunit = 1:numel(goodUnits)
    [ccg, t] = CCG(goodUnits(iunit).spike_times, ones(size(goodUnits(iunit).spike_times)),...
    'binSize', 0.0005, 'duration', 0.1,'norm', 'rate');
    goodUnits(iunit).acg = ccg;
    fit_params_out = fit_ACG(ccg,false);

    goodUnits(iunit).tau_rise = fit_params_out.acg_tau_rise;
end

narrow_idx = find([goodUnits.duration]<=0.45);
wide_idx = find([goodUnits.duration]>0.45 & [goodUnits.tau_rise]>6);
pyr_idx = find(~ismember(1:numel(goodUnits), [narrow_idx,wide_idx]));

[goodUnits.cellType] = deal(nan);
[goodUnits(pyr_idx).cellType] = deal(1);
[goodUnits(narrow_idx).cellType] = deal(2);
[goodUnits(wide_idx).cellType] = deal(3);


%%

dopupilcomp = false;

stat.session = [];
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
stat.onsetCharVec = [];
stat.onset_bestR2Vec = [];
stat.onset_paramsVec = [];
stat.offsetCharVec = [];
stat.offset_bestR2Vec = [];
stat.offset_paramsVec = [];
stat.waveformdur=[];
stat.cellType=[];
stat.meanCorrecVec = [];
stat.peakCorrecVec = [];


run.session = [];
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
run.onsetCharVec = [];
run.onset_bestR2Vec = [];
run.onset_paramsVec = [];
run.offsetCharVec = [];
run.offset_bestR2Vec = [];
run.offset_paramsVec = [];
run.waveformdur=[];
run.cellType=[];
run.meanCorrecVec = [];
run.peakCorrecVec = [];

statSmall.session = [];
statSmall.speedVec = [];
statSmall.stateVec = [];
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
statSmall.waveformdur=[];
statSmall.onsetCharVec = [];
statSmall.onset_bestR2Vec = [];
statSmall.offsetCharVec = [];
statSmall.offset_bestR2Vec = [];


statBig.session = [];
statBig.speedVec = [];
statBig.stateVec = [];
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
statBig.waveformdur=[];
statBig.onsetCharVec = [];
statBig.onset_bestR2Vec = [];
statBig.offsetCharVec = [];
statBig.offset_bestR2Vec = [];

for iunit= 1:numel(goodUnits)

    for ispeed = 1:6

%         psths while stationary
        stat.session = cat(1,stat.session, goodUnits(iunit).session);
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

        % get baseline corrected mean and peak (max/min) FR
        if (goodUnits(iunit).stat_psth(ispeed).meanFR >= goodUnits(iunit).stat_blank.mu)
            stat.meanCorrecVec = cat(1, stat.meanCorrecVec, goodUnits(iunit).stat_psth(ispeed).meanFR - goodUnits(iunit).stat_blank.mu);
            stat.peakCorrecVec = cat(1, stat.peakCorrecVec, goodUnits(iunit).stat_psth(ispeed).maxFR - goodUnits(iunit).stat_blank.mu);
        else
            stat.meanCorrecVec = cat(1, stat.meanCorrecVec, goodUnits(iunit).stat_blank.mu - goodUnits(iunit).stat_psth(ispeed).meanFR);
            stat.peakCorrecVec = cat(1, stat.peakCorrecVec, goodUnits(iunit).stat_blank.mu - goodUnits(iunit).stat_psth(ispeed).minFR);
        end


        stat.onsetCharVec = cat(1, stat.onsetCharVec, goodUnits(iunit).stat_psth(ispeed).onset_char);
        stat.onset_bestR2Vec = cat(1, stat.onset_bestR2Vec, goodUnits(iunit).stat_psth(ispeed).onset_bestR2);
        stat.onset_paramsVec = cat(1, stat.onset_paramsVec, goodUnits(iunit).stat_psth(ispeed).onset_bestParams);
        stat.offsetCharVec = cat(1, stat.offsetCharVec, goodUnits(iunit).stat_psth(ispeed).offset_char);
        stat.offset_bestR2Vec = cat(1, stat.offset_bestR2Vec, goodUnits(iunit).stat_psth(ispeed).offset_bestR2);
        stat.offset_paramsVec = cat(1, stat.offset_paramsVec, goodUnits(iunit).stat_psth(ispeed).offset_bestParams);
        
        stat.waveformdur = cat(1, stat.waveformdur, goodUnits(iunit).duration);
%        stat.cellType = cat(1, stat.cellType, goodUnits(iunit).cellType);



%         psths while running
        run.session = cat(1,run.session, goodUnits(iunit).session);
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

        run.onsetCharVec = cat(1, run.onsetCharVec, goodUnits(iunit).run_psth(ispeed).onset_char);
        run.onset_bestR2Vec = cat(1, run.onset_bestR2Vec, goodUnits(iunit).run_psth(ispeed).onset_bestR2);
        run.onset_paramsVec = cat(1, run.onset_paramsVec, goodUnits(iunit).run_psth(ispeed).onset_bestParams);
        run.offsetCharVec = cat(1, run.offsetCharVec, goodUnits(iunit).run_psth(ispeed).offset_char);
        run.offset_bestR2Vec = cat(1, run.offset_bestR2Vec, goodUnits(iunit).run_psth(ispeed).offset_bestR2);
        run.offset_paramsVec = cat(1, run.offset_paramsVec, goodUnits(iunit).run_psth(ispeed).offset_bestParams);
        
        run.waveformdur = cat(1, run.waveformdur, goodUnits(iunit).duration);
     %    run.cellType = cat(1, run.cellType, goodUnits(iunit).cellType);

        if (goodUnits(iunit).stat_psth(ispeed).meanFR >= goodUnits(iunit).stat_blank.mu)
            run.meanCorrecVec = cat(1, run.meanCorrecVec, goodUnits(iunit).run_psth(ispeed).meanFR - goodUnits(iunit).run_blank.mu);
            run.peakCorrecVec = cat(1, run.peakCorrecVec, goodUnits(iunit).run_psth(ispeed).maxFR - goodUnits(iunit).stat_blank.mu);
        else
            run.meanCorrecVec = cat(1, run.meanCorrecVec, goodUnits(iunit).run_blank.mu - goodUnits(iunit).run_psth(ispeed).meanFR);
            run.peakCorrecVec = cat(1, run.peakCorrecVec, goodUnits(iunit).run_blank.mu - goodUnits(iunit).run_psth(ispeed).minFR);
        end


        if dopupilcomp
%         psths while stationary w/ small pupil
        statSmall.session = cat(1, statSmall.session, goodUnits(iunit).session);
        statSmall.speedVec = cat(1, statSmall.speedVec, ispeed);
        statSmall.stateVec = cat(1, statSmall.stateVec, 1);
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

        statSmall.onsetCharVec = cat(1, statSmall.onsetCharVec, goodUnits(iunit).statSmall_psth(ispeed).onset_char);
        statSmall.onset_bestR2Vec = cat(1, statSmall.onset_bestR2Vec, goodUnits(iunit).statSmall_psth(ispeed).onset_bestR2);
        statSmall.offsetCharVec = cat(1, statSmall.offsetCharVec, goodUnits(iunit).statSmall_psth(ispeed).offset_char);
        statSmall.offset_bestR2Vec = cat(1, statSmall.offset_bestR2Vec, goodUnits(iunit).statSmall_psth(ispeed).offset_bestR2);

%         psths while stationary w/ big pupil
        statBig.session = cat(1, statBig.session, goodUnits(iunit).session);
        statBig.speedVec = cat(1, statBig.speedVec, ispeed);
        statBig.stateVec = cat(1, statBig.stateVec, 1);
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

        statBig.onsetCharVec = cat(1, statBig.onsetCharVec, goodUnits(iunit).statBig_psth(ispeed).onset_char);
        statBig.onset_bestR2Vec = cat(1, statBig.onset_bestR2Vec, goodUnits(iunit).statBig_psth(ispeed).onset_bestR2);
        statBig.offsetCharVec = cat(1, statBig.offsetCharVec, goodUnits(iunit).statBig_psth(ispeed).offset_char);
        statBig.offset_bestR2Vec = cat(1, statBig.offset_bestR2Vec, goodUnits(iunit).statBig_psth(ispeed).offset_bestR2);
        end
    end
end

%% get reliably responsive

rangeThresh = 3;
reliThrsh =-1.645;
zThresh = 3.29;

statGoodVals =stat.reliVec<=reliThrsh & stat.rangeFRVec>=rangeThresh & stat.peakAbsZVec>=zThresh;
runGoodVals = run.reliVec<=reliThrsh & run.rangeFRVec>=rangeThresh & run.peakAbsZVec>=zThresh;
bothGoodVals = statGoodVals & runGoodVals;

statSmallGoodVals =statSmall.reliVec<=reliThrsh & statSmall.rangeFRVec>=rangeThresh & statSmall.peakAbsZVec>=zThresh;
statBigGoodVals =statBig.reliVec<=reliThrsh & statBig.rangeFRVec>=rangeThresh & statBig.peakAbsZVec>=zThresh;
bothStatGoodVals = statSmallGoodVals & statBigGoodVals;

%% plot reliable responses

A = [prop(statGoodVals), prop(runGoodVals)];
I = prop(bothGoodVals);
figure, venn(A,I);

statOnly = prop(statGoodVals & ~bothGoodVals);
runOnly = prop(runGoodVals & ~bothGoodVals);
both = prop(bothGoodVals);
figure, pie([statOnly, runOnly, both, 1-(statOnly+runOnly+both)])
legend({'statOnly', 'runOnly', 'Both', 'Neither'})

for ispeed = 1:6
    statProp(ispeed) = prop(statGoodVals(stat.speedVec==ispeed));
    runProp(ispeed) = prop(runGoodVals(run.speedVec==ispeed)); 
end
figure, hold on
plot(statProp, 'k')
plot(runProp, 'r')

%% loook at metrics across all speeds

metric_name = 'susIndexVec';
binEdges = 0:0.02:1;

statVals = stat.(metric_name);
runVals = run.(metric_name);



figure, hold on
scatter_kde(statVals(bothGoodVals), runVals(bothGoodVals), 'filled', 'MarkerSize', 6)
plot([binEdges(1) binEdges(end)], [binEdges(1) binEdges(end)], 'r')
xlabel('stat'), ylabel('run')


statSmallVals = statSmall.(metric_name);
statBigVals = statBig.(metric_name);

figure, hold on
scatter_kde(statSmallVals(bothStatGoodVals), statBigVals(bothStatGoodVals), 'filled', 'MarkerSize', 6)
plot([binEdges(1) binEdges(end)], [binEdges(1) binEdges(end)], 'r')
xlabel('stat small pupil'), ylabel('stat big pupil')

%% LME for sustainedness index

metric_name = 'susIndexVec';
binEdges = 0:0.02:1;

statVals = stat.(metric_name);
runVals = run.(metric_name);


valsVec = cat(1,statVals(bothGoodVals), runVals(bothGoodVals));
unitVec = categorical(cat(2,1:sum(bothGoodVals), 1:sum(bothGoodVals))');
stateVec = categorical(cat(1,repelem(1,sum(bothGoodVals),1), repelem(2,sum(bothGoodVals),1)));
seshVec = categorical(cat(1,stat.session(bothGoodVals), run.session(bothGoodVals)));

tbl = table(valsVec,unitVec,stateVec,seshVec,...
    'VariableNames',{'vals','unit','state','sesh'});

    f = 'vals ~ state + (1|unit) + (1|sesh)';
      
    lme = fitlme(tbl, f, 'DummyVarCoding', 'reference');

    %% LME for sustainedness index by pupil

metric_name = 'susIndexVec';
binEdges = 0:0.02:1;

 statSmallVals = statSmall.(metric_name);
 statBigVals = statBig.(metric_name);



valsVec = cat(1,statSmallVals(bothStatGoodVals), statBigVals(bothStatGoodVals));
unitVec = categorical(cat(2,1:sum(bothStatGoodVals), 1:sum(bothStatGoodVals))');
stateVec = categorical(cat(1,repelem(1,sum(bothStatGoodVals),1), repelem(2,sum(bothStatGoodVals),1)));
seshVec = categorical(cat(1,stat.session(bothStatGoodVals), run.session(bothStatGoodVals)));

tbl = table(valsVec,unitVec,stateVec,seshVec,...
    'VariableNames',{'vals','unit','state','sesh'});

    f = 'vals ~ state + (1|unit) + (1|sesh)';
      
    lme = fitlme(tbl, f, 'DummyVarCoding', 'reference');


%% changes in mean and peak firing rate

figure, 

metric_name = 'meanCorrecVec';
statVals = stat.(metric_name);
runVals = run.(metric_name);

subplot(221), hold on
scatter_kde(abs(statVals(bothGoodVals)), abs(runVals(bothGoodVals)),'filled','MarkerSize',3)
plot([0 50], [0 50], 'r')
xlim([0 50]), ylim([0 50])
axis equal
defaultAxesProperties(gca, true)
title('Mean Corrected FR')

subplot(223), hold on
meanFR_fracChange = abs(runVals(bothGoodVals))./abs(statVals(bothGoodVals));
histogram(log2(meanFR_fracChange), 'BinEdges',[-4.3219:0.1:4.3219], 'Normalization','Probability');
ylim([0 0.065])
plot([0 0], [0 0.065], 'r')
plot(median(log2(meanFR_fracChange)),0.06, 'v')
ax = gca; ax.XTick = sort([-1.*log2([1 2 5 10]),log2([2 5 10])]); ax.XTickLabel = [0.1 0.2, 0.5 1 2 5 10];
xlim([-4.3219, 4.3219])
defaultAxesProperties(gca, false)

metric_name = 'peakCorrecVec';
statVals = stat.(metric_name);
runVals = run.(metric_name);

subplot(222), hold on
scatter_kde(abs(statVals(bothGoodVals)), abs(runVals(bothGoodVals)),'filled','MarkerSize',3);
plot([0 80], [0 80], 'r')
xlim([0 80]), ylim([0 80])
axis equal
defaultAxesProperties(gca, true)
title('Peak Corrected FR')

subplot(224), hold on
peakFR_fracChange = abs(runVals(bothGoodVals))./abs(statVals(bothGoodVals));
histogram(log2(peakFR_fracChange), 'BinEdges',[-4.3219:0.1:4.3219], 'Normalization','Probability');
ylim([0 0.065])
plot([0 0], [0 0.065], 'r')
plot(median(log2(peakFR_fracChange)),0.06, 'v')
ax = gca; ax.XTick = sort([-1.*log2([1 2 5 10]),log2([2 5 10])]); ax.XTickLabel = [0.1 0.2, 0.5 1 2 5 10];
xlim([-4.3219, 4.3219])
defaultAxesProperties(gca, false)
colormap(viridis)

figure, hold on
scatter_kde(log2(meanFR_fracChange), log2(peakFR_fracChange),'filled','MarkerSize',3)
plot([-4.3219, 4.3219], [-4.3219, 4.3219], 'r')
axis equal
 plot([-4.3219, 4.3219],log2([1 1]),  'k:')
plot(log2([1 1]), [-4.3219, 4.3219], 'k:')

xlim([-4.3219, 4.3219])
ylim([-4.3219, 4.3219])
ax = gca; ax.XTick = sort([-1.*log2([1 2 5 10]),log2([2 5 10])]); ax.YTick = ax.XTick;
ax.XTickLabel = [0.1 0.2, 0.5 1 2 5 10]; ax.YTickLabel = ax.XTickLabel;
defaultAxesProperties(gca, false)
colormap(viridis)


p = signrank(log2(meanFR_fracChange), log2(peakFR_fracChange))





%% onset and offset features plots

doPupilComp=false;
% onset features
statVals = stat.onsetCharVec(statGoodVals);
runVals = run.onsetCharVec(runGoodVals);
statSessions = stat.session(statGoodVals);
runSessions = run.session(runGoodVals);

if doPupilComp
statVals = statSmall.onsetCharVec(statSmallGoodVals);
runVals = statBig.onsetCharVec(statBigGoodVals);
statSessions = statSmall.session(statSmallGoodVals);
runSessions = statBig.session(statBigGoodVals);
end

for ifet = 1:5
    statProp_onset(ifet) = prop(statVals==ifet);
    runProp_onset(ifet) = prop(runVals==ifet);
end

for isession = 1:5
    tempStatVals = statVals(statSessions==isession);
    tempRunVals = runVals(runSessions==isession);
    for ifet = 1:5
        statProp_s_onset(ifet,isession) = prop(tempStatVals==ifet);
        runProp_s_onset(ifet,isession) = prop(tempRunVals==ifet);
    end
end

figure, hold on
% backgpound bar chart
% b1 = bar([statProp_onset; runProp_onset]');
% b1(1).EdgeColor='k'; b1(1).FaceAlpha=0;
% b1(2).EdgeColor='r'; b1(2).FaceAlpha=0;

% data plot
for ifet = 1:5
    for isession = 1:5
    plot(ifet-0.15, statProp_s_onset(ifet,isession), 'ko')
    plot(ifet+0.15, runProp_s_onset(ifet,isession), 'ro')
    plot([ifet-0.15, ifet+0.15],...
        [statProp_s_onset(ifet,isession), runProp_s_onset(ifet,isession)],...
        'Color', [.7 .7 .7])
    plot([ifet-0.3, ifet-0],...
        [mean(statProp_s_onset(ifet,:)), mean(statProp_s_onset(ifet,:))],...
        'Color', 'k')
    plot([ifet+0.0, ifet+0.3],...
        [mean(runProp_s_onset(ifet,:)), mean(runProp_s_onset(ifet,:))],...
        'Color', 'r')
    end
    [onset_p(ifet)] = signrank(statProp_s_onset(ifet,:), runProp_s_onset(ifet,:));

end
title('Onset Features')
ylim([0, 0.8])
defaultAxesProperties(gca, true);



% offset features
statVals = stat.offsetCharVec(statGoodVals);
runVals = run.offsetCharVec(runGoodVals);
statSessions = stat.session(statGoodVals);
runSessions = run.session(runGoodVals);

if doPupilComp
statVals = statSmall.offsetCharVec(statSmallGoodVals);
runVals = statBig.offsetCharVec(statBigGoodVals);
statSessions = statSmall.session(statSmallGoodVals);
runSessions = statBig.session(statBigGoodVals);
end

for ifet = 1:5
    statProp_offset(ifet) = prop(statVals==ifet);
    runProp_offset(ifet) = prop(runVals==ifet);
end

for isession = 1:5
    tempStatVals = statVals(statSessions==isession);
    tempRunVals = runVals(runSessions==isession);
    for ifet = 1:5
        statProp_s_offset(ifet,isession) = prop(tempStatVals==ifet);
        runProp_s_offset(ifet,isession) = prop(tempRunVals==ifet);
    end
end

figure, hold on
% backgpound bar chart
% b1 = bar([statProp_onset; runProp_onset]');
% b1(1).EdgeColor='k'; b1(1).FaceAlpha=0;
% b1(2).EdgeColor='r'; b1(2).FaceAlpha=0;

% data plot
for ifet = 1:5
    for isession = 1:5
    plot(ifet-0.15, statProp_s_offset(ifet,isession), 'ko')
    plot(ifet+0.15, runProp_s_offset(ifet,isession), 'ro')
    plot([ifet-0.15, ifet+0.15],...
        [statProp_s_offset(ifet,isession), runProp_s_offset(ifet,isession)],...
        'Color', [.7 .7 .7])
    plot([ifet-0.3, ifet-0],...
        [mean(statProp_s_offset(ifet,:)), mean(statProp_s_offset(ifet,:))],...
        'Color', 'k')
    plot([ifet+0.0, ifet+0.3],...
        [mean(runProp_s_offset(ifet,:)), mean(runProp_s_offset(ifet,:))],...
        'Color', 'r')
    end
    [offset_p(ifet)] = signrank(statProp_s_offset(ifet,:), runProp_s_offset(ifet,:));

end
title('Offset Features')
ylim([0, 0.8])
defaultAxesProperties(gca, true);



%% onset and offset for pupil + running

%% onset and offset features

% onset features
statSmallVals = statSmall.onsetCharVec(statSmallGoodVals);
statBigVals = statBig.onsetCharVec(statBigGoodVals);
runVals = run.onsetCharVec(runGoodVals);

statSmallSessions = statSmall.session(statSmallGoodVals);
statBigSessions = statBig.session(statBigGoodVals);
runSessions = run.session(runGoodVals);

for ifet = 1:5
    statSmallProp_onset(ifet) = prop(statSmallVals==ifet);
    statBigProp_onset(ifet) = prop(statBigVals==ifet);
    runProp_onset(ifet) = prop(runVals==ifet);
end

for isession = 1:5
    tempStatSmallVals = statSmallVals(statSmallSessions==isession);
    tempStatBigVals = statBigVals(statBigSessions==isession);
    tempRunVals = runVals(runSessions==isession);
    for ifet = 1:5
        statSmallProp_s_onset(ifet,isession) = prop(tempStatSmallVals==ifet);
        statBigProp_s_onset(ifet,isession) = prop(tempStatBigVals==ifet);
        runProp_s_onset(ifet,isession) = prop(tempRunVals==ifet);
    end
end

figure, hold on
% backgpound bar chart
%  b1 = bar([statSmallProp_onset; statBigProp_onset; runProp_onset]');
%  b1(1).EdgeColor='g'; b1(1).FaceAlpha=0;
%  b1(2).EdgeColor='c'; b1(2).FaceAlpha=0;
%  b1(3).EdgeColor='r'; b1(3).FaceAlpha=0;

%
% data plot
for ifet = 1:5
    for isession = 1:5
    plot(ifet-0.15, statSmallProp_s_onset(ifet,isession), 'go')
    plot(ifet, statBigProp_s_onset(ifet,isession), 'co')
    plot(ifet+0.15, runProp_s_onset(ifet,isession), 'ro')
    plot([ifet-0.15, ifet, ifet+0.15],...
        [statSmallProp_s_onset(ifet,isession), statBigProp_s_onset(ifet,isession),...
        runProp_s_onset(ifet,isession)],'Color', [.7 .7 .7])




    plot([ifet-0.2, ifet-0.1],...
        [mean(statSmallProp_s_onset(ifet,:)), mean(statSmallProp_s_onset(ifet,:))],...
        'Color', 'g')
    plot([ifet-0.05, ifet+0.05],...
        [mean(statBigProp_s_onset(ifet,:)), mean(statBigProp_s_onset(ifet,:))],...
        'Color', 'c')
    plot([ifet+0.1, ifet+0.2],...
        [mean(runProp_s_onset(ifet,:)), mean(runProp_s_onset(ifet,:))],...
        'Color', 'r')
     end
    [onset_p(ifet)] = signrank(statProp_s_onset(ifet,:), runProp_s_onset(ifet,:));

end
title('Onset Features')
ylim([0, 0.8])
defaultAxesProperties(gca, true);


%% offset features
statSmallVals = statSmall.offsetCharVec(statSmallGoodVals);
statBigVals = statBig.offsetCharVec(statBigGoodVals);
runVals = run.offsetCharVec(runGoodVals);

statSmallSessions = statSmall.session(statSmallGoodVals);
statBigSessions = statBig.session(statBigGoodVals);
runSessions = run.session(runGoodVals);

for ifet = 1:5
    statSmallProp_onset(ifet) = prop(statSmallVals==ifet);
    statBigProp_onset(ifet) = prop(statBigVals==ifet);
    runProp_onset(ifet) = prop(runVals==ifet);
end

for isession = 1:5
    tempStatSmallVals = statSmallVals(statSmallSessions==isession);
    tempStatBigVals = statBigVals(statBigSessions==isession);
    tempRunVals = runVals(runSessions==isession);
    for ifet = 1:5
        statSmallProp_s_onset(ifet,isession) = prop(tempStatSmallVals==ifet);
        statBigProp_s_onset(ifet,isession) = prop(tempStatBigVals==ifet);
        runProp_s_onset(ifet,isession) = prop(tempRunVals==ifet);
    end
end

figure, hold on
% backgpound bar chart
%  b1 = bar([statSmallProp_onset; statBigProp_onset; runProp_onset]');
%  b1(1).EdgeColor='g'; b1(1).FaceAlpha=0;
%  b1(2).EdgeColor='c'; b1(2).FaceAlpha=0;
%  b1(3).EdgeColor='r'; b1(3).FaceAlpha=0;

%
% data plot
for ifet = 1:5
    for isession = 1:5
    plot(ifet-0.15, statSmallProp_s_onset(ifet,isession), 'go')
    plot(ifet, statBigProp_s_onset(ifet,isession), 'co')
    plot(ifet+0.15, runProp_s_onset(ifet,isession), 'ro')
    plot([ifet-0.15, ifet, ifet+0.15],...
        [statSmallProp_s_onset(ifet,isession), statBigProp_s_onset(ifet,isession),...
        runProp_s_onset(ifet,isession)],'Color', [.7 .7 .7])




    plot([ifet-0.2, ifet-0.1],...
        [mean(statSmallProp_s_onset(ifet,:)), mean(statSmallProp_s_onset(ifet,:))],...
        'Color', 'g')
    plot([ifet-0.05, ifet+0.05],...
        [mean(statBigProp_s_onset(ifet,:)), mean(statBigProp_s_onset(ifet,:))],...
        'Color', 'c')
    plot([ifet+0.1, ifet+0.2],...
        [mean(runProp_s_onset(ifet,:)), mean(runProp_s_onset(ifet,:))],...
        'Color', 'r')
     end
    [onset_p(ifet)] = signrank(statProp_s_onset(ifet,:), runProp_s_onset(ifet,:));

end
title('Offset Features')
ylim([0, 0.8])
defaultAxesProperties(gca, true);

%%

% offset features
statVals = stat.offsetCharVec(statGoodVals);
runVals = run.offsetCharVec(runGoodVals);
statSessions = stat.session(statGoodVals);
runSessions = run.session(runGoodVals);

if doPupilComp
statVals = statSmall.offsetCharVec(statSmallGoodVals);
runVals = statBig.offsetCharVec(statBigGoodVals);
statSessions = statSmall.session(statSmallGoodVals);
runSessions = statBig.session(statBigGoodVals);
end

for ifet = 1:5
    statProp_offset(ifet) = prop(statVals==ifet);
    runProp_offset(ifet) = prop(runVals==ifet);
end

for isession = 1:5
    tempStatVals = statVals(statSessions==isession);
    tempRunVals = runVals(runSessions==isession);
    for ifet = 1:5
        statProp_s_offset(ifet,isession) = prop(tempStatVals==ifet);
        runProp_s_offset(ifet,isession) = prop(tempRunVals==ifet);
    end
end

figure, hold on
% backgpound bar chart
% b1 = bar([statProp_onset; runProp_onset]');
% b1(1).EdgeColor='k'; b1(1).FaceAlpha=0;
% b1(2).EdgeColor='r'; b1(2).FaceAlpha=0;

% data plot
for ifet = 1:5
    for isession = 1:5
    plot(ifet-0.15, statProp_s_offset(ifet,isession), 'ko')
    plot(ifet+0.15, runProp_s_offset(ifet,isession), 'ro')
    plot([ifet-0.15, ifet+0.15],...
        [statProp_s_offset(ifet,isession), runProp_s_offset(ifet,isession)],...
        'Color', [.7 .7 .7])
    plot([ifet-0.3, ifet-0],...
        [mean(statProp_s_offset(ifet,:)), mean(statProp_s_offset(ifet,:))],...
        'Color', 'k')
    plot([ifet+0.0, ifet+0.3],...
        [mean(runProp_s_offset(ifet,:)), mean(runProp_s_offset(ifet,:))],...
        'Color', 'r')
    end
    [offset_p(ifet)] = signrank(statProp_s_offset(ifet,:), runProp_s_offset(ifet,:));

end
title('Offset Features')
ylim([0, 0.8])
defaultAxesProperties(gca, true);








%% STATS for response features

% binomial glme for each feature individually

% feature ~ 1 + state + (1|session)

%%% ONSET
statVals = statSmall.onsetCharVec(statSmallGoodVals);
runVals = run.onsetCharVec(runGoodVals);

stateVec = categorical(cat(1,repelem(1,numel(statVals),1), repelem(2,numel(runVals),1)));
fetVec = categorical(cat(1,statVals,runVals));
statSessions = statSmall.session(statSmallGoodVals);
runSessions = run.session(runGoodVals);
sessionVec = categorical(cat(1,statSessions, runSessions));


basetbl = table(fetVec,stateVec,sessionVec,...
    'VariableNames',{'fet','state','session'});


for ifet = 1:5
    thisTbl = basetbl;
    thisTbl.fet = thisTbl.fet==num2str(ifet);

    f = 'fet ~ state + (1|session)';
      
    glme = fitglme(thisTbl, f, 'Distribution', 'Binomial',...
        'DummyVarCoding','reference','FitMethod','Laplace', 'link', 'logit');

    pVal_onset(ifet) = glme.Coefficients.pValue(2);
end


%%% OFFSET
statVals = statSmall.offsetCharVec(statSmallGoodVals);
runVals = run.offsetCharVec(runGoodVals);

stateVec = categorical(cat(1,repelem(1,numel(statVals),1), repelem(2,numel(runVals),1)));
fetVec = categorical(cat(1,statVals,runVals));
statSessions = statSmall.session(statSmallGoodVals);
runSessions = run.session(runGoodVals);
sessionVec = categorical(cat(1,statSessions, runSessions));


basetbl = table(fetVec,stateVec,sessionVec,...
    'VariableNames',{'fet','state','session'});


for ifet = 1:5
    thisTbl = basetbl;
    thisTbl.fet = thisTbl.fet==num2str(ifet);

    f = 'fet ~ state + (1|session)';
      
    glme = fitglme(thisTbl, f, 'Distribution', 'Binomial',...
        'DummyVarCoding','reference','FitMethod','Laplace', 'link', 'logit');

    pVal_offset(ifet) = glme.Coefficients.pValue(2);
end


%% get basic props
statVals = stat.onsetCharVec(statGoodVals);
runVals = run.onsetCharVec(runGoodVals);

for ifet = 1:5
    propStatOn(ifet) = prop(statVals==ifet);
    propRunOn(ifet) = prop(runVals==ifet);
end

statVals = stat.offsetCharVec(statGoodVals);
runVals = run.offsetCharVec(runGoodVals);


for ifet = 1:5
    propStatOff(ifet) = prop(statVals==ifet);
    propRunOff(ifet) = prop(runVals==ifet);
end



%% STATS for response features, by pupil

% binomial glme for each feature individually

% feature ~ 1 + state + (1|session)

%%% ONSET
statSmallVals = statSmall.onsetCharVec(statSmallGoodVals);
statBigVals = statBig.onsetCharVec(statBigGoodVals);


statSmallSessions = statSmall.session(statSmallGoodVals);
statBigSessions = statBig.session(statBigGoodVals);


stateVec = categorical(cat(1,repelem(1,numel(statSmallVals),1), repelem(2,numel(statBigVals),1)));
fetVec = categorical(cat(1,statSmallVals,statBigVals));
sessionVec = categorical(cat(1,statSmallSessions, statBigSessions));


basetbl = table(fetVec,stateVec,sessionVec,...
    'VariableNames',{'fet','state','session'});


for ifet = 1:5
    thisTbl = basetbl;
    thisTbl.fet = thisTbl.fet==num2str(ifet);

    f = 'fet ~ state + (1|session)';
      
    glme = fitglme(thisTbl, f, 'Distribution', 'Binomial',...
        'DummyVarCoding','reference','FitMethod','Laplace', 'link', 'logit');

    pVal_onset(ifet) = glme.Coefficients.pValue(2);
end


%%% ONSET
statSmallVals = statSmall.offsetCharVec(statSmallGoodVals);
statBigVals = statBig.offsetCharVec(statBigGoodVals);

stateVec = categorical(cat(1,repelem(1,numel(statSmallVals),1), repelem(2,numel(statBigVals),1)));
fetVec = categorical(cat(1,statSmallVals,statBigVals));
sessionVec = categorical(cat(1,statSmallSessions, statBigSessions));


basetbl = table(fetVec,stateVec,sessionVec,...
    'VariableNames',{'fet','state','session'});


for ifet = 1:5
    thisTbl = basetbl;
    thisTbl.fet = thisTbl.fet==num2str(ifet);

    f = 'fet ~ state + (1|session)';
      
    glme = fitglme(thisTbl, f, 'Distribution', 'Binomial',...
        'DummyVarCoding','reference','FitMethod','Laplace', 'link', 'logit');

    pVal_offset(ifet) = glme.Coefficients.pValue(2);
end

%% get basic props
statSmallVals = statSmall.onsetCharVec(statSmallGoodVals);
statBigVals = statBig.onsetCharVec(statBigGoodVals);

for ifet = 1:5
    propSmallOn(ifet) = prop(statSmallVals==ifet);
    propBigOn(ifet) = prop(statBigVals==ifet);
end

statSmallVals = statSmall.offsetCharVec(statSmallGoodVals);
statBigVals = statBig.offsetCharVec(statBigGoodVals);

for ifet = 1:5
    propSmallOff(ifet) = prop(statSmallVals==ifet);
    propBigOff(ifet) = prop(statBigVals==ifet);
end


%% plot example PSTH w/ onset/offset feature

figure
idx = 726;

gaussFun =  @(params,xdata) params(1) + params(2).*exp(-(((xdata-params(3)).^2)/(2*(params(4).^2))));


allP = stat.psthVec(statGoodVals,:);
ps = normalize(allP(idx,:),'range');

onsetParams = stat.onset_paramsVec(statGoodVals,:);
offsetParams = stat.offset_paramsVec(statGoodVals,:);

plot(ps, 'k'), hold on

plot(21:51, feval(gaussFun,[onsetParams(idx,:)], 1:31), 'r-')
plot(121:151, feval(gaussFun,[offsetParams(idx,:)], 1:31), 'r-')
defaultAxesProperties(gca, true);
ax = gca; ax.XTick = [1:20:200, 200];


