%% analyse allen PSTHs

areas = {'VISp', 'VISl', 'VISal', 'VISrl', 'VISam', 'VISpm', 'LGd', 'LP'};
minTrials = 10;

%% load data
dataDirectory = 'E:\AllenInstituteAnalyses\AllenDataAnalysis\Data\rawMatFiles';
dataFiles = dir(fullfile(dataDirectory, 'psth_session*'));

for isession = 1:numel(dataFiles)


    load(fullfile([dataDirectory,'\' dataFiles(isession).name]))
    [units.session] = deal(isession);

    session(isession).units = units;

end

%% get good units
allUnits = cat(2,session.units);
goodUnits = allUnits([allUnits.isi_violations]<=0.1...
    & [allUnits.amplitude_cutoff]<=0.1 & [allUnits.waveform_amplitude]>=50 & [allUnits.firing_rate]>=0);


%% get metrics

stat.session = [];
stat.unit = [];
stat.speedVec = [];
stat.areaVec = [];
stat.dirVec = [];
stat.stateVec = [];
stat.psthVec = [];
stat.reliVec = [];
stat.rangeFRVec = [];
stat.peakAbsZVec = [];

stat.susIndexVec = [];
stat.onsetCharVec = [];
stat.onset_bestR2Vec = [];
stat.onset_paramsVec = [];
stat.offsetCharVec = [];
stat.offset_bestR2Vec = [];
stat.offset_paramsVec = [];


run.session = [];
run.unit = [];
run.speedVec = [];
run.areaVec = [];
run.dirVec = [];
run.stateVec = [];
run.psthVec = [];
run.reliVec = [];
run.rangeFRVec = [];
run.peakAbsZVec = [];

run.susIndexVec = [];
run.onsetCharVec = [];
run.onset_bestR2Vec = [];
run.onset_paramsVec = [];
run.offsetCharVec = [];
run.offset_bestR2Vec = [];
run.offset_paramsVec = [];


for iunit= 1:numel(goodUnits)
    if ismember(goodUnits(iunit).ecephys_structure_acronym, areas)
        for idir = 1:4
            for ispeed = 1:6
                % stat
                if goodUnits(iunit).statPSTH(idir,ispeed).nTrials>=minTrials

                    % stat
                    stat.session = cat(1,stat.session, goodUnits(iunit).session);
                    stat.unit = cat(1,stat.unit, goodUnits(iunit).ID);
                    stat.speedVec = cat(1,stat.speedVec, ispeed);
                    stat.areaVec = cat(1, stat.areaVec, goodUnits(iunit).ecephys_structure_acronym);
                    stat.dirVec = cat(1,stat.dirVec, idir);
                    stat.stateVec = cat(1,stat.stateVec, 1);
                    stat.psthVec = cat(1,stat.psthVec, goodUnits(iunit).statPSTH(idir,ispeed).psth);
                    stat.reliVec = cat(1,stat.reliVec,goodUnits(iunit).statPSTH(idir,ispeed).reliability_z);
                    stat.rangeFRVec = cat(1,stat.rangeFRVec, goodUnits(iunit).statPSTH(idir,ispeed).rangeFR);
                    stat.peakAbsZVec = cat(1,stat.peakAbsZVec, goodUnits(iunit).statPSTH(idir,ispeed).peakAbsZ);

                    stat.susIndexVec = cat(1,stat.susIndexVec, goodUnits(iunit).statPSTH(idir,ispeed).susIdx);
                    stat.onsetCharVec = cat(1,stat.onsetCharVec, goodUnits(iunit).statPSTH(idir,ispeed).onset_char);
                    stat.onset_bestR2Vec = cat(1,stat.onset_bestR2Vec, goodUnits(iunit).statPSTH(idir,ispeed).onset_bestR2);
                    stat.onset_paramsVec = cat(1,stat.onset_paramsVec, goodUnits(iunit).statPSTH(idir,ispeed).onset_bestParams);
                    stat.offsetCharVec = cat(1,stat.offsetCharVec, goodUnits(iunit).statPSTH(idir,ispeed).offset_char);
                    stat.offset_bestR2Vec = cat(1,stat.offset_bestR2Vec, goodUnits(iunit).statPSTH(idir,ispeed).offset_bestR2);
                    stat.offset_paramsVec = cat(1,stat.offset_paramsVec, goodUnits(iunit).statPSTH(idir,ispeed).offset_bestParams);
                end
                if goodUnits(iunit).runPSTH(idir,ispeed).nTrials>=minTrials
                    % run
                    run.session = cat(1,run.session, goodUnits(iunit).session);
                    run.unit = cat(1,run.unit, goodUnits(iunit).ID);
                    run.speedVec = cat(1,run.speedVec, ispeed);
                    run.areaVec = cat(1, run.areaVec, goodUnits(iunit).ecephys_structure_acronym);

                    run.dirVec = cat(1,run.dirVec, idir);
                    run.stateVec = cat(1,run.stateVec, 1);
                    run.psthVec = cat(1,run.psthVec, goodUnits(iunit).runPSTH(idir,ispeed).psth);
                    run.reliVec = cat(1,run.reliVec,goodUnits(iunit).runPSTH(idir,ispeed).reliability_z);
                    run.rangeFRVec = cat(1,run.rangeFRVec, goodUnits(iunit).runPSTH(idir,ispeed).rangeFR);
                    run.peakAbsZVec = cat(1,run.peakAbsZVec, goodUnits(iunit).runPSTH(idir,ispeed).peakAbsZ);

                    run.susIndexVec = cat(1,run.susIndexVec, goodUnits(iunit).runPSTH(idir,ispeed).susIdx);
                    run.onsetCharVec = cat(1,run.onsetCharVec, goodUnits(iunit).runPSTH(idir,ispeed).onset_char);
                    run.onset_bestR2Vec = cat(1,run.onset_bestR2Vec, goodUnits(iunit).runPSTH(idir,ispeed).onset_bestR2);
                    run.onset_paramsVec = cat(1,run.onset_paramsVec, goodUnits(iunit).runPSTH(idir,ispeed).onset_bestParams);
                    run.offsetCharVec = cat(1,run.offsetCharVec, goodUnits(iunit).runPSTH(idir,ispeed).offset_char);
                    run.offset_bestR2Vec = cat(1,run.offset_bestR2Vec, goodUnits(iunit).runPSTH(idir,ispeed).offset_bestR2);
                    run.offset_paramsVec = cat(1,run.offset_paramsVec, goodUnits(iunit).runPSTH(idir,ispeed).offset_bestParams);
                end
            end
        end
    end
end

%% get reliably responsive

rangeThresh = 3;
reliThrsh =-1.645;
zThresh = 3.29;

statGoodVals =stat.reliVec<=reliThrsh & stat.rangeFRVec>=rangeThresh & stat.peakAbsZVec>=zThresh;
runGoodVals = run.reliVec<=reliThrsh & run.rangeFRVec>=rangeThresh & run.peakAbsZVec>=zThresh;

%% mean responses

for iarea = 1:8

    statVals = stat.psthVec(statGoodVals & strcmp(stat.areaVec, areas(iarea)),:);
    runVals = run.psthVec(runGoodVals & strcmp(run.areaVec, areas(iarea)),:);

    figure, hold on
    shadedErrorBar(-195:10:1795, mean(statVals,1),sem(statVals,1))
    shadedErrorBar(-195:10:1795, mean(runVals,1),sem(runVals,1), 'lineProps','r')
    title(areas(iarea))
    ylim([0 24])
    defaultAxesProperties(gca, true)
    ax = gca; ax.XTick = -200:200:1800; ax.YTick = [0:4:24];
end




%% sustainedness index
figure
for iarea = 1:8

    statVals = stat.susIndexVec(statGoodVals & strcmp(stat.areaVec, areas(iarea)));
    runVals = run.susIndexVec(runGoodVals & strcmp(run.areaVec, areas(iarea)));

    subplot(2,4,iarea), hold on

    cols_cellArray = {[0 0 0], [1 0 1]};
    allVals = {statVals, runVals};
    xvals = [1, 2];

    % violin plot
    distributionPlot(allVals,'globalNorm',3,'histOpt',1,'divFactor',3,...
        'xValues', xvals, 'addSpread', 0, 'distWidth', 0.8, 'Color', cols_cellArray,...
        'showMM',6)
    ylim([0 1])
    defaultAxesProperties(gca, true)




    % bar([1 2], [mean(statVals(:)), mean(runVals(:))]);
    % errorbar([1 2], [mean(statVals(:)), mean(runVals(:))], [sem(statVals(:)), sem(runVals(:))])

end



%% plot as scatter plot

areaCols = tab20;
areaCols = areaCols(1:2:17,:);

figure, hold on, 
xlim([0 1])
ylim([0 1])
plot([0 1],[0 1],'k:')
for iarea = 1:8

    statVals = stat.susIndexVec(statGoodVals & strcmp(stat.areaVec, areas(iarea)));
    runVals = run.susIndexVec(runGoodVals & strcmp(run.areaVec, areas(iarea)));

    meanStat = mean(statVals);
    semStat = sem(statVals);

    meanRun=mean(runVals);
    semRun=sem(runVals);

    errorbar(meanStat,meanRun,semRun,semRun,semStat,semStat, 'Color',areaCols(iarea,:))

  

end

%% onset and offset features


plotFlag = false;

for iarea = 1:8
    figure
    statVals = stat.onsetCharVec(statGoodVals & strcmp(stat.areaVec, areas(iarea)));
    runVals = run.onsetCharVec(runGoodVals & strcmp(run.areaVec, areas(iarea)));

    statSessions = stat.session(statGoodVals & strcmp(stat.areaVec, areas(iarea)));
    runSessions = run.session(runGoodVals & strcmp(run.areaVec, areas(iarea)));



        % stats
stateVec = categorical(cat(1,repelem(1,numel(statVals),1), repelem(2,numel(runVals),1)));
fetVec = categorical(cat(1,statVals,runVals));
sessionVec = categorical(cat(1,statSessions, runSessions));

basetbl = table(fetVec,stateVec,sessionVec,...
    'VariableNames',{'fet','state','session'});

for ifet = 1:5
    thisTbl = basetbl;
    thisTbl.fet = thisTbl.fet==num2str(ifet);

    f = 'fet ~ state + (1|session)';
      
    glme = fitglme(thisTbl, f, 'Distribution', 'Binomial',...
        'DummyVarCoding','reference','FitMethod','Laplace', 'link', 'logit');

    ta(iarea).pVal_onset(ifet) = glme.Coefficients.pValue(2);
end



        ta(iarea).nStat = numel(statVals);
    ta(iarea).nRun = numel(runVals);

    statProp_onset = [];
    runProp_onset = [];


    for ifet = 1:5
        statProp_onset(ifet) = prop(statVals==ifet);
        runProp_onset(ifet) = prop(runVals==ifet);
    end

    ta(iarea).statProp_onset = statProp_onset;
    ta(iarea).runProp_onset = runProp_onset;

    subplot(121), hold on
    % backgpound bar chart
    b1 = bar([statProp_onset; runProp_onset]');
    b1(1).EdgeColor='k'; b1(1).FaceAlpha=0;
    b1(2).EdgeColor='r'; b1(2).FaceAlpha=0;

    statProp_onsetSesh = [];
    runProp_onsetSesh = [];
    for isession = 1:numel(dataFiles)
        statVals = stat.onsetCharVec(statGoodVals & strcmp(stat.areaVec, areas(iarea)) & [stat.session]==isession);
        runVals = run.onsetCharVec(runGoodVals & strcmp(run.areaVec, areas(iarea))& [run.session]==isession);
        for ifet = 1:5
            statProp_onsetSesh(isession,ifet) = prop(statVals==ifet);
            runProp_onsetSesh(isession, ifet) = prop(runVals==ifet);
        end
    end




    


    ta(iarea).nStatSesh = sum(~isnan(sum(statProp_onsetSesh,2)));
    ta(iarea).nRunSesh = sum(~isnan(sum(runProp_onsetSesh,2)));


    if plotFlag

    for ifet = 1:5
        for isession = 1:24
            plot(ifet-0.1429, statProp_onsetSesh(isession,ifet),'ko');
            plot(ifet+0.1429, runProp_onsetSesh(isession,ifet),'ro');
        end
    end
    end


    ylim([0, 0.65])

    title(areas(iarea))
    defaultAxesProperties(gca, false)


    statVals = stat.offsetCharVec(statGoodVals & strcmp(stat.areaVec, areas(iarea)));
    runVals = run.offsetCharVec(runGoodVals & strcmp(run.areaVec, areas(iarea)));

            % stats
stateVec = categorical(cat(1,repelem(1,numel(statVals),1), repelem(2,numel(runVals),1)));
fetVec = categorical(cat(1,statVals,runVals));
sessionVec = categorical(cat(1,statSessions, runSessions));

basetbl = table(fetVec,stateVec,sessionVec,...
    'VariableNames',{'fet','state','session'});

for ifet = 1:5
    thisTbl = basetbl;
    thisTbl.fet = thisTbl.fet==num2str(ifet);

    f = 'fet ~ state + (1|session)';
      
    glme = fitglme(thisTbl, f, 'Distribution', 'Binomial',...
        'DummyVarCoding','reference','FitMethod','Laplace', 'link', 'logit');

    ta(iarea).pVal_offset(ifet) = glme.Coefficients.pValue(2);
end

    statProp_offset = [];
    runProp_offset = [];

    for ifet = 1:5
        statProp_offset(ifet) = prop(statVals==ifet);
        runProp_offset(ifet) = prop(runVals==ifet);
    end

    ta(iarea).statProp_offset = statProp_onset;
    ta(iarea).runProp_offset = runProp_onset;

    subplot(122), hold on
    % backgpound bar chart
    b1 = bar([statProp_offset; runProp_offset]');
    b1(1).EdgeColor='k'; b1(1).FaceAlpha=0;
    b1(2).EdgeColor='r'; b1(2).FaceAlpha=0;
    ylim([0, 0.65])
    title(['#stat = ', num2str(numel(statVals)),',',' #run = ', num2str(numel(runVals))])
    defaultAxesProperties(gca, false)


    statProp_offsetSesh = [];
    runProp_offsetSesh = [];
    for isession = 1:numel(dataFiles)
        statVals = stat.offsetCharVec(statGoodVals & strcmp(stat.areaVec, areas(iarea)) & [stat.session]==isession);
        runVals = run.offsetCharVec(runGoodVals & strcmp(run.areaVec, areas(iarea))& [run.session]==isession);
        for ifet = 1:5
            statProp_offsetSesh(isession,ifet) = prop(statVals==ifet);
            runProp_offsetSesh(isession, ifet) = prop(runVals==ifet);
        end
    end

    if plotFlag
    for ifet = 1:5
        for isession = 1:24
            plot(ifet-0.1429, statProp_offsetSesh(isession,ifet),'ko');
            plot(ifet+0.1429, runProp_offsetSesh(isession,ifet),'ro');
        end
    end
    end

end


