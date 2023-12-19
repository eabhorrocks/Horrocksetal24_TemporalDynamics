%% Analyse PSTH metrics by cell-type

%% onset and offset features by cellType

plotFlag=false;

for itype = 1:3

%%%%%%%%%%%% onset features %%%%%%%%%

statVals = stat.onsetCharVec(statGoodVals & stat.cellType==itype);
runVals = run.onsetCharVec(runGoodVals & run.cellType==itype);
statSessions = stat.session(statGoodVals & stat.cellType==itype);
runSessions = run.session(runGoodVals & run.cellType==itype);

cellType(itype).nStat = numel(statVals);
cellType(itype).nRun = numel(runVals);


statProp_onset = [];
runProp_onset = [];

for ifet = 1:5
    statProp_onset(ifet) = prop(statVals==ifet);
    runProp_onset(ifet) = prop(runVals==ifet);
end

cellType(itype).statProp_onset = statProp_onset;
cellType(itype).runProp_onset = runProp_onset;


for isession = 1:5
    tempStatVals = statVals(statSessions==isession);
    tempRunVals = runVals(runSessions==isession);
    for ifet = 1:5
        statProp_s_onset(ifet,isession) = prop(tempStatVals==ifet);
        runProp_s_onset(ifet,isession) = prop(tempRunVals==ifet);
    end
end


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

    cellType(itype).pVal_onset(ifet) = glme.Coefficients.pValue(2);
end



if plotFlag
figure, hold on
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
end



%%%%%%%%%%%% offset features %%%%%%%%%
statVals = stat.offsetCharVec(statGoodVals & stat.cellType==itype);
runVals = run.offsetCharVec(runGoodVals & run.cellType==itype);
statSessions = stat.session(statGoodVals & stat.cellType==itype);
runSessions = run.session(runGoodVals & run.cellType==itype);

statProp_offset = [];
runProp_offset = [];

for ifet = 1:5
    statProp_offset(ifet) = prop(statVals==ifet);
    runProp_offset(ifet) = prop(runVals==ifet);
    
end

cellType(itype).statProp_offset = statProp_offset;
cellType(itype).runProp_offset = runProp_offset;

for isession = 1:5
    tempStatVals = statVals(statSessions==isession);
    tempRunVals = runVals(runSessions==isession);
    for ifet = 1:5
        statProp_s_offset(ifet,isession) = prop(tempStatVals==ifet);
        runProp_s_offset(ifet,isession) = prop(tempRunVals==ifet);
    end
end

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

    cellType(itype).pVal_offset(ifet) = glme.Coefficients.pValue(2);
end


if plotFlag
figure, hold on

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
end


end

%% sustainedness index by cellType

for itype = 1:3

metric_name = 'susIndexVec';
statVals = stat.(metric_name);
runVals = run.(metric_name);

valsVec = cat(1,statVals(bothGoodVals & stat.cellType==itype), runVals(bothGoodVals & run.cellType==itype));
unitVec = categorical(cat(2,1:sum(bothGoodVals & stat.cellType==itype), 1:sum(bothGoodVals & run.cellType==itype))');
stateVec = categorical(cat(1,repelem(1,sum(bothGoodVals & stat.cellType==itype),1), repelem(2,sum(bothGoodVals & run.cellType==itype),1)));
seshVec = categorical(cat(1,stat.session(bothGoodVals & stat.cellType==itype), run.session(bothGoodVals & run.cellType==itype)));

tbl = table(valsVec,unitVec,stateVec,seshVec,...
    'VariableNames',{'vals','unit','state','sesh'});

    f = 'vals ~ state + (1|unit) + (1|sesh)';
      
    lme = fitlme(tbl, f, 'DummyVarCoding', 'reference');

    cellType(itype).susIndexP = lme.Coefficients.pValue(2);

    cellType(itype).susIndexQuant_stat = quantile(statVals(bothGoodVals & stat.cellType==itype),3);
    cellType(itype).susIndexQuant_run = quantile(runVals(bothGoodVals & run.cellType==itype),3);
    
    cellType(itype).nSusIdx = sum(bothGoodVals & stat.cellType==itype);

end


cellTypes = {'pyramidal', 'narrow', 'wide'};
cellCols = {'r', 'b', 'c'};


%% plot sustainedness index by cell-type
cellTypes = {'pyramidal', 'narrow', 'wide'};
cellCols = {'r', 'b', 'c'};



metric_name = 'susIndexVec';
binEdges = 0:0.05:1;

statVals = stat.(metric_name);
runVals = run.(metric_name);
figure
for itype=1:3
    thisStat = statVals(bothGoodVals & stat.cellType==itype);
    thisRun  = runVals(bothGoodVals & run.cellType==itype);

    subplot(1,3,itype), hold on
    plot(thisStat, thisRun, 'o', 'MarkerSize', 4, 'MarkerFaceColor', cellCols{itype}, 'MarkerEdgeColor', 'none')
    plot([0 1], [0 1], 'k')
    xlim([0 1]), ylim([0 1])
    title(cellTypes{itype})
    defaultAxesProperties(gca,true)
end


