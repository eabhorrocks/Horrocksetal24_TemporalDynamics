%% Analyse PSTH metrics by speed

%% onset and offset features by speed

plotFlag=false;

for ispeed = 1:6

%%%%%%%%%%%% onset features %%%%%%%%%

statVals = stat.onsetCharVec(statGoodVals & stat.speedVec==ispeed);
runVals = run.onsetCharVec(runGoodVals & run.speedVec==ispeed);
statSessions = stat.session(statGoodVals & stat.speedVec==ispeed);
runSessions = run.session(runGoodVals & run.speedVec==ispeed);

speed(ispeed).nStat = numel(statVals);
speed(ispeed).nRun = numel(runVals);


statProp_onset = [];
runProp_onset = [];

for ifet = 1:5
    statProp_onset(ifet) = prop(statVals==ifet);
    runProp_onset(ifet) = prop(runVals==ifet);
end

speed(ispeed).statProp_onset = statProp_onset;
speed(ispeed).runProp_onset = runProp_onset;


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

    speed(ispeed).pVal_onset(ifet) = glme.Coefficients.pValue(2);
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
statVals = stat.offsetCharVec(statGoodVals & stat.speedVec==ispeed);
runVals = run.offsetCharVec(runGoodVals & run.speedVec==ispeed);
statSessions = stat.session(statGoodVals & stat.speedVec==ispeed);
runSessions = run.session(runGoodVals & run.speedVec==ispeed);

statProp_offset = [];
runProp_offset = [];

for ifet = 1:5
    statProp_offset(ifet) = prop(statVals==ifet);
    runProp_offset(ifet) = prop(runVals==ifet);
    
end

speed(ispeed).statProp_offset = statProp_offset;
speed(ispeed).runProp_offset = runProp_offset;

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

    speed(ispeed).pVal_offset(ifet) = glme.Coefficients.pValue(2);
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

%% sustainedness index by speed

for ispeed = 1:6

metric_name = 'susIndexVec';
statVals = stat.(metric_name);
runVals = run.(metric_name);

valsVec = cat(1,statVals(bothGoodVals & stat.speedVec==ispeed), runVals(bothGoodVals & run.speedVec==ispeed));
unitVec = categorical(cat(2,1:sum(bothGoodVals & stat.speedVec==ispeed), 1:sum(bothGoodVals & run.speedVec==ispeed))');
stateVec = categorical(cat(1,repelem(1,sum(bothGoodVals & stat.speedVec==ispeed),1), repelem(2,sum(bothGoodVals & run.speedVec==ispeed),1)));
seshVec = categorical(cat(1,stat.session(bothGoodVals & stat.speedVec==ispeed), run.session(bothGoodVals & run.speedVec==ispeed)));

tbl = table(valsVec,unitVec,stateVec,seshVec,...
    'VariableNames',{'vals','unit','state','sesh'});

    f = 'vals ~ state + (1|unit) + (1|sesh)';
      
    lme = fitlme(tbl, f, 'DummyVarCoding', 'reference');

    speed(ispeed).susIndexP = lme.Coefficients.pValue(2);

    speed(ispeed).susIndexQuant_stat = quantile(statVals(bothGoodVals & stat.speedVec==ispeed),3);
    speed(ispeed).susIndexQuant_run = quantile(runVals(bothGoodVals & run.speedVec==ispeed),3);
    
    speed(ispeed).nSusIdx = sum(bothGoodVals & stat.speedVec==ispeed);

end


