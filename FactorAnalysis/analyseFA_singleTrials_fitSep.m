%% analyse FA single trials (fit separattely)

%% load dataa
sessionTags = {'M22027', '20220517';...
    'M22029', '20220607';...
    'M22031', '20220524';
    'M22032', '20220621';
    'M22033', '20220706'};

dataDir = 'C:\Users\edward.horrocks\Documents\Code\V1Dynamics\Data\basic_111022';

s = struct;

for isession = 1:5
    tic
    fname = [sessionTags{isession,1},'_', sessionTags{isession,2},'_FA_fitSep.mat'];

    load(fullfile(dataDir,fname))

    s(isession).session = session;
end

%% recreate cond vaara
for isesh=1:5
for icond = 1:6
    s(isesh).session.s.cond(icond) = s(isesh).session.stat.s.cond(icond);
end
for icond = 7:12
    s(isesh).session.s.cond(icond) = s(isesh).session.run.s.cond(icond-6);
end
end


%% Local tangling


for isession = 1:5

    cond = s(isession).session.s.cond;
    nTrials = size(cond(1).catData,3);

    clear sesh Qstat Qrun

    qOpt_stat = s(isesh).session.stat.s.qOpt;
    qOpt_run = s(isesh).session.run.s.qOpt;
    qOpt = min([qOpt_stat, qOpt_run]);
    deltaT = 10;
    epsVal = 0.1;
    binSize = 10;
    windowSize = 200/binSize;
    prctileVal= 90;

    for icond = 1:12
        for itrial = 1:nTrials
            itime = 1:200;
            thisTraj =  cond(icond).catData(:,1:qOpt,itrial);
            Q=[];

            for it = 2:numel(itime)
                for it2 = 2:numel(itime)

                    xt = thisTraj(itime(it),:);
                    xt2 = thisTraj(itime(it2),:);

                    xderiv_t = (thisTraj(itime(it),:) - thisTraj(itime(it-1),:))/deltaT;
                    xderiv_t2 = (thisTraj(itime(it2),:) - thisTraj(itime(it2-1),:))/deltaT;

                    Q(it,it2) = (norm(xderiv_t-xderiv_t2)^2) / (norm(xt-xt2)^2 + epsVal);

                end
            end

            cond(icond).trial(itrial).Q = Q;

            for itime = 1:180
                tq = Q(itime:itime+windowSize, itime:itime+windowSize);
                cond(icond).trial(itrial).localTangling(itime) = ...
                    prctile(tq(1+windowSize/2,:),prctileVal,'all');

            end
        end
    end
    s(isession).session.cond = cond;
end



%% plot average of single-trial tangling

allStatTangling = [];
allRunTangling = [];


for isession = 1:5
    nTrials = size(s(isession).session.cond(1).catData,3);
    for icond = 1:6
        for itrial = 1:nTrials
            allStatTangling=cat(1,allStatTangling,...
                s(isession).session.cond(icond).trial(itrial).localTangling);

            allRunTangling=cat(1,allRunTangling,...
                s(isession).session.cond(icond+6).trial(itrial).localTangling);            
        end
    end
end

figure, hold on
shadedErrorBar(-95:10:1695, mean(allStatTangling,1), sem(allStatTangling,1))
hold on
shadedErrorBar(-95:10:1695, mean(allRunTangling,1), sem(allRunTangling,1), 'lineProps','r')


%  plot speed-wise local tangling
figure

for ispeed = 1:6
    allStat =[];
    allRun=[];
    for isesh = 1:5
    allStat=cat(1,allStat, cat(1,s(isesh).session.cond(ispeed).trial.localTangling));
    allRun=cat(1,allRun, cat(1,s(isesh).session.cond(ispeed+6).trial.localTangling));
    end
    subplot(2,3,ispeed), hold on
    shadedErrorBar(-95:10:1695, mean(allStat,1), sem(allStat,1))
    shadedErrorBar(-95:10:1695, mean(allRun,1), sem(allRun,1),...
        'lineProps','r')

end



%% stats for single trial tangling

subjVec=[];
speedVec=[];
allStatMeanTangling = [];
allRunMeanTangling = [];

for isession = 1:5
    for icond = 1:6
        
        allStatMeanTangling = cat(1,allStatMeanTangling, mean(cat(1,(s(isession).session.cond(icond).trial.localTangling)),1));
        allRunMeanTangling = cat(1,allRunMeanTangling, mean(cat(1,(s(isession).session.cond(icond+6).trial.localTangling)),1));

            subjVec = cat(1,subjVec,repelem(isession,180,1));
            speedVec = cat(1,speedVec,repelem(icond,180,1));

    end
end

% % RM-anova for tangling
% 
% time = 1:180;
% statVals = allStatMeanTangling';
% runVals = allRunMeanTangling';
% 
% allVals = cat(1,statVals(:),runVals(:));
% 
% timeVec = categorical(repmat(1:180,1,size(allVals,1)/180))';
% speedVec = categorical(cat(1,speedVec,speedVec));
% subjVec = categorical(cat(1,subjVec,subjVec));
% stateVec = categorical(repelem(1:2,1,numel(statVals))');
% 
% % [p,tbl,stats,terms] = anovan(allVals,{timeVec,stateVec,subjVec},'model',2,'random',3,'varnames',{'Time','State','Subj'});
% 
% 
% [p,tbl,stats,terms] = anovan(allVals,{timeVec,stateVec,subjVec,speedVec},'model','full','varnames',{'Time','State','Subj','Speed'});


%% get speed and acceleration of trajectories


for isession = 1:5
    qOpt_stat = s(isession).session.stat.s.qOpt;
    qOpt_run = s(isession).session.run.s.qOpt;
    qOpt = min([qOpt_stat, qOpt_run]);  
    nTrials = size(s(isession).session.s.cond(1).catData,3);
    nFactors = qOpt;
    weights = repelem(1, 1, qOpt);

    for icond = 1:12
        for itrial = 1:nTrials
            thisTraj =  s(isession).session.s.cond(icond).catData(:,1:qOpt,itrial);

            for itime = 1:200-1 % speed
                pt1 = thisTraj(itime,:);
                pt2 = thisTraj(itime+1,:);
                cond(icond).trial(itrial).speed(itime) = ...
                    (vecnorm(pt2-pt1)*(1000/10)) / nFactors;

            end
            cond(icond).trial(itrial).speedAcc = diff(cond(icond).trial(itrial).speed);

        cond(icond).trial(itrial).maxSpeedOnset = max(cond(icond).trial(itrial).speed(1:61));
        cond(icond).trial(itrial).maxSpeedOffset = max(cond(icond).trial(itrial).speed(100:161));


        end
    end

    s(isession).session.cond = cond;
end




%% plot trial-average of trajectory speed and acceleration

allStatSpeed = [];
allRunSpeed = [];

allStatVelAcc = [];
allRunVelAcc = [];

for isession = 1:5
    nTrials = size(s(isession).session.cond(1).catData,3);
    for icond = 1:6
        for itrial = 1:nTrials
            allStatSpeed=cat(1,allStatSpeed,...
                s(isession).session.cond(icond).trial(itrial).speed);

            allRunSpeed=cat(1,allRunSpeed,...
                s(isession).session.cond(icond+6).trial(itrial).speed);

            allStatVelAcc=cat(1,allStatVelAcc,...
                s(isession).session.cond(icond).trial(itrial).speedAcc);

            allRunVelAcc=cat(1,allRunVelAcc,...
                s(isession).session.cond(icond+6).trial(itrial).speedAcc);
        end
    end
end

figure, hold on
shadedErrorBar(-190:10:1790, mean(allStatSpeed,1), sem(allStatSpeed,1))
hold on
shadedErrorBar(-190:10:1790, mean(allRunSpeed,1), sem(allRunSpeed,1), 'lineProps','r')
xlim([-200 1800])
ax=gca; ax.XTick = -200:200:1800;
defaultAxesProperties(gca, true)
ylabel('Speed')

figure, hold on
shadedErrorBar(-180:10:1790, mean(allStatVelAcc.*100,1), sem(allStatVelAcc.*100,1))
hold on
shadedErrorBar(-180:10:1790, mean(allRunVelAcc.*100,1), sem(allRunVelAcc.*100,1), 'lineProps','r')
xlim([-200 1800])
ax=gca; ax.XTick = -200:200:1800;
defaultAxesProperties(gca, true)
ylabel('Speed acceleration')


%% stats for speeds of trajectories

speedVec=[];
seshVec=[];
onsetMaxSpdVec=[];
offsetMaxSpdVec=[];
stateVec=[];

for isesh = 1:5
    for ispeed = 1:6
        for itrial = 1:numel(s(isesh).session.cond(ispeed).trial)
        speedVec=cat(1,speedVec,ispeed);
        seshVec=cat(1,seshVec,isesh);
        stateVec=cat(1,stateVec,1);
        onsetMaxSpdVec=cat(1,onsetMaxSpdVec,s(isesh).session.cond(ispeed).trial(itrial).maxSpeedOnset);
        offsetMaxSpdVec=cat(1,offsetMaxSpdVec,s(isesh).session.cond(ispeed).trial(itrial).maxSpeedOffset);
        end
    end
end



for isesh = 1:5
    for ispeed = 7:12
        for itrial = 1:numel(s(isesh).session.cond(ispeed).trial)
        speedVec=cat(1,speedVec,ispeed-6);
        seshVec=cat(1,seshVec,isesh);
        stateVec=cat(1,stateVec,2);
        onsetMaxSpdVec=cat(1,onsetMaxSpdVec,s(isesh).session.cond(ispeed).trial(itrial).maxSpeedOnset);
        offsetMaxSpdVec=cat(1,offsetMaxSpdVec,s(isesh).session.cond(ispeed).trial(itrial).maxSpeedOffset);
        end
    end
end


valsVec = offsetMaxSpdVec;

tbl = table(valsVec,speedVec,stateVec,seshVec,...
    'VariableNames',{'vals','speed','state','sesh'});

    f = 'vals ~ state + (1|speed) + (1|sesh)';
      
    lme = fitlme(tbl, f, 'DummyVarCoding', 'reference')


    [mean(valsVec(1:30)), sem(valsVec(1:30)); mean(valsVec(31:60)), sem(valsVec(31:60))]
