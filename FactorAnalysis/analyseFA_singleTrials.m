%% analyse FA single-trials

sessionTags = {'M22027', '20220517';...
    'M22029', '20220607';...
    'M22031', '20220524';
    'M22032', '20220621';
    'M22033', '20220706'};

dataDir = 'C:\Users\edward.horrocks\Documents\Code\V1Dynamics\Data\basic_111022';

s = struct;

for isession = 1:5
    tic
    fname = [sessionTags{isession,1},'_', sessionTags{isession,2},'_FAsmooth_Mar23.mat'];

    load(fullfile(dataDir,fname))

    s(isession).session = session;
end

speedcols=inferno(6);

%% plot decoding results

% plot average decoding performance
bmp = -177:10:1775; 

allStat = [];
allRun = [];

for isession = 1:5
    allStat = cat(1,allStat, s(isession).session.stat.meanPerf);
    allRun = cat(1,allRun, s(isession).session.run.meanPerf);
end

figure,  hold on
plot([-0.2 1.8], [1/6 1/6], 'k:')
shadedErrorBar(bmp, mean(allStat,1), sem(allStat,1))
shadedErrorBar(bmp, mean(allRun,1), sem(allRun,1),'lineProps', 'r')
ax = gca; ax.XTick = -200:200:1800;
xlim([-200 1800])
defaultAxesProperties(gca,true)


%% plot single trial examples, for multiple speeds

speedcols = inferno(6);

nTrials = 10;
isession = 2;
dims = [1 3];

speeds2plot = [2 4];
time2plot = 21:50;



figure
% stat
subplot(121), hold on
trials2plot = randsample(size(s(isession).session.cond(ispeed1).catData,3),nTrials);

for itr = 1:nTrials

    itrial = trials2plot(itr);

    for ispeed = speeds2plot
        % ispeed 1
        plot(s(isession).session.cond(ispeed).catData(time2plot,dims(1),itrial),...
            s(isession).session.cond(ispeed).catData(time2plot,dims(2),itrial),'Color', speedcols(ispeed,:))

        plot(s(isession).session.cond(ispeed).catData(time2plot(1),dims(1),itrial),...
            s(isession).session.cond(ispeed).catData(time2plot(1),dims(2),itrial),'Color', speedcols(ispeed,:),'Marker','v')

        plot(s(isession).session.cond(ispeed).catData(time2plot(end),dims(1),itrial),...
            s(isession).session.cond(ispeed).catData(time2plot(end),dims(2),itrial),'Color', speedcols(ispeed,:),'Marker','o')
    end
end

% run
subplot(122), hold on

for itr = 1:nTrials

    itrial = trials2plot(itr);

    for ispeed = speeds2plot
        plot(s(isession).session.cond(ispeed+6).catData(time2plot,dims(1),itrial),...
            s(isession).session.cond(ispeed+6).catData(time2plot,dims(2),itrial),'Color', speedcols(ispeed,:))

        plot(s(isession).session.cond(ispeed+6).catData(time2plot(1),dims(1),itrial),...
            s(isession).session.cond(ispeed+6).catData(time2plot(1),dims(2),itrial),'Color', speedcols(ispeed,:),'Marker','v')

        plot(s(isession).session.cond(ispeed+6).catData(time2plot(end),dims(1),itrial),...
            s(isession).session.cond(ispeed+6).catData(time2plot(end),dims(2),itrial),'Color', speedcols(ispeed,:),'Marker','o')
    end
end

%% 3d plots of single trials for 2 speeds


speedcols = inferno(6);

nTrials = 10;
isession = 2;
dims = [1 3 4];

speeds2plot = [2 4];
time2plot = 21:70;



figure
% stat
subplot(121), hold on
% trials2plot = randsample(size(s(isession).session.cond(ispeed1).catData,3),nTrials);

for itr = 1:nTrials

    itrial = trials2plot(itr);

    for ispeed = speeds2plot
        % ispeed 1
        plot3(s(isession).session.cond(ispeed).catData(time2plot,dims(1),itrial),...
            s(isession).session.cond(ispeed).catData(time2plot,dims(2),itrial),...
            s(isession).session.cond(ispeed).catData(time2plot,dims(3),itrial),'Color', speedcols(ispeed,:))

        plot3(s(isession).session.cond(ispeed).catData(time2plot(1),dims(1),itrial),...
            s(isession).session.cond(ispeed).catData(time2plot(1),dims(2),itrial),...
            s(isession).session.cond(ispeed).catData(time2plot(1),dims(3),itrial),'Color', speedcols(ispeed,:),'Marker','v')

        plot3(s(isession).session.cond(ispeed).catData(time2plot(end),dims(1),itrial),...
            s(isession).session.cond(ispeed).catData(time2plot(end),dims(2),itrial),...
            s(isession).session.cond(ispeed).catData(time2plot(end),dims(3),itrial),'Color', speedcols(ispeed,:),'Marker','o')
    end
end

% run
subplot(122), hold on

for itr = 1:nTrials

    itrial = trials2plot(itr);

    for ispeed = speeds2plot
        plot3(s(isession).session.cond(ispeed+6).catData(time2plot,dims(1),itrial),...
            s(isession).session.cond(ispeed+6).catData(time2plot,dims(2),itrial),...
            s(isession).session.cond(ispeed+6).catData(time2plot,dims(3),itrial),'Color', speedcols(ispeed,:))

        plot3(s(isession).session.cond(ispeed+6).catData(time2plot(1),dims(1),itrial),...
            s(isession).session.cond(ispeed+6).catData(time2plot(1),dims(2),itrial),...
            s(isession).session.cond(ispeed+6).catData(time2plot(1),dims(3),itrial),'Color', speedcols(ispeed,:),'Marker','v')

        plot3(s(isession).session.cond(ispeed+6).catData(time2plot(end),dims(1),itrial),...
            s(isession).session.cond(ispeed+6).catData(time2plot(end),dims(2),itrial),...
            s(isession).session.cond(ispeed+6).catData(time2plot(end),dims(3),itrial),'Color', speedcols(ispeed,:),'Marker','o')
    end
end



%% Local tangling


for isession = 1:5

    cond = s(isession).session.cond;
    nTrials = size(cond(1).catData,3);

    clear sesh Qstat Qrun

    qOpt = s(isession).session.s.qOpt;
    deltaT = 10;
    epsVal = 0.1;
    binSize = 10;
    windowSize = 200/binSize;
    prctileVal= 90;

    for icond = 1:12
        for itrial = 1:nTrials
            itime = 1:200;
            thisTraj =  s(isession).session.cond(icond).catData(:,1:qOpt,itrial);
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
allStat =[];
allRun=[];
for icond = 1:6
    allStat=cat(1,allStat, cat(1,cond(icond).trial.localTangling));
    allRun=cat(1,allRun, cat(1,cond(icond+6).trial.localTangling));
    subplot(2,3,icond), hold on
    shadedErrorBar(-95:10:1695, mean(cat(1,cond(icond).trial.localTangling),1), sem(cat(1,cond(icond).trial.localTangling),1))
    shadedErrorBar(-95:10:1695, mean(cat(1,cond(icond+6).trial.localTangling),1), sem(cat(1,cond(icond+6).trial.localTangling),1),...
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

% RM-anova for tangling

time = 1:180;
statVals = allStatMeanTangling';
runVals = allRunMeanTangling';

allVals = cat(1,statVals(:),runVals(:));

timeVec = categorical(repmat(1:180,1,size(allVals,1)/180))';
speedVec = categorical(cat(1,speedVec,speedVec));
subjVec = categorical(cat(1,subjVec,subjVec));
stateVec = categorical(repelem(1:2,1,numel(statVals))');

% [p,tbl,stats,terms] = anovan(allVals,{timeVec,stateVec,subjVec},'model',2,'random',3,'varnames',{'Time','State','Subj'});


[p,tbl,stats,terms] = anovan(allVals,{timeVec,stateVec,subjVec,speedVec},'model','full','varnames',{'Time','State','Subj','Speed'});


%% get speed and acceleration of trajectories


for isession = 1:5
    qOpt = s(isession).session.s.qOpt;
    cond= s(isession).session.cond;
    nTrials = size(cond(1).catData,3);
    nFactors = size(s(1).session.s.loadings,2);
    weights = repelem(1, 1, qOpt);

    for icond = 1:12
        for itrial = 1:nTrials
            thisTraj =  s(isession).session.cond(icond).catData(:,1:qOpt,itrial);

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





%% plot some single trial speed profiles

isession = 2;
icond = 2;

figure, hold on
for itrial = 1:numel(s(isession).session.cond(icond).trial)
    plot(s(isession).session.cond(icond).trial(itrial).speed,'k')
end
for itrial = 1:numel(s(isession).session.cond(icond).trial)
    plot(s(isession).session.cond(icond+6).trial(itrial).speed,'r')
end



%% get geometric measures of trials


% set parameters
bin_N = 1; % bin every N elements
binSize = 10; % bin size in ms
binVector = -200:binSize*bin_N:(1800-binSize*bin_N);
timeBinVector = round(binVector,2); % time associated with each index (1st edge of bin)
weights = repelem(1,qOpt,1); weights = weights(:)'; % weights are frac. of shared variance explained by each component
nWins = 10; % number of consecutive windows required to reach steady state
initialSSIdx = 1:find(timeBinVector==0)-1;
finalSSIdx = find(timeBinVector==500):find(timeBinVector==1000)-1;
nTrials = size(cond(1).catData,3);

for isession = 1:5

    cond = s(isession).session.cond;
    nTrials = size(cond(1).catData,3);


    for icond = 1:12
        for itrial = 1:nTrials
            neuralTrajectory = cond(icond).catData(:,1:qOpt,itrial);
            cond(icond).trial(itrial).geo_onset = calcPopulationGeoMetrics(neuralTrajectory, weights, timeBinVector, initialSSIdx, finalSSIdx, nWins);
        end
    end

    s(isession).session.cond = cond;

end


% offset trajectories
initialSSIdx = find(timeBinVector==500):find(timeBinVector==1000)-1;
finalSSIdx = find(timeBinVector==1500):find(timeBinVector==1800-binSize*bin_N);

for isession = 1:5

    cond = s(isession).session.cond;
    nTrials = size(cond(1).catData,3);


    for icond = 1:12
        for itrial = 1:nTrials
            neuralTrajectory = cond(icond).catData(:,1:qOpt,itrial);
            cond(icond).trial(itrial).geo_offset = calcPopulationGeoMetrics(neuralTrajectory, weights, timeBinVector, initialSSIdx, finalSSIdx, nWins);
        end
    end

    s(isession).session.cond = cond;

end


for isession = 1:5
    nTrials = size(s(isession).session.cond(1).catData,3);
    for icond = 1:12
        tempVals = [];
        for itrial = 1:nTrials
            tempVals = cat(1,tempVals,s(isession).session.cond(icond).trial(itrial).geo_onset.distanceRatio);
        end
        s(isession).session.cond(icond).meanDRon = nanmean(tempVals);
    end
end

for isession = 1:5
    nTrials = size(s(isession).session.cond(1).catData,3);
    for icond = 1:12
        tempVals = [];
        for itrial = 1:nTrials
            tempVals = cat(1,tempVals,s(isession).session.cond(icond).trial(itrial).geo_offset.distanceRatio);
        end
        s(isession).session.cond(icond).meanDRoff = nanmean(tempVals);
    end
end

figure,
subplot(121), hold on
plot(log2([1 10]), log2([1 10]), 'k:')
for icond = 1:6
    for isession = 1:5
        plot(log2(s(isession).session.cond(icond).meanDRon), log2(s(isession).session.cond(icond+6).meanDRon),...
            'o', 'MarkerEdgeColor','none', 'MarkerFaceColor',speedcols(icond,:))
    end
end
xlim(log2([1 10]))
ylim(log2([1 10]))
ax = gca; ax.XTick = log2([1 2 4 8]); ax.YTick = ax.XTick;
ax.XTickLabel = [1 2 4 8]; ax.YTickLabel = ax.XTickLabel;


subplot(122), hold on
plot(log2([1 10]), log2([1 10]), 'k:')
for icond = 1:6
    for isession = 1:5
        plot(log2(s(isession).session.cond(icond).meanDRoff), log2(s(isession).session.cond(icond+6).meanDRoff),...
            'o', 'MarkerEdgeColor','none', 'MarkerFaceColor',speedcols(icond,:))
    end
end
xlim(log2([1 10]))
ylim(log2([1 10]))
ax = gca; ax.XTick = log2([1 2 4 8]); ax.YTick = ax.XTick;
ax.XTickLabel = [1 2 4 8]; ax.YTickLabel = ax.XTickLabel;





%% stats on distance ratios

statVals = [];
runVals = [];

speedVec = [];
subjVec = [];

for isession = 1:5
    for ispeed = 1:6
            statVals = cat(1,statVals, s(isession).session.cond(ispeed).meanDRoff);
            runVals = cat(1,runVals, s(isession).session.cond(ispeed+6).meanDRoff);

            speedVec = cat(1,speedVec,ispeed);
            subjVec = cat(1,subjVec, isession);
    end
end



valsVec = log2(cat(1,statVals(:),runVals(:)));

sessionVec = cat(1,subjVec,subjVec);
speedVec = cat(1,speedVec,speedVec);
stateVec = categorical(cat(1,repelem(1,numel(statVals),1), repelem(2,numel(runVals),1)));

tbl = table(valsVec,speedVec,stateVec,sessionVec,...
    'VariableNames',{'vals','speed','state','sesh'});

    f = 'vals ~ state + (1|speed) + (1|sesh)';
      
    lme = fitlme(tbl, f, 'DummyVarCoding', 'reference')


    [mean(statVals(:)), sem(statVals(:)); mean(runVals(:)), sem(runVals(:))]


%% angle of approach analysis


tic
for isession = 1:5
    nTrials = size(s(isession).session.cond(1).catData,3);


    for icond = 1:12
        bin_N = 1; % bin every N elements
        binSize = 10; % bin size in ms
        binVector = -200:binSize*bin_N:(1800-binSize*bin_N);
        timeBinVector = round(binVector,2); % time associated with each index (1st edge of bin)
        nWins = 10; % number of consecutive windows required to reach steady state

        % set parameters
        % weights = s.propSharedVariance(:)'; % weights are frac. of shared variance explained by each component
        weights = repelem(1, 1, s(isession).session.s.qOpt);

        % onset trajectories
        initialSSIdx = 1:find(timeBinVector==0)-1;
        finalSSIdx = find(timeBinVector==500):find(timeBinVector==1000)-1;
        finalBaseIdx = find(timeBinVector==1500):find(timeBinVector==1790)-1;

        for itrial = 1:nTrials

            traj = s(isession).session.cond(icond).catData(:,1:s(isession).session.s.qOpt,itrial);


            ss_mean = mean(traj(finalSSIdx,:),1);
            ini_mean = mean(traj(initialSSIdx,:),1);
            final_mean = mean(traj(finalBaseIdx,:),1);

            refAxis = ss_mean-ini_mean;
            refAxisOffset = final_mean-ss_mean;


            % onset
            angle=[];
            angleSub=[];
            angleDir=[];
            v1full=[];
            for itime = 1:200
                v1full(itime,:) = ss_mean-traj(itime,:);
                v1 = ss_mean-traj(itime,:); v1=v1(:);
                v2 = refAxis; v2=v2(:);
                angle(itime) = dot(v1 / norm(v1), v2 / norm(v2));
                distFromSS(itime) = sum(abs(traj(itime,:)-ini_mean));
                %      angleSub(itime) = subspace(v1,v2);


            end

            % offset
            angleOffset=[];
            distFromSSOffset=[];

            for itime = 1:200
                v1full(itime,:) = final_mean-traj(itime,:);
                v1 = final_mean-traj(itime,:); v1=v1(:);
                v2 = refAxisOffset; v2=v2(:);
                angleOffset(itime) = dot(v1 / norm(v1), v2 / norm(v2));
                distFromSSOffset(itime) = sum(abs(traj(itime,:)-ini_mean));
                %      angleSub(itime) = subspace(v1,v2);
            end

            s(isession).session.cond(icond).trial(itrial).angle = acosd(angle);%mod(acosd(angle).*sign(angle),360);
            s(isession).session.cond(icond).trial(itrial).angleOffset = acosd(angleOffset);%mod(acosd(angle).*sign(angle),360);

            s(isession).session.cond(icond).trial(itrial).distFromStim = distFromSS;
            s(isession).session.cond(icond).trial(itrial).distFromStim = distFromSSOffset;

            s(isession).session.cond(icond).trial(itrial).sumDiffAngle = sum(abs(diff(s(isession).session.cond(icond).trial(itrial).angle)));
            s(isession).session.cond(icond).trial(itrial).sumDiffAngleOffset = sum(abs(diff(s(isession).session.cond(icond).trial(itrial).angleOffset)));

            s(isession).session.cond(icond).trial(itrial).v1full = v1full;
            % s(isession).session.cond(icond).angleSub = angleSub;

        end


    end
end

toc

%% plot angle over time
allStat=[];
allRun=[];

for isession = 1:5
    nTrials = size(s(isession).session.cond(1).catData,3);
    for ispeed = 1:6
        for itrial = 1:nTrials
            allStat=cat(1,allStat, s(isession).session.cond(ispeed).trial(itrial).angle);
            allRun=cat(1,allRun, s(isession).session.cond(ispeed+6).trial(itrial).angle);
        end
    end
end

figure
shadedErrorBar(-195:10:1795,mean(allStat,1), sem(allStat,1), 'lineProps','k')
shadedErrorBar(-195:10:1795,mean(allRun,1), sem(allRun,1), 'lineProps','r')
defaultAxesProperties(gca, true)
title('Onset')

figure
subplot(121)
imagesc(allStat)
colorbar
caxis([0 160])
subplot(122)
imagesc(allRun)
colorbar
caxis([0 160])
title('Onset')


allVals = allRun(:,21:120);
maxStretch =5; % bins

nPSTH = size(allVals,1);
dm = nan([nPSTH]);

% nest loop version
tic
for ipsth1 = 1:nPSTH
    psth_temp = allVals(ipsth1,:);
    parfor ipsth2 = ipsth1:nPSTH
        dm_temp(:,ipsth2) = dtw(psth_temp,allVals(ipsth2,:), maxStretch);
    end
    dm(ipsth1,:) = dm_temp;
    dm_temp = [];
end
toc
full_dm = triu(dm)+triu(dm,1)';

% get dendrogram
Z_sub = linkage(full_dm, 'average');
% get optimal leaf ordering
tic
leafOrder = optimalleaforder(Z_sub,full_dm);
toc
[denHandle, leafnode_idx, outpermNodeOrder] = dendrogram(Z_sub, nPSTH, 'reorder', leafOrder, 'orientation', 'right');

newPSTHorder = [];
speedVecReordered = [];

tic
for ileaf = 1:numel(outpermNodeOrder)
    idx = find(leafnode_idx==outpermNodeOrder(ileaf));
    newPSTHorder = vertcat(newPSTHorder, idx(:));
end
toc

allStatReordered = allStat(newPSTHorder,:);
allRunReordered = allRun(newPSTHorder,:);

figure
subplot(121)
imagesc(allStatReordered);
caxis([0 160])
colorbar
subplot(122)
imagesc(allRunReordered)
caxis([0 160])
colorbar
colormap(viridis)

%% offset

allStat=[];
allRun=[];

for isession = 1:5
    nTrials = size(s(isession).session.cond(1).catData,3);
    for ispeed = 1:6
        for itrial = 1:nTrials

            allStat=cat(1,allStat, s(isession).session.cond(ispeed).trial(itrial).angleOffset);
            allRun=cat(1,allRun, s(isession).session.cond(ispeed+6).trial(itrial).angleOffset);
        end
    end
end

figure
shadedErrorBar(-195:10:1795,mean(allStat,1), sem(allStat,1), 'lineProps','k')
shadedErrorBar(-195:10:1795,mean(allRun,1), sem(allRun,1), 'lineProps','r')
ylim([0 140])
defaultAxesProperties(gca, true)
title('Offset')


figure
subplot(121)
imagesc(allStat)
colorbar
caxis([0 160])

subplot(122)
imagesc(allRun)
colorbar
caxis([0 160])
title('Offset')

allVals = allRun(:,121:end);
maxStretch = 5; % bins

nPSTH = size(allVals,1);
dm = nan([nPSTH]);

% nest loop version
tic
for ipsth1 = 1:nPSTH
    psth_temp = allVals(ipsth1,:);
    parfor ipsth2 = ipsth1:nPSTH
        dm_temp(:,ipsth2) = dtw(psth_temp,allVals(ipsth2,:), maxStretch);
    end
    dm(ipsth1,:) = dm_temp;
    dm_temp = [];
end
toc
full_dm = triu(dm)+triu(dm,1)';

% get dendrogram
Z_sub = linkage(full_dm, 'average');
% get optimal leaf ordering
tic
leafOrder = optimalleaforder(Z_sub,full_dm);
toc
[denHandle, leafnode_idx, outpermNodeOrder] = dendrogram(Z_sub, nPSTH, 'reorder', leafOrder, 'orientation', 'right');

newPSTHorder = [];
speedVecReordered = [];

tic
for ileaf = 1:numel(outpermNodeOrder)
    idx = find(leafnode_idx==outpermNodeOrder(ileaf));
    newPSTHorder = vertcat(newPSTHorder, idx(:));
end
toc

allStatReordered = allStat(newPSTHorder,:);





allRunReordered = allRun(newPSTHorder,:);

figure
subplot(121)
imagesc(allStatReordered);
caxis([0 160])
colorbar
subplot(122)
imagesc(allRunReordered)
caxis([0 160])
colorbar
colormap(viridis)
title('offset')


%% stats on angle of approach

% sum change in angle

speedVec=[];
seshVec=[];
onsetSumAngle=[];
offsetSumAngle=[];
stateVec=[];

for isesh = 1:5
    for ispeed = 1:6
        for itrial=1:size(s(isesh).session.cond(ispeed).catData,3)
        speedVec=cat(1,speedVec,ispeed);
        seshVec=cat(1,seshVec,isesh);
        stateVec=cat(1,stateVec,1);
        onsetSumAngle=cat(1,onsetSumAngle,s(isesh).session.cond(ispeed).trial(itrial).sumDiffAngle);
        offsetSumAngle=cat(1,offsetSumAngle,s(isesh).session.cond(ispeed).trial(itrial).sumDiffAngleOffset);
    
        end
    end
end



for isesh = 1:5
    for ispeed = 7:12
        for itrial=1:size(s(isesh).session.cond(ispeed).catData,3)
        speedVec=cat(1,speedVec,ispeed-6);
        seshVec=cat(1,seshVec,isesh);
        stateVec=cat(1,stateVec,2);
        onsetSumAngle=cat(1,onsetSumAngle,s(isesh).session.cond(ispeed).trial(itrial).sumDiffAngle);
        offsetSumAngle=cat(1,offsetSumAngle,s(isesh).session.cond(ispeed).trial(itrial).sumDiffAngleOffset);
   
        end
    end
end


valsVec = offsetSumAngle;

tbl = table(valsVec,speedVec,stateVec,seshVec,...
    'VariableNames',{'vals','speed','state','sesh'});

    f = 'vals ~ state + (1|speed) + (1|sesh)';
      
    lme = fitlme(tbl, f, 'DummyVarCoding', 'reference')


    [mean(valsVec(1:30)), sem(valsVec(1:30)); mean(valsVec(31:60)), sem(valsVec(31:60))]