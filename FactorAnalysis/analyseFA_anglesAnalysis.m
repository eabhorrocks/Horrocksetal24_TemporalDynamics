%% FA angle analysis

%% load data
sessionTags = {'M22027', '20220517';...
    'M22029', '20220607';...
    'M22031', '20220524';
    'M22032', '20220621';
    'M22033', '20220706'};

dataDir = 'C:\Users\edward.horrocks\Documents\Code\V1Dynamics\Data\basic_111022';
dataDir = 'E:\V1Data\Data\v1_fromC24';
s = struct;

for isession = 1:5
    tic
    fname = [sessionTags{isession,1},'_', sessionTags{isession,2},'_FA_normal.mat'];

    load(fullfile(dataDir,fname))

    s(isession).session = session;
end
%%



speedcols=inferno(6);
%% angle of trajectory vs sustained period


for isession = 1:5
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


        traj = s(isession).session.cond(icond).meanTrajectory(:, 1:s(isession).session.s.qOpt);


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


        s(isession).session.cond(icond).angle = acosd(angle);%mod(acosd(angle).*sign(angle),360);
        s(isession).session.cond(icond).angleOffset = acosd(angleOffset);%mod(acosd(angle).*sign(angle),360);

        s(isession).session.cond(icond).distFromStim = distFromSS;
        s(isession).session.cond(icond).distFromStim = distFromSSOffset;

        s(isession).session.cond(icond).sumDiffAngle = sum(abs(diff(s(isession).session.cond(icond).angle(1:111))));
        s(isession).session.cond(icond).sumDiffAngleOffset = sum(abs(diff(s(isession).session.cond(icond).angleOffset(111:end))));

        s(isession).session.cond(icond).v1full = v1full;
        % s(isession).session.cond(icond).angleSub = angleSub;



    end
end



%% sum change in angle

speedVec=[];
seshVec=[];
onsetSumAngle=[];
offsetSumAngle=[];
stateVec=[];

for isesh = 1:5
    for ispeed = 1:6
        speedVec=cat(1,speedVec,ispeed);
        seshVec=cat(1,seshVec,isesh);
        stateVec=cat(1,stateVec,1);
        onsetSumAngle=cat(1,onsetSumAngle,s(isesh).session.cond(ispeed).sumDiffAngle);
        offsetSumAngle=cat(1,offsetSumAngle,s(isesh).session.cond(ispeed).sumDiffAngleOffset);
    end
end



for isesh = 1:5
    for ispeed = 7:12
        speedVec=cat(1,speedVec,ispeed-6);
        seshVec=cat(1,seshVec,isesh);
        stateVec=cat(1,stateVec,2);
        onsetSumAngle=cat(1,onsetSumAngle,s(isesh).session.cond(ispeed).sumDiffAngle);
        offsetSumAngle=cat(1,offsetSumAngle,s(isesh).session.cond(ispeed).sumDiffAngleOffset);
    end
end


valsVec = offsetSumAngle;

tbl = table(valsVec,speedVec,stateVec,seshVec,...
    'VariableNames',{'vals','speed','state','sesh'});

    f = 'vals ~ state + (1|speed) + (1|sesh)';
      
    lme = fitlme(tbl, f, 'DummyVarCoding', 'reference')


    [mean(valsVec(1:30)), sem(valsVec(1:30)); mean(valsVec(31:60)), sem(valsVec(31:60))]


    
%% plot overall change in angle
figure
subplot(121), hold on
for ispeed = 1:6
    for isession = 1:5

        plot(s(isession).session.cond(ispeed).sumDiffAngle, s(isession).session.cond(ispeed+6).sumDiffAngle,'o', 'MarkerEdgeColor','none','MarkerFaceColor',speedcols(ispeed,:));
    end
end

subplot(122), hold on
for ispeed = 1:6
    for isession = 1:5

        plot(s(isession).session.cond(ispeed).sumDiffAngleOffset, s(isession).session.cond(ispeed+6).sumDiffAngleOffset,'o', 'MarkerEdgeColor','none','MarkerFaceColor',speedcols(ispeed,:));
    end
end



subplot(121)
plot([0 600],[0 600],'k:')
xlim([0 600])
ylim([0 600])
axis equal
defaultAxesProperties(gca, true)

subplot(122)
plot([0 600],[0 600],'k:')
xlim([0 600])
ylim([0 600])
axis equal
defaultAxesProperties(gca, true)

%% plot average angle over time

allStat=[];
allRun=[];

for isession = 1:5
    for ispeed = 1:6
        allStat=cat(1,allStat, s(isession).session.cond(ispeed).angle);
        allRun=cat(1,allRun, s(isession).session.cond(ispeed+6).angle);
    end
end

figure
shadedErrorBar(-195:10:1795,mean(allStat,1), sem(allStat,1), 'lineProps','k')
shadedErrorBar(-195:10:1795,mean(allRun,1), sem(allRun,1), 'lineProps','r')
defaultAxesProperties(gca, true)
title('Onset')

allStat=[];
allRun=[];

for isession = 1:5
    for ispeed = 1:6
        allStat=cat(1,allStat, s(isession).session.cond(ispeed).angleOffset);
        allRun=cat(1,allRun, s(isession).session.cond(ispeed+6).angleOffset);
    end
end

figure
shadedErrorBar(-195:10:1795,mean(allStat,1), sem(allStat,1), 'lineProps','k')
shadedErrorBar(-195:10:1795,mean(allRun,1), sem(allRun,1), 'lineProps','r')
ylim([0 140])
defaultAxesProperties(gca, true)
title('Offset')


%% plot by speed

figure
% plot by speed
for ispeed = 1:6
    allStat=[];
    allRun=[];
for isession=1:5
    allStat = cat(1,allStat,s(isession).session.cond(ispeed).angle);
    allRun = cat(1,allRun,s(isession).session.cond(ispeed+6).angle);
end
subplot(221), hold on
plot(-195:10:1795, mean(allStat,1),'Color', speedcols(ispeed,:))
subplot(223), hold on
plot(-195:10:1795, mean(allRun,1),'Color', speedcols(ispeed,:))
end



% plot by speed
for ispeed = 1:6
    allStat=[];
    allRun=[];
for isession=1:5
    allStat = cat(1,allStat,s(isession).session.cond(ispeed).angleOffset);
    allRun = cat(1,allRun,s(isession).session.cond(ispeed+6).angleOffset);
end
subplot(222), hold on
plot(-195:10:1795, mean(allStat,1),'Color', speedcols(ispeed,:))
subplot(224), hold on
plot(-195:10:1795, mean(allRun,1),'Color', speedcols(ispeed,:))
end

subplot(221)
xlim([-200 1800])
ylim([0 160])
defaultAxesProperties(gca, true)

subplot(223)
xlim([-200 1800])
ylim([0 160])
defaultAxesProperties(gca, true)

subplot(222)
xlim([-200 1800])
ylim([0 160])
defaultAxesProperties(gca, true)

subplot(224)
xlim([-200 1800])
ylim([0 160])
defaultAxesProperties(gca, true)

%%

% plot example
ispeed = 3;
isession = 2;

figure, subplot(211),hold on
plot(-195:10:1795,s(isession).session.cond(ispeed).angle,'k')
plot(-195:10:1795,s(isession).session.cond(ispeed+6).angle,'r')
xlabel('Time'), ylabel('Angle from direct path between baseline and stimulus steady-state')
defaultAxesProperties(gca, false)

subplot(212),hold on
plot(-195:10:1795,s(isession).session.cond(ispeed).distFromStim,'k')
plot(-195:10:1795,s(isession).session.cond(ispeed+6).distFromStim,'r')
xlabel('Time'), ylabel('Distance from stimulus steady-state')
defaultAxesProperties(gca, false)


speedcols = inferno(6);
statVals=[];
runVals=[];

figure, hold on
plot([300 900], [300, 900], 'k:')
xlabel('stat'), ylabel('run'), title('Absolute sum of dtheta/dt')
for isession = 1:5
    for ispeed = 1:6
        statVals=cat(1,statVals, s(isession).session.cond(ispeed).sumDiffAngle);
        runVals=cat(1,runVals, s(isession).session.cond(ispeed+6).sumDiffAngle);

        plot(statVals(end), runVals(end), 'o', 'MarkerEdgeColor','none',...
            'MarkerFaceColor', speedcols(ispeed,:))
    end


end
defaultAxesProperties(gca, false)



%% plot example 2d trajectories
% session =1, dims = [2,3], speed = 2
% session =2, dims = [2,4], speed = 2
% session = 2, dims = [1 2], speed = 5
% 3, [1 2], 5
% 2     3     5     5
%     2     1     3     5
% 
%     4     1     3     3
% 4 1 3 4
% 4 2 3 4
for isession = 4
for ispeed =4
    figure

dims = [2 3];

[isession, dims, ispeed]
time = 1:200;%find(ismember(binVector,0:));
time2plot=21:70;
cols = magma(numel((time2plot)));

% stat
refTraj = s(isession).session.cond(ispeed).v1full;
u = refTraj(time, dims(1));
v = refTraj(time, dims(2));
ss_mean = mean([u(finalSSIdx),v(finalSSIdx)],1);
iniAngle = 90-rad2deg(cart2pol(u(1), v(1)));
ss_mean = mean([u(finalSSIdx),v(finalSSIdx)],1);

% run
refTraj2 = s(isession).session.cond(ispeed+6).v1full;
u2 = refTraj2(time, dims(1));
v2 = refTraj2(time, dims(2));
ss_mean2 =  mean([u2(finalSSIdx),v2(finalSSIdx)],1);
iniAngle2 = 90-rad2deg(cart2pol(u2(1), v2(1)));


% view(-iniAngle,90)
% axis equal
subplot(121)
plot(u(time2plot),v(time2plot),'k','Marker','.'),
hold on,
plot(u(time2plot(1)),v(time2plot(1)),'kv'),
plot(u(time2plot(end)),v(time2plot(end)),'ko')
plot(ss_mean(1), ss_mean(2),'k*')
% view(-iniAngle,90)
xlabel(['Dim: ', num2str(dims(1))]), ylabel(['Dim: ', num2str(dims(2))])
title([isession, ispeed])


subplot(122)
plot(u2(time2plot),v2(time2plot),'r','Marker','.'),
hold on,
plot(u2(time2plot(1)),v2(time2plot(1)),'rv'),
plot(u2(time2plot(end)), v2(time2plot(end)),'ro')
plot(ss_mean2(1), ss_mean2(2),'r*')

xlabel(['Dim: ', num2str(dims(1))]), ylabel(['Dim: ', num2str(dims(2))])
% view(-iniAngle2,90)

end
end


%% plot example 2d with single trials

for isession = 4
for ispeed =3
    figure

dims = [1 3];

[isession, dims, ispeed]
time = 1:200;%find(ismember(binVector,0:));
time2plot=21:70;
cols = magma(numel((time2plot)));

% stat
refTraj = s(isession).session.cond(ispeed).meanTrajectory;
u = refTraj(time, dims(1));
v = refTraj(time, dims(2));

subplot(121), hold on
% single trials
for itrial = 1:31

    plot(s(isession).session.cond(ispeed).catData(time2plot(1),dims(1),itrial),s(isession).session.cond(ispeed).catData(time2plot(1),dims(2),itrial), 'Color', [.7 .7 .7], 'Marker', 'v')
    plot(s(isession).session.cond(ispeed).catData(time2plot(end),dims(1),itrial),s(isession).session.cond(ispeed).catData(time2plot(end),dims(2),itrial), 'Color', [.7 .7 .7], 'Marker', 'o')
    plot(s(isession).session.cond(ispeed).catData(time2plot,dims(1),itrial),s(isession).session.cond(ispeed).catData(time2plot,dims(2),itrial), 'Color', [.7 .7 .7], 'LineWidth',0.1)
end
% average
plot(u(time2plot),v(time2plot),'k','Marker','.'),
hold on,
plot(u(time2plot(1)),v(time2plot(1)),'kv'),
plot(u(time2plot(end)),v(time2plot(end)),'ko')
% view(-iniAngle,90)
xlabel(['Dim: ', num2str(dims(1))]), ylabel(['Dim: ', num2str(dims(2))])
title([isession, ispeed])

% plot single trials



% run
refTraj2 = s(isession).session.cond(ispeed+6).meanTrajectory;
u2 = refTraj2(time, dims(1));
v2 = refTraj2(time, dims(2));
ss_mean2 =  mean([u2(finalSSIdx),v2(finalSSIdx)],1);
iniAngle2 = 90-rad2deg(cart2pol(u2(1), v2(1)));


% view(-iniAngle,90)
% axis equal
subplot(122), hold on
for itrial = 1:31

    plot(s(isession).session.cond(ispeed+6).catData(time2plot(1),dims(1),itrial),s(isession).session.cond(ispeed+6).catData(time2plot(1),dims(2),itrial), 'Color', [.7 .7 .7], 'Marker', 'v')
    plot(s(isession).session.cond(ispeed+6).catData(time2plot(end),dims(1),itrial),s(isession).session.cond(ispeed+6).catData(time2plot(end),dims(2),itrial), 'Color', [.7 .7 .7], 'Marker', 'o')
    plot(s(isession).session.cond(ispeed+6).catData(time2plot,dims(1),itrial),s(isession).session.cond(ispeed+6).catData(time2plot,dims(2),itrial), 'Color', [.7 .7 .7], 'LineWidth', 0.1)
end

plot(u2(time2plot),v2(time2plot),'r','Marker','.'),
hold on,
plot(u2(time2plot(1)),v2(time2plot(1)),'rv'),
plot(u2(time2plot(end)), v2(time2plot(end)),'ro')

xlabel(['Dim: ', num2str(dims(1))]), ylabel(['Dim: ', num2str(dims(2))])
% view(-iniAngle2,90)

end
end

%% 

seshVec = [1 2 2 3 2 2 4 4 4];
dimCell = {[2 3], [2, 4], [1, 2], [1, 2], [3, 5], [1, 3], [1, 3], [1, 3], [2, 3]};
speedVec = [2 2 5 5 5 5 3 4 4];
time2plot = 21:70;

for iex = 1:numel(seshVec)

isession = seshVec(iex);
dims = dimCell{iex};
ispeed = speedVec(iex);

figure, hold on
plot(s(isession).session.cond(ispeed).meanTrajectory(time2plot,dims(1)),s(isession).session.cond(ispeed).meanTrajectory(time2plot,dims(2)),'k','Marker','.');
plot(s(isession).session.cond(ispeed).meanTrajectory(time2plot(1),dims(1)),s(isession).session.cond(ispeed).meanTrajectory(time2plot(1),dims(2)),'kv');
plot(s(isession).session.cond(ispeed).meanTrajectory(time2plot(end),dims(1)),s(isession).session.cond(ispeed).meanTrajectory(time2plot(end),dims(2)),'ko');

plot(s(isession).session.cond(ispeed+6).meanTrajectory(time2plot,dims(1)),s(isession).session.cond(ispeed+6).meanTrajectory(time2plot,dims(2)),'r','Marker','.');
plot(s(isession).session.cond(ispeed+6).meanTrajectory(time2plot(1),dims(1)),s(isession).session.cond(ispeed+6).meanTrajectory(time2plot(1),dims(2)),'rv');
plot(s(isession).session.cond(ispeed+6).meanTrajectory(time2plot(end),dims(1)),s(isession).session.cond(ispeed+6).meanTrajectory(time2plot(end),dims(2)),'ro');

xlabel(dims(1)), ylabel(dims(2))
title(['Mouse: ', num2str(isession), ', Speed: ', num2str(ispeed)])
defaultAxesProperties(gca, false)
end




%% compass plot

% session =1, dims = [2,3], speed = 2
% session =2, dims = [2,4], speed = 2
% session = 2, dims = [1 2], speed = 5
% 3, [1 2], 5
% 2     3     5     5
%     2     1     3     5
% 
%     4     1     3     3
% 4 1 3 4
% 4 2 3 4
for isession =2
for ispeed =5
    figure

dims = [3 5 ];

[isession, dims, ispeed]
time = 1:200;%find(ismember(binVector,0:));
time2plot=21:70;
cols = magma(numel((time2plot)));


refTraj = s(isession).session.cond(ispeed).v1full;
u = refTraj(time, dims(1));
v = refTraj(time, dims(2));
ss_mean = mean([u(finalSSIdx),v(finalSSIdx)],1);
iniAngle = 90-rad2deg(cart2pol(u(1), v(1)));

refTraj2 = s(isession).session.cond(ispeed+6).v1full;
u2 = refTraj2(time, dims(1));
v2 = refTraj2(time, dims(2));
ss_mean2 =  mean([u2(finalSSIdx),v2(finalSSIdx)],1);
iniAngle2 = 90-rad2deg(cart2pol(u2(1), v2(1)));


subplot(231)
sh = compass(u(time2plot), v(time2plot));
%  view(-iniAngle,90)
% axis equal
subplot(232)
plot(u(time2plot),v(time2plot),'k','Marker','.'),
hold on,
plot(u(time2plot(1)),v(time2plot(1)),'kv'),
plot(u(time2plot(end)),v(time2plot(end)),'ko')
plot(ss_mean(1), ss_mean(2),'k*')
 view(-iniAngle,90)
xlabel(['Dim: ', num2str(dims(1))]), ylabel(['Dim: ', num2str(dims(2))])
subplot(233)
plot(binVector(time2plot),s(isession).session.cond(ispeed).angle(time2plot),'k')
xlim(binVector([time2plot(1), time2plot(end)]))

xlabel('Time'), ylabel('Angle from direct path in full factor space')

subplot(234)
rh = compass(u2(time2plot),v2(time2plot));
%  view(-iniAngle2,90)
% axis equal
subplot(235)
plot(u2(time2plot),v2(time2plot),'r','Marker','.'),
hold on,
plot(u2(time2plot(1)),v2(time2plot(1)),'rv'),
plot(u2(time2plot(end)), v2(time2plot(end)),'ro')
plot(ss_mean2(1), ss_mean2(2),'r*')

xlabel(['Dim: ', num2str(dims(1))]), ylabel(['Dim: ', num2str(dims(2))])
 view(-iniAngle2,90)

subplot(236)
plot(binVector(time2plot),s(isession).session.cond(ispeed+6).angle(time2plot),'r')
xlim(binVector([time2plot(1), time2plot(end)]))
xlabel('Time'), ylabel('Angle from direct path in full factor space')

for itime = 1:numel(time2plot)
    sh(itime).Color = cols(itime,:);
    rh(itime).Color = cols(itime,:);
end

pause

end

end
%% animation plots

figure
for itime = 1:200
    subplot(121)
    hold off
    max_lim = 3;
    x_fake=[0 max_lim 0 -max_lim];
    y_fake=[max_lim 0 -max_lim 0];
    h_fake=compass(x_fake,y_fake);
    hold on
    compass(u(itime),v(itime))
    set(h_fake,'Visible','off');
    title(num2str(binVector(itime)))

    subplot(122)
    hold off
    max_lim=6;
    x_fake=[0 max_lim 0 -max_lim];
    y_fake=[max_lim 0 -max_lim 0];
    h_fake=compass(x_fake,y_fake);
    hold on
    compass(u2(itime),v2(itime))
    set(h_fake,'Visible','off');
    title(num2str(binVector(itime)))
    pause(0.1)
end



%% rotation of response over time



for isession = 1:5
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


        traj = s(isession).session.cond(icond).meanTrajectory(:, 1:s(isession).session.s.qOpt);


        ss_mean = mean(traj(finalSSIdx,:),1);
        ini_mean = mean(traj(initialSSIdx,:),1);

        refAxis = ss_mean-ini_mean;

        angle=[];

        for itime = 2:199
            v1 = traj(itime,:)-traj(itime-1,:);
            v2 = traj(itime+1,:)-traj(itime,:);
           angle(itime-1) = acosd(calcCosineTheta(v1,v2));
        end

        s(isession).session.cond(icond).velAngle = angle;

    end
end

%% plot vel-angle
allStat=[];
allRun=[];

for isession = 1:5
    for ispeed = 1:6
        allStat=cat(1,allStat, s(isession).session.cond(ispeed).velAngle);
        allRun=cat(1,allRun, s(isession).session.cond(ispeed+6).velAngle);
    end
end

figure
shadedErrorBar(-180:10:1790,mean(allStat,1), sem(allStat,1), 'lineProps','k')
shadedErrorBar(-180:10:1790,mean(allRun,1), sem(allRun,1), 'lineProps','r')
defaultAxesProperties(gca, true)


%% plot individual factors


for isession = 1:5
    for ispeed = 1:6
isession
ispeed
        figure
        for idim = 1:8
        subplot(2,4,idim), hold on
        plot(s(isession).session.cond(ispeed).meanTrajectory(21:70,idim),'k')
        plot(s(isession).session.cond(ispeed+6).meanTrajectory(21:70,idim),'r')
        end
        pause
    end
end

