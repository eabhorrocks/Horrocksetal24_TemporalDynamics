%% response shape similarity using Procrustes analysis

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
    fname = [sessionTags{isession,1},'_', sessionTags{isession,2},'_FAsmooth_Mar23.mat'];

    load(fullfile(dataDir,fname))

    s(isession).session = session;
end

%%

time_idx = 21:50; % first 300ms of stim
doNorm = false;

for isession = 1:5
    dims2use = 1:s(isession).session.s.qOpt;

    for ispeed1 = 1:6

        % behavioural state comparison of responses
        Y = s(isession).session.cond(ispeed1).meanTrajectory(time_idx, dims2use);
        X = s(isession).session.cond(ispeed1+6).meanTrajectory(time_idx, dims2use);

        if doNorm
            X = zscore(X,1,2);
            Y = zscore(Y,1,2);
        end

        sesh(isession).state(ispeed1).X = X;
        sesh(isession).state(ispeed1).Y = Y;
        [sesh(isession).state(ispeed1).d,sesh(isession).state(ispeed1).Z,sesh(isession).state(ispeed1).transform] =...
            procrustes(X,Y,'Scaling',true,'reflection','best');
        sesh(isession).state(ispeed1).normTrans = norm(sesh(isession).state(ispeed1).transform.c(1,:));


        for ispeed2 = 1:6
            % speed comparison for stationary trials
            X = s(isession).session.cond(ispeed1).meanTrajectory(time_idx, dims2use);
            Y = s(isession).session.cond(ispeed2).meanTrajectory(time_idx, dims2use);

            if doNorm
                X = zscore(X,1,2);
                Y = zscore(Y,1,2);
            end

            sesh(isession).speed_stat(ispeed1,ispeed2).X = X;
            sesh(isession).speed_stat(ispeed1,ispeed2).Y = Y;
            [sesh(isession).speed_stat(ispeed1,ispeed2).d,sesh(isession).speed_stat(ispeed1,ispeed2).Z,sesh(isession).speed_stat(ispeed1,ispeed2).transform] =...
                procrustes(X,Y,'Scaling',true,'reflection','best');
            sesh(isession).speed_stat(ispeed1,ispeed2).normTrans = norm(sesh(isession).speed_stat(ispeed1,ispeed2).transform.c(1,:));

            % speed comparison for run trials
            X = s(isession).session.cond(ispeed1+6).meanTrajectory(time_idx, dims2use);
            Y = s(isession).session.cond(ispeed2+6).meanTrajectory(time_idx, dims2use);

            if doNorm
                X = zscore(X,1,2);
                Y = zscore(Y,1,2);
            end

            sesh(isession).speed_run(ispeed1,ispeed2).X = X;
            sesh(isession).speed_run(ispeed1,ispeed2).Y = Y;
            [sesh(isession).speed_run(ispeed1,ispeed2).d,sesh(isession).speed_run(ispeed1,ispeed2).Z,sesh(isession).speed_run(ispeed1,ispeed2).transform] =...
                procrustes(X,Y,'Scaling',true,'reflection','best');
            sesh(isession).speed_run(ispeed1,ispeed2).normTrans = norm(sesh(isession).speed_run(ispeed1,ispeed2).transform.c(1,:));

            % cross-state speed comparison
            Y = s(isession).session.cond(ispeed1).meanTrajectory(time_idx, dims2use);
            X = s(isession).session.cond(ispeed2+6).meanTrajectory(time_idx, dims2use);

            if doNorm
                X = zscore(X,1,2);
                Y = zscore(Y,1,2);
            end

            sesh(isession).speed_cross(ispeed1,ispeed2).X = X;
            sesh(isession).speed_cross(ispeed1,ispeed2).Y = Y;
            [sesh(isession).speed_cross(ispeed1,ispeed2).d,sesh(isession).speed_cross(ispeed1,ispeed2).Z,sesh(isession).speed_cross(ispeed1,ispeed2).transform] =...
                procrustes(X,Y,'Scaling',true,'reflection','best');
            sesh(isession).speed_cross(ispeed1,ispeed2).normTrans = norm(sesh(isession).speed_cross(ispeed1,ispeed2).transform.c(1,:));

        end
    end

end


%% get proscrustes distances
for isession = 1:5
    sesh(isession).speed_stat_d_array = setUpperTri2NaN(reshape(cat(1,sesh(isession).speed_stat.d),6,6));
    sesh(isession).speed_run_d_array = setUpperTri2NaN(reshape(cat(1,sesh(isession).speed_run.d),6,6));
    sesh(isession).state_d_array = [sesh(isession).state.d];
    sesh(isession).speed_cross_d_array = reshape(cat(1,sesh(isession).speed_cross.d),6,6);
    sesh(isession).speed_cross_d_array(1:size(sesh(isession).speed_cross_d_array,1)+1:end) = nan;

    sesh(isession).mean_speed_d_stat = mean(sesh(isession).speed_stat_d_array(:),'omitnan');
    sesh(isession).mean_speed_d_run = mean(sesh(isession).speed_run_d_array(:),'omitnan');
    sesh(isession).mean_state_d = mean(sesh(isession).state_d_array,'omitnan');
    sesh(isession).mean_crossStateSpeed_d = mean(sesh(isession).speed_cross_d_array(:),'omitnan');

end



% figure of proscrutes distances

allStatSpeed = [sesh.speed_stat_d_array];
allRunSpeed = [sesh.speed_run_d_array];
allState = [sesh.state_d_array];
allStateSpeed = [sesh.speed_cross_d_array];

%scatJit(vec, jitFactor, col, circleSize, color

figure, hold on
scatJit(allStatSpeed(:), 0.3, 1, 10, 'k')
scatJit(allRunSpeed(:), 0.3, 2, 10, 'r')
scatJit(allState(:), 0.3, 3, 10, 'm')
scatJit(allStateSpeed(:), 0.3, 4, 10, 'b')

ax=gca; ax.XTick=0:5; ax.XTickLabel = {'','Speed (stat)', 'Speed (run)', 'State', 'State and Speed',''};
ylabel('Procrustes Distance')
ax.YLim(1) = 0;
defaultAxesProperties(ax, true)
xlim([0 5])



%% what are the properties of proscrutes transformations?

% translation, scaling, rotation/reflection

%% scaling of responses
for isession = 1:5
    sesh(isession).speedScalingArray_stat = nan(6);
    sesh(isession).speedScalingArray_run = nan(6);
    sesh(isession).stateSpeedScalingArray_stat = nan(6);
    sesh(isession).stateScalingArray = nan(1,6);



    for ispeed1 = 1:6
        sesh(isession).stateScalingArray(ispeed1) = sesh(isession).state(ispeed1).transform.b;

        for ispeed2 = 1:6
            if ispeed1~=ispeed2
                sesh(isession).speedScalingArray_stat(ispeed1,ispeed2) = ...
                    sesh(isession).speed_stat(ispeed1,ispeed2).transform.b;

                sesh(isession).speedScalingArray_run(ispeed1,ispeed2) = ...
                    sesh(isession).speed_run(ispeed1,ispeed2).transform.b;

                sesh(isession).stateSpeedScalingArray(ispeed1,ispeed2) = ...
                    sesh(isession).speed_cross(ispeed1,ispeed2).transform.b;


            end
        end
    end

end

for isession = 1:5
    sesh(isession).speed_stat_b_array = sesh(isession).speedScalingArray_stat;
    sesh(isession).speed_run_b_array = sesh(isession).speedScalingArray_run;
    sesh(isession).state_b_array = sesh(isession).stateScalingArray;
    sesh(isession).speed_cross_b_array = sesh(isession).stateSpeedScalingArray;
    sesh(isession).speed_cross_b_array(1:size(sesh(isession).speed_cross_b_array,1)+1:end) = nan;

    sesh(isession).mean_speed_d_stat = mean(sesh(isession).speed_stat_b_array(:),'omitnan');
    sesh(isession).mean_speed_d_run = mean(sesh(isession).speed_run_b_array(:),'omitnan');
    sesh(isession).mean_state_d = mean(sesh(isession).state_b_array,'omitnan');
    sesh(isession).mean_crossStateSpeed_d = mean(sesh(isession).speed_cross_b_array(:),'omitnan');

end

allStatSpeed = [sesh.speed_stat_b_array];
allRunSpeed = [sesh.speed_run_b_array];
allState = [sesh.state_b_array];
allStateSpeed = [sesh.speed_cross_b_array];


figure, hold on
plot([0.5 4.5], log2([1 1]), 'k')
scatJit(log2(allStatSpeed(:)), 0.2, 1, 10, 'k')
scatJit(log2(allRunSpeed(:)), 0.2, 2, 10, 'r')
scatJit(log2(allState(:)), 0.2, 3, 10, 'm')
scatJit(log2(allStateSpeed(:)), 0.2, 4, 10, 'b')

ax=gca; ax.XTick=0:5; ax.XTickLabel = {'','Speed (stat)', 'Speed (run)', 'State', 'State and Speed',''};
ylabel('Scaling Factor')
ax.YTick = sort([[-1*log2([ 1.25 1.5, 2, 3, 4])], [log2([1 1.25 1.5, 2, 3, 4])]]);
ax.YTickLabel = sort([[1 1.25 1.5, 2, 3, 4],-1*[ 1.25 1.5, 2, 3, 4]]);
ax.YLim = [-1*log2(4), log2(4)];
% ax.YLim(1) = 0;
defaultAxesProperties(ax, true)
xlim([0 5])


%% translation of responses

for isession = 1:5
    sesh(isession).speed_stat_nT_array = setUpperTri2NaN(reshape(cat(1,sesh(isession).speed_stat.normTrans),6,6));
    sesh(isession).speed_run_nT_array = setUpperTri2NaN(reshape(cat(1,sesh(isession).speed_run.normTrans),6,6));
    sesh(isession).state_nT_array = [sesh(isession).state.normTrans];
    sesh(isession).speed_cross_nT_array = reshape(cat(1,sesh(isession).speed_cross.normTrans),6,6);
    sesh(isession).speed_cross_nT_array(1:size(sesh(isession).speed_cross_d_array,1)+1:end) = nan;

    sesh(isession).mean_speed_nT_stat = mean(sesh(isession).speed_stat_d_array(:),'omitnan');
    sesh(isession).mean_speed_nT_run = mean(sesh(isession).speed_run_d_array(:),'omitnan');
    sesh(isession).mean_state_nT = mean(sesh(isession).state_d_array,'omitnan');
    sesh(isession).mean_crossStateSpeed_nT = mean(sesh(isession).speed_cross_d_array(:),'omitnan');

end

allStatSpeed = [sesh.speed_stat_nT_array];
allRunSpeed = [sesh.speed_run_nT_array];
allState = [sesh.state_nT_array];
allStateSpeed = [sesh.speed_cross_nT_array];

%scatJit(vec, jitFactor, col, circleSize, color

figure, hold on
scatJit(allStatSpeed(:), 0.2, 1, 10, 'k')
scatJit(allRunSpeed(:), 0.2, 2, 10, 'r')
scatJit(allState(:), 0.2, 3, 10, 'm')
scatJit(allStateSpeed(:), 0.2, 4, 10, 'b')

ax=gca; ax.XTick=0:5; ax.XTickLabel = {'','Speed (stat)', 'Speed (run)', 'State', 'State and Speed',''};
ylabel('Norm of translation vector')
% ax.YLim(1) = 0;
defaultAxesProperties(ax, true)
xlim([0 5])


%% rotation of responses

for isession = 1:5
    sesh(isession).speedRotArray_stat = nan(6);
    sesh(isession).speedRotArray_run = nan(6);
    sesh(isession).stateSpeedRotArray = nan(6);
    sesh(isession).stateRotArray = nan(1,6);



    for ispeed1 = 1:6
        sesh(isession).stateRotArray(ispeed1) = norm(sesh(isession).state(ispeed1).transform.T-eye(size(sesh(isession).state(ispeed1).transform.T)),"fro");

        for ispeed2 = 1:6
            if ispeed1~=ispeed2
                sesh(isession).speedRotArray_stat(ispeed1,ispeed2) = ...
                 norm(sesh(isession).speed_stat(ispeed1,ispeed2).transform.T-eye(size(sesh(isession).speed_stat(ispeed1,ispeed2).transform.T)),"fro");

                sesh(isession).speedRotArray_run(ispeed1,ispeed2) = ...
                 norm(sesh(isession).speed_run(ispeed1,ispeed2).transform.T-eye(size(sesh(isession).speed_run(ispeed1,ispeed2).transform.T)),"fro");

                sesh(isession).stateSpeedRotArray(ispeed1,ispeed2) = ...
                 norm(sesh(isession).speed_cross(ispeed1,ispeed2).transform.T-eye(size(sesh(isession).speed_cross(ispeed1,ispeed2).transform.T)),"fro");


            end
        end
    end

end

for isession = 1:5
    sesh(isession).speed_stat_rot_array = sesh(isession).speedRotArray_stat;
    sesh(isession).speed_run_rot_array = sesh(isession).speedRotArray_run;
    sesh(isession).state_rot_array = sesh(isession).stateRotArray;
    sesh(isession).speed_cross_rot_array = sesh(isession).stateSpeedRotArray;
    sesh(isession).speed_cross_rot_array(1:size(sesh(isession).speed_cross_rot_array,1)+1:end) = nan;
end

allStatSpeed = [sesh.speed_stat_rot_array];
allRunSpeed = [sesh.speed_run_rot_array];
allState = [sesh.state_rot_array];
allStateSpeed = [sesh.speed_cross_rot_array];

%scatJit(vec, jitFactor, col, circleSize, color

figure, hold on
scatJit(allStatSpeed(:), 0.2, 1, 10, 'k')
scatJit(allRunSpeed(:), 0.2, 2, 10, 'r')
scatJit(allState(:), 0.2, 3, 10, 'm')
scatJit(allStateSpeed(:), 0.2, 4, 10, 'b')

ax=gca; ax.XTick=0:5; ax.XTickLabel = {'','Speed (stat)', 'Speed (run)', 'State', 'State and Speed',''};
ylabel('Norm of translation vector')
% ax.YLim(1) = 0;
defaultAxesProperties(ax, true)
xlim([0 5])


%% Examples of analysis

speedcols = inferno(6);
time_idx = 21:50;
isession = 4;
ispeed1 =3; ispeed2 = 4;
dims2use = 1:s(isession).session.s.qOpt;

Xstat = s(isession).session.cond(ispeed1).meanTrajectory(time_idx,dims2use);
Ystat = s(isession).session.cond(ispeed2).meanTrajectory(time_idx,dims2use);

Xrun = s(isession).session.cond(ispeed1).meanTrajectory(time_idx,dims2use);
Yrun = s(isession).session.cond(ispeed2+6).meanTrajectory(time_idx,dims2use);
%Yrun = s(isession).session.cond(ispeed1+6).meanTrajectory(time_idx,dims2use);

[d_run_speed,Z_run_speed,transform_run_speed] = procrustes(Xrun,Yrun,'Scaling',true,'reflection','best');

% figure, hold on
% plot3(Xrun(1,1),Xrun(1,2),Xrun(1,3),'o','Color', speedcols(ispeed1,:))
% plot3(Xrun(:,1),Xrun(:,2),Xrun(:,3),'Color', speedcols(ispeed1,:))
% 
% plot3(Yrun(1,1),Yrun(1,2),Yrun(1,3),'o','Color', speedcols(ispeed2,:))
% plot3(Yrun(:,1),Yrun(:,2),Yrun(:,3),'Color', speedcols(ispeed2,:))
% 
% xlabel('Factor 1'), ylabel('Factor 2'), zlabel('Factor 3')
% grid on
% view(-120, 30)

%
figure,

% rotate
Yrun_rot = Yrun*transform_run_speed.T;

ax(1) = subplot(131); hold on

plot3(Xrun(1,1),Xrun(1,2),Xrun(1,3),'^','Color', speedcols(ispeed1,:))
plot3(Xrun(end,1),Xrun(end,2),Xrun(end,3),'o','Color', speedcols(ispeed1,:))
plot3(Xrun(:,1),Xrun(:,2),Xrun(:,3),'Color', speedcols(ispeed1,:))

plot3(Yrun(1,1),Yrun(1,2),Yrun(1,3),'^','Color', speedcols(ispeed2,:))
plot3(Yrun(end,1),Yrun(end,2),Yrun(end,3),'o','Color', speedcols(ispeed2,:))
plot3(Yrun(:,1),Yrun(:,2),Yrun(:,3),'Color', speedcols(ispeed2,:))

plot3(Yrun_rot(1,1),Yrun_rot(1,2),Yrun_rot(1,3),'^','Color', speedcols(ispeed2,:))
plot3(Yrun_rot(end,1),Yrun_rot(end,2),Yrun_rot(end,3),'o','Color', speedcols(ispeed2,:))
plot3(Yrun_rot(:,1),Yrun_rot(:,2),Yrun_rot(:,3),'Color', 'm', 'LineStyle', '-')
xlabel('Factor 1'), ylabel('Factor 2'), zlabel('Factor 3')
grid on
title('Rotate')


% scale
Yrun_scale = transform_run_speed.b*Yrun*transform_run_speed.T;

ax(2)= subplot(132); hold on

plot3(Xrun(1,1),Xrun(1,2),Xrun(1,3),'^','Color', speedcols(ispeed1,:))
plot3(Xrun(end,1),Xrun(end,2),Xrun(end,3),'o','Color', speedcols(ispeed1,:))
plot3(Xrun(:,1),Xrun(:,2),Xrun(:,3),'Color', speedcols(ispeed1,:))

plot3(Yrun_rot(1,1),Yrun_rot(1,2),Yrun_rot(1,3),'^','Color', speedcols(ispeed2,:))
plot3(Yrun_rot(end,1),Yrun_rot(end,2),Yrun_rot(end,3),'o','Color', speedcols(ispeed2,:))
plot3(Yrun_rot(:,1),Yrun_rot(:,2),Yrun_rot(:,3),'Color', speedcols(ispeed2,:))

plot3(Yrun_scale(1,1),Yrun_scale(1,2),Yrun_scale(1,3),'^','Color', speedcols(ispeed2,:))
plot3(Yrun_scale(end,1),Yrun_scale(end,2),Yrun_scale(end,3),'o','Color', speedcols(ispeed2,:))
plot3(Yrun_scale(:,1),Yrun_scale(:,2),Yrun_scale(:,3),'Color', 'm', 'LineStyle', '-')
xlabel('Factor 1'), ylabel('Factor 2'), zlabel('Factor 3')
grid on
view(-145, 20)
title('Scale')

% translate
ax(3) = subplot(133); hold on
Yrun_trans = Yrun_scale + transform_run_speed.c;

plot3(Xrun(1,1),Xrun(1,2),Xrun(1,3),'^','Color', speedcols(ispeed1,:))
plot3(Xrun(end,1),Xrun(end,2),Xrun(end,3),'o','Color', speedcols(ispeed1,:))
plot3(Xrun(:,1),Xrun(:,2),Xrun(:,3),'Color', speedcols(ispeed1,:))

plot3(Yrun_scale(1,1),Yrun_scale(1,2),Yrun_scale(1,3),'^','Color', speedcols(ispeed2,:))
plot3(Yrun_scale(end,1),Yrun_scale(end,2),Yrun_scale(end,3),'o','Color', speedcols(ispeed2,:))
plot3(Yrun_scale(:,1),Yrun_scale(:,2),Yrun_scale(:,3),'Color', speedcols(ispeed2,:))

plot3(Yrun_trans(1,1),Yrun_trans(1,2),Yrun_trans(1,3),'^','Color', speedcols(ispeed2,:))
plot3(Yrun_trans(end,1),Yrun_trans(end,2),Yrun_trans(end,3),'o','Color', speedcols(ispeed2,:))
plot3(Yrun_trans(:,1),Yrun_trans(:,2),Yrun_trans(:,3),'Color', 'm', 'LineStyle', '-')
xlabel('Factor 1'), ylabel('Factor 2'), zlabel('Factor 3')
grid on
title('Translate')


for iplot = 1:3
    subplot(1,3,iplot)
    view(-150, 30)

    allYLim = get(ax, {'YLim'});
    allYLim = cat(2, allYLim{:});
    set(ax, 'YLim', [min(allYLim), max(allYLim)]);

    allXLim = get(ax, {'XLim'});
    allXLim = cat(2, allXLim{:});
    set(ax, 'XLim', [min(allXLim), max(allXLim)]);

    allZLim = get(ax, {'ZLim'});
    allZLim = cat(2, allZLim{:});
    set(ax, 'ZLim', [min(allZLim), max(allZLim)]);
end

   
