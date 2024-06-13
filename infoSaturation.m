%% info saturation estimation based on kafashan et al

% to do
% - load basic.mat
% - partition trials by state (runFlag)
% - get binned spike counts
% - loop through each speed combination (just neighbouring?), w/ and w/o
% shuffling ?one set of trials


%% loada data


sessionTags = {'M22027', '20220517';...
    'M22029', '20220607';...
    'M22031', '20220524';
    'M22032', '20220621';
    'M22033', '20220706'};

dataDir = 'C:\Users\edward.horrocks\Documents\Code\V1Dynamics\Data\basic_111022';


for isession = 1:size(sessionTags,1)
    fname = [sessionTags{isession,1},'_', sessionTags{isession,2},'_basic.mat']; %PSTH_noEye %psth3rds

    load(fullfile(dataDir,fname))

    isession

%% pre-process data - classify trials by behavioural state

% remove any trials that didnt run for correcgt amount of time
trialDur = 1;
tolerance = 0.05; % within 50ms
tsd = trials.Speed2D;
trialdurs = [tsd.PDend]-[tsd.PDstart];
invalidDurs_idx = abs(trialdurs-1)>tolerance;
tsd(invalidDurs_idx)=[];

tsd = tsd([tsd.Contrast1]==1 & [tsd.numDots1]==573); % only using full contrast trials


stateTrialType = 'normal'; % 'normal', 'strict', 'changepoints'

switch stateTrialType
    case 'normal'

        % normal paper criteria
        stat_idx = find(cellfun(@(x) prop(x<3)>=0.75 & mean(x)<0.5, {tsd.WheelSpeed}));
        run_idx = find(cellfun(@(x) prop(x>0.5)>=0.75 & mean(x)>3, {tsd.WheelSpeed}));

    case 'strict'
        % stricter criteria
        stat_idx = find(cellfun(@(x) prop(x<0.5)>=1 & mean(x)<0.5, {tsd.WheelSpeed}));
        run_idx = find(cellfun(@(x) prop(x>0.5)>=0.9 & mean(x)>3, {tsd.WheelSpeed}));

    case 'changepoints'


        wheelOn=[]; wheelOff=[];
        for iwheel = 1:numel(wheel)
            wheel(iwheel).rawSpeedInterp = cat(1,0, wheel(iwheel).rawSpeedInterp);
            wheel(iwheel).eTimeInterp = cat(1,wheel(iwheel).eTimeInterp(1)-0.01, wheel(iwheel).eTimeInterp);
            wheel(iwheel).rawSpeedInterp(end+1) = 0;
            wheel(iwheel).eTimeInterp(end+1) = wheel(iwheel).eTimeInterp(end)+0.01;


            wheelOn(iwheel) = wheel(iwheel).eTimeInterp(1);
            wheelOff(iwheel) = wheel(iwheel).eTimeInterp(end);
        end

        wheelSpeed = cat(1,wheel.rawSpeedInterp);
        wheelZSpd = zscore(wheelSpeed);
        wheelTime = cat(1,wheel.eTimeInterp);
        data = wheelZSpd;
        zThres = 0.005; % moving standard deviations exceeded/fell below an empirical threshold of 0.005
        timestamps = wheelTime;
        inSampleRate = 100;
        smoothWin = 2; % moving standard deviation of speed (2s in Lohani)
        changeDur=5; % minimum duration of the state change in seconds (5s in Lohani)
        timeBetween=0.5; % minimum time between off and the next on in seconds

        % lohani use 3s buffer from start/end points for sustained


        [~,OnTStamp ,OffTStamp ] =changepoints(data, zThres,timestamps,inSampleRate,...
            smoothWin,changeDur,timeBetween);

        % check if intervals occur inbetween distinct recordings - if so, set off
        % to end time on wheel/on to 1st wheel time
        nEpochs = numel(OnTStamp);

        for iepoch = 1:nEpochs
            for iwheel = 1:numel(wheel)
                % if epoch goes passed where wheel recording ends (for this stim set)
                % set its end time to wheel end + create a new epoch that starts
                % with next wheel
                if OnTStamp(iepoch)<wheelOff(iwheel) && OffTStamp(iepoch)>wheelOff(iwheel)
                    tempOff = OffTStamp(iepoch);
                    % end this epoch at end of wheel
                    OffTStamp(iepoch)=wheelOff(iwheel);
                    % create new epoch at start of next wheel and with original end
                    OnTStamp(numel(OnTStamp)+1) = wheelOn(iwheel+1);
                    OffTStamp(numel(OffTStamp)+1) = tempOff;
                end
            end

        end

        locoInterval_meanSpeed=[];
        locoInterval_duration=[];
        locoInterval_timeSinceLast=[];

        for i=1:numel(OnTStamp)
            locoInterval_meanSpeed(i) = mean(wheelSpeed(find(wheelTime==OnTStamp(i)):find(wheelTime==OffTStamp(i))));
            locoInterval_duration(i) = OffTStamp(i)-OnTStamp(i);
            if i>1
                locoInterval_timeSinceLast(i) = OnTStamp(i)-OffTStamp(i-1);
            end
        end

        % remove intervals that don't meet criteria
        toDelIdx = find([locoInterval_meanSpeed]<3 | [locoInterval_duration]<5);

        OnTStamp(toDelIdx)=[];
        OffTStamp(toDelIdx)=[];

        % find trials that start and end within a locomotion epoch
        % loop through each locomotion epoch and find valid trials
        nEpochs = numel(OnTStamp);

        allStartTimes = cat(1,tsd.PDstart)-0.2;
        allEndTimes = cat(1,tsd.PDend)+0.8;

        runIdx = [];
        % must start and finish within the time range

        OnTStamp = OnTStamp+0.5; % remove smaall buffer from start and end of epochs
        OffTStamp = OffTStamp-0.5;

        for iepoch = 1:nEpochs
            runIdx_temp = find(allStartTimes>OnTStamp(iepoch) & allStartTimes<OffTStamp(iepoch)...
                & allEndTimes>OnTStamp(iepoch) & allEndTimes<OffTStamp(iepoch));

            runIdx=cat(1,runIdx,runIdx_temp);
        end
        run_idx = runIdx;
        % use stricter stationary trial criteria
        stat_idx = find(cellfun(@(x) prop(x<0.5)>=1 & mean(x)<0.5, {tsd.WheelSpeed}));
end


[tsd.runFlag] = deal([nan]);

[tsd(stat_idx).runFlag] = deal([0]);
[tsd(run_idx).runFlag] = deal([1]);
tsd = tsd(~isnan([tsd.runFlag]));



%% get binned spikes for all trials
dmtrials = tsd;

for itrial =1 :numel(dmtrials)
    dmtrials(itrial).start_time = dmtrials(itrial).PDstart;
    dmtrials(itrial).absVel = abs(dmtrials(itrial).VelX1);
end

for iunit = 1:numel(units)
    units(iunit).spiketimes = units(iunit).spike_times;
end

options.intervalStart = -0.2;
options.intervalEnd = 1.8;
options.binSpacing=0.01; % 10ms bin spacing.

[anUnits, cond] = getBinnedSpikeCounts(dmtrials, units, {'absVel', 'runFlag'}, options);

for iunit = 1:numel(units)
    units(iunit).allSpikes = anUnits(iunit).allSpikes;
end

%%
minTrial = min(cellfun(@(x) size(x,2), units(1).allSpikes(:)));

for iunit = 1:numel(units)
    units(iunit).allSpikes = cellfun(@(x) x(:,1:minTrial), units(iunit).allSpikes, 'UniformOutput', false);
end

%% get time-avaerage values for now

time_idx = 21:121;

for iunit = 1:numel(units)
    units(iunit).avS = cellfun(@(x) mean(x(time_idx),1), units(iunit).allSpikes,'UniformOutput',false);
   allSpikes = vertcat(units(iunit).allSpikes{:});
    units(iunit).dm_fr = mean(allSpikes(:)).*100;

end

  
%%
goodUnits = units([units.isi_viol]<=0.1...
    & [units.amplitude_cutoff]<=0.1 & [units.amplitude]>=50 & [units.dm_fr]>=1);

% remove low spiking units?
%%
% sr is spike rate (nTrial x nCell)
% ori is nTrial vector of orientations
% con is nTrial vector of contrasts (fixed for m25a at least)
% oris(ori1id), oris(ori2id) are the actual orientations (e.g 45, 90)
% subdims is number of cells (could be optionally less if wanted?)
% min([T Tb]) is inf in the case I tested (nTrials related?)



Nmax = 25;
clear state
for istate = 1:2
    for ishuf = 1:2
        if ishuf==2
            shuf = true;
        else
            shuf = false;
        end

        subdims = 999999;
        popsamples = 10000;

        sr = [];
        for iunit = 1:numel(goodUnits)
            allResps = [goodUnits(iunit).avS{:,istate}];
            sr(:,iunit) = allResps;
        end

        velVec = repelem(1:6,cellfun(@numel, goodUnits(1).avS(:,istate)));
        velVec = velVec(:);

        %
        %     for ivel = 1:5
        for ivel1 = 1:5
            %for ivel2 = 1:6
            v1 = ivel1;
            v2=v1+1; 
            ivel2 = 1;
            %v2 = ivel2;

            sr1 = sr(velVec==v1, :); % pop vector for all stim1 trials
            sr2 = sr(velVec==v2, :); % pop vector for all stim2 trials

            if shuf
                sr2temp = [];
                for iunit = 1:size(sr2,2)
                    sr2temp(:,iunit) = sr2(randperm(size(sr2,1)));
                end
                sr2 = sr2temp;
            end

            % get activity from T trials
            T = min([size(sr1,1),size(sr2,1)]); % or downsample?
            sr1 = sr1(1:T, :);
            sr2 = sr2(1:T, :);


            % ds is difference in speeds (log2 deg/s?) important for units of Fisher
            % Info
            ds = abs(v1 - v2);

            mu = (mean(sr1, 1) - mean(sr2, 1)) / ds;
            S = 0.5 * cov(sr1) + 0.5 * cov(sr2);
            N = length(mu);
           % Nmax = min(2*T-3, N);  % avoid infinite var estimator
            % subdimensions, if requested
            if subdims < N
                % order eigenvectors by decreasing order of eigenvalues
                [Q,Lam] = eig(S);
                [~,i] = sort(diag(Lam),'descend');
                Q = Q(:,i);
                % map into lower-dimensional subspace
                QM = Q(:,1:subdims) * Q(:,1:subdims)';
                mu = mu * QM;
                S = QM * S * QM;
                Nmax = min(subdims, Nmax);
            end


            % compute moments
            % inputs are the outputs of datamoments. popsamples is the number
            % of random orders to add cells

            fprintf('Computing information estimates (N=%d, Nmax=%d) ...\n', N, Nmax);
            [Iincr_mu, Iincr_var, Iincr_samples] = ...
                estIincrMoments(mu, S, T, ds, popsamples, Nmax);

            I_var = cumsum(Iincr_var);
            I_mu = cumsum(Iincr_mu);

            state(istate).shuf(ishuf).vel(ivel1,ivel2).I_mu = I_mu;
            state(istate).shuf(ishuf).vel(ivel1,ivel2).I_var = I_var;

            end
        %end

    end
end


s(isession).state = state;
end

%
%     end
% end

% % plot result
% 
% % plot moments
% figure('Color', 'white');
% subplot(2, 1, 1);  hold on;
% patch([1:Nmax fliplr(1:Nmax)], ...
%     [(Iincr_mu+sqrt(Iincr_var)) fliplr(Iincr_mu-sqrt(Iincr_var))],1,...
%     'FaceColor', [0.8 0.8 0.8], 'EdgeColor', 'none');
% plot(1:Nmax, Iincr_mu, 'k-', 'LineWidth', 2);
% plot([1 Nmax], [0 0], 'k--');
% ylabel('Info increase estimate');
% 
% subplot(2,1,2)
% hold on;
% I_var = cumsum(Iincr_var);
% I_mu = cumsum(Iincr_mu);
% patch([1:Nmax fliplr(1:Nmax)], ...
%     [(I_mu+sqrt(I_var)) fliplr(I_mu-sqrt(I_var))],1,...
%     'FaceColor', [0.8 0.8 0.8], 'EdgeColor', 'none');
% plot(1:Nmax, I_mu, 'k-', 'LineWidth', 2);
% xlabel('N');
% ylabel('Info estimate');
% 
% end
% end
% 



%% plot results

plotFlag = true
if plotFlag
for isession =1:5
    state = s(isession).state;

% for ivel = 1:6
%     state(1).shuf(1).vel(ivel,ivel).I_mu = nan(size(state(1).shuf(1).vel(ivel,ivel).I_mu));
%     state(1).shuf(2).vel(ivel,ivel).I_mu = nan(size(state(1).shuf(1).vel(ivel,ivel).I_mu));
%     state(2).shuf(1).vel(ivel,ivel).I_mu = nan(size(state(1).shuf(1).vel(ivel,ivel).I_mu));
%     state(2).shuf(2).vel(ivel,ivel).I_mu = nan(size(state(1).shuf(1).vel(ivel,ivel).I_mu));
% end

figure
shadedErrorBar(1:Nmax, nanmean(cat(1,state(1).shuf(1).vel.I_mu),1), nansem(cat(1,state(1).shuf(1).vel.I_mu),1),'lineProps','k')
shadedErrorBar(1:Nmax, nanmean(cat(1,state(1).shuf(2).vel.I_mu),1), nansem(cat(1,state(1).shuf(2).vel.I_mu),1),'lineProps','k:')
shadedErrorBar(1:Nmax, nanmean(cat(1,state(2).shuf(1).vel.I_mu),1), nansem(cat(1,state(1).shuf(1).vel.I_mu),1),'lineProps','r')
shadedErrorBar(1:Nmax, nanmean(cat(1,state(2).shuf(2).vel.I_mu),1), nansem(cat(1,state(1).shuf(2).vel.I_mu),1),'lineProps','r:')


% get ratios of shuffled vs normal info

allStatNorm = cat(1,state(1).shuf(1).vel.I_mu);
allStatShuf = cat(1,state(1).shuf(2).vel.I_mu);

allRunNorm = cat(1,state(2).shuf(1).vel.I_mu);
allRunShuf = cat(1,state(2).shuf(2).vel.I_mu);

% figure, hold on
% plot(allStatNorm(:,Nmax), allStatShuf(:,Nmax),'ko')
% plot(allRunNorm(:,Nmax), allRunShuf(:,Nmax),'ro')
% plot([0 10], [0 10], 'k:')

%
allStatNorm(allStatNorm<0) = nan;
allStatShuf(allStatShuf<0) = nan;
allRunNorm(allRunNorm<0) = nan;
allRunShuf(allRunShuf<0) = nan;

statRat = allStatShuf./allStatNorm;
locoRat = allRunShuf./allRunNorm;

% figure, hold on
% shadedErrorBar(1:Nmax, nanmean(statRat,1), nansem(statRat,1),'lineProps','k')
% shadedErrorBar(1:Nmax, nanmean(locoRat,1), nansem(locoRat,1),'lineProps','r')
% 



s(isession).allStatNorm = cat(1,state(1).shuf(1).vel.I_mu);
s(isession).allStatShuf = cat(1,state(1).shuf(2).vel.I_mu);

s(isession).allRunNorm = cat(1,state(2).shuf(1).vel.I_mu);
s(isession).allRunShuf = cat(1,state(2).shuf(2).vel.I_mu);


s(isession).statRat = allStatShuf./allStatNorm;
s(isession).locoRat = allRunShuf./allRunNorm;

end

%%


% combine results and plot

figure
shadedErrorBar(1:Nmax, nanmean(cat(1,s.allStatNorm),1), nansem(cat(1,s.allStatNorm),1),'lineProps','k')
shadedErrorBar(1:Nmax, nanmean(cat(1,s.allStatShuf),1), nansem(cat(1,s.allStatShuf),1),'lineProps','k:')

shadedErrorBar(1:Nmax, nanmean(cat(1,s.allRunNorm),1), nansem(cat(1,s.allRunNorm),1),'lineProps','r')
shadedErrorBar(1:Nmax, nanmean(cat(1,s.allRunShuf),1), nansem(cat(1,s.allRunShuf),1),'lineProps','r:')

figure
shadedErrorBar(1:Nmax, nanmean(cat(1,s.statRat),1), nansem(cat(1,s.statRat),1),'lineProps','k')
shadedErrorBar(1:Nmax, nanmean(cat(1,s.locoRat),1), nansem(cat(1,s.locoRat),1),'lineProps','r')



% plot stat vs run fisher info

allStatNorm = cat(1,s.allStatNorm);
allRunNorm = cat(1,s.allRunNorm);

figure, hold on
plot(allStatNorm(:,25), allRunNorm(:,25),'ko')
plot([0 10], [0 10],'k:')
xlim([0 10]), ylim([0 10])



end

