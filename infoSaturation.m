%% info saturation estimation based on kafashan et al

% to do
% - load basic.mat
% - partition trials by state (runFlag)
% - get binned spike counts
% - loop through each speed combination (just neighbouring?), w/ and w/o
% shuffling ?one set of trials


%% loada data

load('M22027_20220517_basic.mat')


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

%% get time-avaerage values for now

for iunit = 1:numel(units)
    units(iunit).avS = cellfun(@(x) mean(x,1).*100, units(iunit).allSpikes,'UniformOutput',false);
end

%%
goodUnits = units([units.isi_viol]<=0.1...
& [units.amplitude_cutoff]<=0.1 & [units.amplitude]>=50);
 

%%
% sr is spike rate (nTrial x nCell)
% ori is nTrial vector of orientations
% con is nTrial vector of contrasts (fixed for m25a at least)
% oris(ori1id), oris(ori2id) are the actual orientations (e.g 45, 90)
% subdims is number of cells (could be optionally less if wanted?)
% min([T Tb]) is inf in the case I tested (nTrials related?)


istate = 1;
shuf = false;
subdims = 20;

sr = [];
for iunit = 1:numel(goodUnits)
    allResps = [goodUnits(iunit).avS{:,istate}];
    sr(:,iunit) = allResps;
end

velVec = repelem(1:6,cellfun(@numel, goodUnits(1).avS(:,istate)));
velVec = velVec(:);

%
% for istate = 1:2
%     for ivel = 1:5
ivel = 2;
v1 = ivel;
v2 = ivel+1;

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
Nmax = min(2*T-3, N);  % avoid infinite var estimator

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
popsamples = 1000;

fprintf('Computing information estimates (N=%d, Nmax=%d) ...\n', N, Nmax);
[Iincr_mu, Iincr_var, Iincr_samples] = ...
    estIincrMoments(mu, S, T, ds, popsamples, Nmax);

% 
%     end
% end

% plot result

 % plot moments
        figure('Color', 'white');
        subplot(2, 1, 1);  hold on;
        patch([1:Nmax fliplr(1:Nmax)], ...
            [(Iincr_mu+sqrt(Iincr_var)) fliplr(Iincr_mu-sqrt(Iincr_var))],1,...
            'FaceColor', [0.8 0.8 0.8], 'EdgeColor', 'none');
        plot(1:Nmax, Iincr_mu, 'k-', 'LineWidth', 2);
        plot([1 Nmax], [0 0], 'k--');
        ylabel('Info increase estimate');

        subplot(2, 1, 2);  hold on;
        I_var = cumsum(Iincr_var);
        I_mu = cumsum(Iincr_mu);
        patch([1:Nmax fliplr(1:Nmax)], ...
            [(I_mu+sqrt(I_var)) fliplr(I_mu-sqrt(I_var))],1,...
            'FaceColor', [0.8 0.8 0.8], 'EdgeColor', 'none');
        plot(1:Nmax, I_mu, 'k-', 'LineWidth', 2);
        xlabel('N');
        ylabel('Info estimate');

