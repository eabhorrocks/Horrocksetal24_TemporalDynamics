function [psth, blank, binMidPoints] = makePSTH(spikeTimes,intervalsCell,blankIntervals,options)


%% options

if ~exist('options','var'),                 options=struct;                 end

% timing and smoothing options
if ~isfield(options,'binWidth'),            options.binWidth=10;            end
if ~isfield(options,'smoothWidth'),         options.smoothWidth=1;          end
if ~isfield(options,'smoothType'),          options.smoothType='gaussian';  end
if ~isfield(options,'preTime'),             options.preTime=500;            end
if ~isfield(options,'postTime'),            options.postTime=500;           end

% reliability options
if ~isfield(options,'getReliability'),      options.getReliability = false; end
if ~isfield(options,'distType'),            options.distType = 'euclidean'; end % any pdist arg, including custom fun
if ~isfield(options,'kfold'),               options.kfold = 2;              end % put kfold>=nTrials to do LOOCV
if ~isfield(options,'nReliPerms'),          options.nReliPerms = 10;        end % nReps for cross-val
if ~isfield(options,'nShuffle'),            options.nShuffle = 50;          end % nShuffles per k-fold cross val

% plot options
if ~isfield(options,'psthSize'),            options.psthSize=1;             end
if ~isfield(options,'plot'),                options.plot=true;              end
if ~isfield(options,'plotFRrange'),         options.plotFRrange = [];       end
if ~isfield(options,'plot95CI'),            options.plot95CI=false;         end
if ~isfield(options,'plotNTrials'),         options.plotNTrials=false;      end
if ~isfield(options,'cols'),                options.cols = lines;           end
if ~isfield(options,'patchCol'),            options.patchCol = [.9,.9,.9];  end

%% get basic info and initialise output fields

numPSTH = numel(intervalsCell); % number of stimulus conditions
allints = vertcat(intervalsCell{:}); % all intervals concatenated
nTrials = size(allints,1); % number of trials in total
hz_convert = 1000/options.binWidth; % conversion factor for PSTH in spikes/s

blankPSTH = [];
blank_mu = [];
blank_std = [];
blank.psth = blankPSTH;
blank.mu = blank_mu;
blank.std = blank_std;
binMidPoints = [];
for ipsth =1:numPSTH
    psth(ipsth).psth = [];
    psth(ipsth).zpsth = [];
end

if nTrials==0 % if there are not trials use default output args
    for ipsth =1:numPSTH
    psth(ipsth).nTrials = 0;
    end
    psth = reshape(psth,size(intervalsCell));

    return % finish function
end

trialLength = round(mean(allints(:,2)-allints(:,1)),-1,'decimals'); % round to 10ms
binedges = -1*options.preTime:options.binWidth:trialLength+options.postTime; % for discretising psth
binMidPoints = movmean(binedges,2,'Endpoints', 'discard'); % for plotting, latencies, etc.


%% use blank trials to generate spontaneous rate

% if blankIntervals provided:, otherwise use preTime baseline period later
% on
if ~isempty(blankIntervals)
tints = blankIntervals; % these intervals

stimOnTimes = tints(:,1);
stimOffTimes = tints(:,2);
intStarts = stimOnTimes-options.preTime; % add pre and post-time for interval
intStops = stimOnTimes+trialLength+options.postTime;

% loop through intervals (usually trials) and bin spikes
for itrial = 1:size(tints,1)
    
    % spikes that occur during this interval
    relevantSpikeTimes = spikeTimes(spikeTimes>=intStarts(itrial) &...
        spikeTimes <= intStops(itrial));
    % 0-centre to stimulus onset
    t(itrial).trialCentricSpikes = relevantSpikeTimes - stimOnTimes(itrial);
    % discretise spikes
    t(itrial).binnedSpikes = histcounts(t(itrial).trialCentricSpikes, binedges);
end

% do normal smoothing
meanBinSpikes = mean(vertcat(t.binnedSpikes),1);
smthBinSpikes = smoothdata(meanBinSpikes, options.smoothType, options.smoothWidth/options.binWidth);
blankPSTH = smthBinSpikes.*hz_convert;
blank_mu = mean(blankPSTH);
blank_std = std(blankPSTH);
blank.psth = blankPSTH;
blank.mu = blank_mu;
blank.std = blank_std;

end


clear t;


%% loop through stimulus conditions and create psths

for ipsth = 1:numPSTH
    % compute PSTH
    
    % get set of intervals for this PSTH
    tints = intervalsCell{ipsth};
    if isempty(tints)
            psth(ipsth).nTrials = 0;
         continue
    end
    
    stimOnTimes = tints(:,1);
    stimOffTimes = tints(:,2);
    intStarts = stimOnTimes-options.preTime;
    intStops = stimOnTimes+trialLength+options.postTime;
    
    % loop through intervals (usually trials) and bin spikes
    for itrial = 1:size(tints,1)
        
        % spikes that occur during this interval
        relevantSpikeTimes = spikeTimes(spikeTimes>=intStarts(itrial) &...
            spikeTimes <= intStops(itrial));
        % 0-centre to stimulus onset. save these for later raster plot
        p(ipsth).trial(itrial).trialCentricSpikes = relevantSpikeTimes - stimOnTimes(itrial);
        p(ipsth).trial(itrial).nSpikes = numel(p(ipsth).trial(itrial).trialCentricSpikes);
        % discretise spikes
        p(ipsth).trial(itrial).binnedSpikes = histcounts(p(ipsth).trial(itrial).trialCentricSpikes, binedges);
    end
    
    % do normal smoothing
    meanBinSpikes = mean(vertcat(p(ipsth).trial.binnedSpikes),1);
    smthBinSpikes = smoothdata(meanBinSpikes, options.smoothType, options.smoothWidth/options.binWidth);
    psth(ipsth).psth = smthBinSpikes.*hz_convert;

    % compute z-score PSTH using no-stimulus intervals (blanks or ITI)
    if ~isempty(blankIntervals) % if blank trials available
    psth(ipsth).zpsth = z_score(psth(ipsth).psth, blank_mu, blank_std);
    else % otherwise use pre-stimulus baseline period of PSTH
        baselineStart = find(binMidPoints>-options.preTime,1,'first');
        baselineEnd = find(binMidPoints<0,1, 'last');
        [baseline_z, b_mu, b_sig] = zscore(psth(ipsth).psth(baselineStart:baselineEnd));
        psth(ipsth).zpsth = z_score(psth(ipsth).psth, b_mu, b_sig);
    end

    psth(ipsth).nTrials = numel(p(ipsth).trial);
end

%% cross-val reliability of PSTH shape

% add option of getting reliability for specific periods (e.g. ignore
% pre-stim)
% use fixed number of training trials instead of k-fold


if options.getReliability
    for ipsth = 1:numel(p)
        
        nTrials = numel(p(ipsth).trial);
        if nTrials<5
            psth(ipsth).reliabilityTrial_z = nan;
            psth(ipsth).reliability_z = nan;
            continue
        end
        % check of kfold>nTrials. If so, we're doing LOOCV.
        if options.kfold>=nTrials
            options.kfold = nTrials;
            nReliPerms = 1; % there's no point repeating LOOCV.
            distArray = nan(nTrials,1); % save one dist per test trial
            distArray_shuffle = nan(options.kfold,options.nShuffle,nReliPerms); % save all the shuffles
            
        else % otherwise we're doing kfold with nReliPerm reps
            nReliPerms = options.nReliPerms;
            distArray = nan(nReliPerms,1); % save one dist per nReliPerms
            distArray_shuffle = nan(options.kfold,options.nShuffle,nReliPerms); % save all the shuffles
        end
        
        for iperm = 1:nReliPerms % loop through cv perms
            
            cvobj = cvpartition(nTrials, 'kfold', options.kfold); % generate train and test sets
            
            for ik = 1:options.kfold
                tempDist = [];
                
                % split trials into train/test
                trainIdx = training(cvobj,ik);
                testIdx = test(cvobj,ik);
                
                % make training and test psth
                trainPSTH = smoothdata(mean(vertcat(p(ipsth).trial(trainIdx).binnedSpikes).*hz_convert,1),...
                    options.smoothType, options.smoothWidth/options.binWidth);
                
                testPSTH = smoothdata(mean(vertcat(p(ipsth).trial(testIdx).binnedSpikes).*hz_convert,1),...
                    options.smoothType, options.smoothWidth/options.binWidth);
                
                tempDist = pdist([trainPSTH;testPSTH], options.distType); % get distance between trainPSTH and testPSTH
                %tempDist = dtw(trainPSTH, testPSTH, 20);
                
                % now make shuffled PSTHs from test trials and compare to trainPSTH
                tempShuffDist = [];
                for ishuffle = 1:options.nShuffle
                    totalUseTime = options.preTime+trialLength+options.postTime; % total amount of time to construct psth during
                    totalSpikes = sum([p(ipsth).trial(testIdx).nSpikes]); % number of spikes from unit summed over trials
                    randSpikes = sort(rand(1,totalSpikes)*totalUseTime)-options.preTime; % generate randomly timed spikes
                    randSpikes_binned = histcounts(randSpikes, binedges); % bin randomly timed spikes
                    shufflePSTH = smoothdata(((randSpikes_binned.*hz_convert)./sum(testIdx)),... % make averaged, smoothed PSTH
                        options.smoothType, options.smoothWidth/options.binWidth);
                    
                    tempShuffDist(ishuffle) = pdist([trainPSTH;shufflePSTH], options.distType);  % get distance betwen shuffled and trainPSTH
                    
                end
                
                if options.kfold>=nTrials % if doing LOOCV
                    distArray(testIdx) = tempDist; % save distance to relevant trial
                    distArray_shuffle(testIdx,:,iperm) = tempShuffDist; % save shuffle distances
                else
                    tempDistArray(ik) = tempDist; % if doing normal k-fold, save tempDist under this k-index
                    distArray_shuffle(ik, :, iperm) = tempShuffDist; % save shuffle distances
                end
                
            end
            
            if options.kfold<nTrials % if not doing LOOCV, average dtw dist over k-folds
                distArray(iperm) = mean(tempDistArray);
                tempDistArray=[];
            end
            
        end
        
        % z-score the mean real dist (distArray) against the distribution
        % of shuffled dists
        
        shuffle_mu = nanmean(distArray_shuffle(:));
        shuffle_std = nanstd(distArray_shuffle(:));
        psth(ipsth).reliability_z = z_score(nanmean(distArray(:)), shuffle_mu, shuffle_std);
        
        % if doing LOOCV, additionally get reliability_z for each trial
        if options.kfold>=nTrials
            psth(ipsth).reliabilityTrial_z = z_score(distArray, shuffle_mu, shuffle_std);
        else
            psth(ipsth).reliabilityTrial_z = nan(nTrials,1);
        end
    end
end

%% reshape psth into original shape of intervalsCell
psth = reshape(psth,size(intervalsCell));



%% plotting

if options.plot
    hold on
    % scale PSTHs for plotting
    nTrials = size(allints,1);
    desiredSize = nTrials*options.psthSize;
    gap = ceil(nTrials*0.05); % 5% of ntrials gap...
    PSTHstart = nTrials + gap;
    
    % get scaling factor for plotting at appropriate size
    if isempty(options.plotFRrange)
        
        minfr = floor(min([psth.psth])); % lowest fr of all psths
        maxfr = ceil(max([psth.psth])); % highest fr of all psths
    else % user-custom FR range set, e.g. to compare responses produced by different function calls
        if numel(options.plotFRrange)~=2
            error('options.plotFRrange must be [minFR, maxFR]')
        end
        minfr = options.plotFRrange(1);
        maxfr = options.plotFRrange(2);
    end
    
    PSTHrange = maxfr-minfr;
    scalingFactor = desiredSize/PSTHrange;
    
    for ipsth = 1:numPSTH
        p(ipsth).PSTHscaled = (psth(ipsth).psth-minfr)*scalingFactor;
    end
    
    psthconcatscaled = [p.PSTHscaled];
    maxPt = ceil(PSTHstart+maxfr*scalingFactor);
    
    % plot background trial reference patch
    v = [0 0.2; trialLength 0.2; trialLength maxPt+gap; 0 maxPt+gap];
    f = [1 2 3 4];
    patch('Faces',f,'Vertices',v,...
        'EdgeColor','none','FaceColor',options.patchCol,'LineWidth',2);
    
    
    % plot the PSTHs
    
    if options.plot95CI
        
        for ipsth = 1:numPSTH
            stdBinSpikes = std(vertcat(p(ipsth).trial.binnedSpikes),1);
            smthStdSpikes = smoothdata(stdBinSpikes, options.smoothType, options.smoothWidth);
            scaledSTD = smthStdSpikes*scalingFactor;
            CI95 = scaledSTD*1.96;
            %plot 95% CI
            plotCIs(PSTHstart, p(ipsth).PSTHscaled, CI95, binMidPoints, options.cols(ipsth,:))
        end
    end
    
    % plot means
    for ipsth = 1:numPSTH
        plot(binMidPoints,PSTHstart+p(ipsth).PSTHscaled, '-', 'LineWidth', 1.5,...
            'Color', options.cols(ipsth,:))
    end
    
    % plot spike rasters as dots
    startPoint = 0;
    for ipsth = 1:numPSTH
        for itrial = 1:size(intervalsCell{ipsth},1)
            plot(p(ipsth).trial(itrial).trialCentricSpikes, startPoint+repelem(itrial, p(ipsth).trial(itrial).nSpikes),...
                '.', 'MarkerSize', 3, 'Color', options.cols(ipsth,:)*0.7)
        end
        startPoint = startPoint + size(intervalsCell{ipsth},1);
    end
    
    %% set axes properties

    tickJumps = round(PSTHrange/5); % good?
    tickJumpsAxis = tickJumps * scalingFactor;
    nTrialsVec = cellfun(@numel, {p.trial});
    ax = gca;
    if options.plotNTrials
        cumsumnTrialsVec = cumsum(nTrialsVec);
        ax.YTick = [cumsumnTrialsVec, PSTHstart:tickJumpsAxis:maxPt];
        ax.YTickLabel = {nTrialsVec, minfr:tickJumps:minfr+tickJumps*numel(ax.YTick)};
    else
        ax.YTick = PSTHstart:tickJumpsAxis:maxPt;
        ax.YTickLabel = {minfr:tickJumps:minfr+tickJumps*numel(ax.YTick)};
    end
    xlim([-options.preTime, trialLength+options.postTime])
    ylim([-1 maxPt+gap]);
    xlabel('Trial time (ms)');
    ylabel('Firing Rate (Hz)')
    ylabh = get(gca,'ylabel');
    set(ylabh,'Units','normalized');
    shift = (options.psthSize/(1+options.psthSize))/2;
    set(ylabh,'position',[-0.0657 0.5+shift 0]);
    
    
end

end


function plotCIs(PSTHstart, PSTHscaled, CI95, binMidPoints, patchCol)
% adapted patch code from shadedErrorBar
faceAlpha=0.2;
patchColor=patchCol;
%Calculate the error bars
x = binMidPoints;
uE=PSTHstart+PSTHscaled+CI95;
lE=PSTHstart+PSTHscaled-CI95;
%Make the patch (the shaded error bar)
yP=[lE,fliplr(uE)];
xP=[x,fliplr(x)];
%remove nans otherwise patch won't work
xP(isnan(yP))=[];
yP(isnan(yP))=[];

H.patch=patch(xP,yP,1);
set(H.patch,'facecolor',patchColor, ...
    'edgecolor','none', ...
    'facealpha',faceAlpha, ...
    'HandleVisibility', 'off', ...
    'Tag', 'shadedErrorBar_patch')
end


