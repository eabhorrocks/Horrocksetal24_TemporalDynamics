function [pout, info] = spikeTimesPSTH(spikeTimes,intervalsCell,binWidth,options)


%% Input vars
% - spikeTimes is vector of spike times for the unit of interest
% - intervalsCell is a cell array of arrays (size: nTrials x 2) where the
% first column is the start time and 2nd column is end time of the period
% of interest.
% - binWidth is the window over which spikes are summed
% - options is a struct which contains extra options for smoothing,
% altering the trial-window, and plotting options.


%% Example usage
% load('pro_session_821695405.mat')
% dmtrials = trials(strcmp([trials.stimulus_name],'dot_motion'));
% spikeTimes = units(104).spiketimes;
% speeds = unique([dmtrials.Speed])
% for ispeed = 1:numel(speeds)
% t = dmtrials([dmtrials.Speed]==speeds(ispeed));
% sints{ispeed}= [vertcat(t.start_time), vertcat(t.stop_time)];
% end
% binWidth = 0.01;
% options.smoothWidth = 10;
% options.cols = jet(7);
% [PSTH, info] = spikeTimesPSTH(spikeTimes,sints,binWidth, options)


%% options

if ~exist('options','var'),                 options=struct;                 end
if ~isfield(options,'smoothWidth'),         options.smoothWidth=1;          end
if ~isfield(options,'maxSmooth'),           options.maxSmooth = inf;        end
if ~isfield(options,'smoothType'),          options.smoothType='gaussian';  end
if ~isfield(options,'preTime'),             options.preTime=500;            end
if ~isfield(options,'postTime'),            options.postTime=500;           end
if ~isfield(options,'preUseTime'),          options.preUseTime=500;         end
if ~isfield(options,'postUseTime'),         options.postUseTime=500;        end
if ~isfield(options,'psthSize'),            options.psthSize=1;             end
if ~isfield(options,'plot'),                options.plot=true;              end
if ~isfield(options,'plotFRrange'),         options.plotFRrange = [];       end
if ~isfield(options,'plot95CI'),            options.plot95CI=false;         end
if ~isfield(options,'plotNTrials'),         options.plotNTrials=false;      end
if ~isfield(options,'cols'),                options.cols = lines;           end
if ~isfield(options,'patchCol'),            options.patchCol = [.9,.9,.9];  end

if ~isfield(options,'respThresh'),          options.respThresh=5;           end 
if ~isfield(options,'psthRangeThresh'),     options.psthRangeThresh=5;      end

if ~isfield(options,'getReliability'),      options.getReliability = false; end
if ~isfield(options,'maxStretch'),          options.maxStretch = 20;     end
if ~isfield(options,'nReliPerms'),          options.nReliPerms = 10;        end
if ~isfield(options,'nShuffle'),            options.nShuffle = 10;          end
if ~isfield(options,'getOnsetOffsetReli'),  options.getOnsetOffsetReli = false;          end
if ~isfield(options,'onset_idx'),           options.onset_idx = 1;          end
if ~isfield(options,'offset_idx'),          options.offset_idx = 1;          end



%%
numPSTH = numel(intervalsCell);
allints = vertcat(intervalsCell{:});
nTrials = size(allints,1);
hz_convert = 1000/binWidth;
% initalise output fields
pout = struct; info = struct;
ipsthVals = 1:numPSTH;
[pout, info] = initOutputFields(ipsthVals, pout, info);

if nTrials==0 % if no cells have an interval...
    return % finish function
end

trialLength = round(mean(allints(:,2)-allints(:,1)),2); % assumes equal length trials
trialLength = trialLength;

binedges = -1*options.preTime:binWidth:trialLength+options.postTime;
binMidPoints = movmean(binedges,2,'Endpoints', 'discard');
baselineStart = find(binMidPoints>-options.preUseTime,1,'first');
baselineEnd = find(binMidPoints<0,1, 'last');

activeStart = baselineEnd+1;
activeEnd = find(binMidPoints<trialLength+options.postUseTime, 1, 'last');

startPoint = 0;
cols = options.cols;

for ipsth = 1:numPSTH
    %% compute PSTH
    
    % get set of intervals for this PSTH
    tints = intervalsCell{ipsth};
    if isempty(tints)
        continue
    end
    
    
    stimOnTimes = tints(:,1);
    stimOffTimes = tints(:,2);
    meanDuration = mean(stimOffTimes-stimOnTimes);
    intStarts = stimOnTimes-options.preTime;
    intStops = stimOffTimes+options.postTime;
    p(ipsth).intervals = [intStarts intStops];
    
    % loop through intervals (usually trials) and bin spikes
    p(ipsth).t(size(p(ipsth).intervals,1)).spikes = [];
    for iint = 1:size(p(ipsth).intervals,1)
        % spikes that occur during this interval
        p(ipsth).t(iint).rawspikes = spikeTimes(spikeTimes>=p(ipsth).intervals(iint,1) &...
            spikeTimes <= p(ipsth).intervals(iint,2));
        p(ipsth).t(iint).nSpikes = numel(p(ipsth).t(iint).rawspikes);
        % interval-centric spike times
        p(ipsth).t(iint).tspikes = p(ipsth).t(iint).rawspikes - p(ipsth).intervals(iint,1) - options.preTime;
        % spikes binned into predefined binedges
        p(ipsth).t(iint).binSpikes = histcounts(p(ipsth).t(iint).tspikes, binedges);
        p(ipsth).t(iint).nBaseLineSpikes = ...
            sum(p(ipsth).t(iint).binSpikes(baselineStart:baselineEnd));
        p(ipsth).t(iint).baselineFR = p(ipsth).t(iint).nBaseLineSpikes/options.preUseTime * 1000;
        p(ipsth).t(iint).nActiveSpikes =...
            sum(p(ipsth).t(iint).binSpikes(activeStart:activeEnd));
        p(ipsth).t(iint).activeFR =  p(ipsth).t(iint).nActiveSpikes/trialLength+options.postUseTime * 1000;
    end
    
    
    if numel([p(ipsth).t(:).tspikes])<2 % if less than 2 spikes, use gauss smoothing
        options.smoothType = 'gaussian';
        swidth = options.maxSmooth;
    end
    if strcmp(options.smoothType,'ssv') % adaptive smoothing
        normCoeff = nanmean(cellfun(@numel, {p(ipsth).t(:).tspikes})); % undo ssv normalisation
        [y,~,optw] = ssvkernel([p(ipsth).t(:).tspikes], binMidPoints);
        PSTH_temp = y*normCoeff;
    elseif strcmp(options.smoothType,'ss')
        %normCoeff = nanmean(cellfun(@numel, {p(ipsth).t(:).tspikes}));
        [~,~,optw] = sskernel([p(ipsth).t(:).tspikes], binMidPoints); % get smth width
        swidth = optw*10; % in ms
        if swidth > options.maxSmooth
            swidth = options.maxSmooth;
        end
        
        % do normal gauss smoothing
        meanBinSpikes = mean(vertcat(p(ipsth).t.binSpikes),1);
        smthBinSpikes = smoothdata(meanBinSpikes, 'gaussian', swidth/binWidth);
        PSTH_temp = smthBinSpikes .* hz_convert; % convert to Hz fr
        
    else
        % get mean for each timepoint, over all intervals, smooth if req'd.
        meanBinSpikes = mean(vertcat(p(ipsth).t.binSpikes),1);
        smthBinSpikes = smoothdata(meanBinSpikes, options.smoothType, options.smoothWidth/binWidth);
        PSTH_temp = smthBinSpikes .* hz_convert; % convert to Hz fr
        optw = nan;
        swidth = options.smoothWidth;
    end
    
    if any(isnan(PSTH_temp)) % if ssv/ss method failed
        meanBinSpikes = mean(vertcat(p(ipsth).t.binSpikes),1);
        smthBinSpikes = smoothdata(meanBinSpikes, 'gaussian', options.smoothWidth);
        PSTH_temp = smthBinSpikes .* hz_convert; % convert to Hz fr
        optw = nan;
    end
    
    %PSTH_start_idx = find(binMidPoints>-options.preUseTime,1,'first');
    %PSTH_end_idx = find(binMidPoints<meanDuration+options.postUseTime,1,'last');
    PSTH=PSTH_temp(baselineStart:activeEnd);
    
    
    if numel(optw)>1
        pout(ipsth).optw = optw(baselineStart:activeEnd);
    else
        pout(ipsth).optw = optw;
    end
    %try
    pout(ipsth).swidth = swidth;
    info(ipsth).swidth = swidth;
    %catch
    debug=1;
    %end
    %minfr = floor(min(PSTH)); % get min fr
    p(ipsth).psth = PSTH;
    pout(ipsth).psth = PSTH;
    pout(ipsth).zPSTH = nan;
    pout(ipsth).BMP = binMidPoints(baselineStart:activeEnd);
    
end

%% Metrics

% update binMidPoints to encompass 'use' time, update event bin indexes
% temp fix

binMidPoints = binMidPoints(baselineStart:activeEnd);
baselineStart = 1;
baselineEnd = find(binMidPoints<0,1, 'last');
activeStart = baselineEnd+1;
activeEnd = numel(binMidPoints);

for ipsth = 1:numPSTH
    if all(isnan(pout(ipsth).psth))
        continue
    end
    % initialise info struct fields
    info(ipsth).nTrials = size(p(ipsth).intervals,1);
    info(ipsth).responsive = false;
    info(ipsth).zresponsive = false;
    info(ipsth).excited = false;
    info(ipsth).inhibited = false;
    info(ipsth).mixed = false;
    
    info(ipsth).baselineMean = nan;
    info(ipsth).baselineSigma = nan;
    
    info(ipsth).PSTHrange = nan;
    info(ipsth).peakSign = nan;
    info(ipsth).susIdx = nan;
    info(ipsth).onsetLatency = nan;
    info(ipsth).offsetLatency = nan;
    info(ipsth).peakLatency = nan;
    info(ipsth).peak90latency = nan;
    
    info(ipsth).dtw_dist = nan;
    info(ipsth).reliability_p = nan;
    info(ipsth).reliability_trial_p = nan;
    
    PSTH = pout(ipsth).psth; % assign PSTH using ipsth
    while numel(PSTH) < 200
        PSTH(end+1) = PSTH(end);
    end
    
    %%% baseline firing rate (mean and sigma)
    info(ipsth).meanActiveFR = mean(pout(ipsth).psth(activeStart:activeEnd));
    info(ipsth).minActiveFR = min(pout(ipsth).psth(activeStart:activeEnd));
    info(ipsth).maxActiveFR = max(pout(ipsth).psth(activeStart:activeEnd));
    
    info(ipsth).minBaselineFR = min(PSTH(baselineStart:baselineEnd));
    info(ipsth).maxBaselineFR = max(PSTH(baselineStart:baselineEnd));
    
    info(ipsth).PSTHrange = max(PSTH)-min(PSTH); % dynamic range of PSTH
    
    
    info(ipsth).onsetRange = max(PSTH(baselineEnd+options.onset_idx))-min(PSTH(baselineEnd+options.onset_idx));
   % try
    info(ipsth).offsetRange = max(PSTH(baselineEnd+options.offset_idx))-min(PSTH(baselineEnd+options.offset_idx));
    %catch
    %        info(ipsth).offsetRange = max(PSTH(baselineEnd+options.offset_idx))-min(PSTH(baselineEnd+options.offset_idx(1):end));

    %end
    
    
    %%% z-scoring %%%
    % z-score the baseline, get the normal dist, use to z-score the PSTH
    [baseline_z, b_mu, b_sig] = zscore(PSTH(baselineStart:baselineEnd));
    info(ipsth).baselineMean = b_mu;
    info(ipsth).baselineSigma = b_sig;
    PSTH_stim = PSTH(activeStart:activeEnd);
    
    z_PSTH = z_score(PSTH,b_mu,b_sig); % z-score psth according to baseline
    pout(ipsth).zPSTH = z_PSTH;
    
    info(ipsth).maxZ = max(z_PSTH(activeStart:activeEnd));
    info(ipsth).minZ = min(z_PSTH(activeStart:activeEnd));
    
    %%% determine if a unit is responsive %%%
    % check if excited
    if info(ipsth).maxZ > options.respThresh || (info(ipsth).meanActiveFR)>2*info(ipsth).baselineMean
        info(ipsth).excited = true; % if exceeds z-threshold, 'excited'
        % if median([p(ipsth).t.activeFR]) > 2*median([p(ipsth).t.baselineFR])
        %     if median([p(ipsth).t.activeFR]) > 2
        if info(ipsth).PSTHrange >= options.psthRangeThresh
            info(ipsth).responsive = true; % only responsive if passes furher criteria
        end
        %     end
        % end
    end
    
    % check if inhibited
    if info(ipsth).minZ < -options.respThresh || info(ipsth).meanActiveFR <0.5*info(ipsth).baselineMean
        info(ipsth).inhibited = true;
        %if median([p(ipsth).t.activeFR]) < 2*median([p(ipsth).t.baselineFR])
        %    if median([p(ipsth).t.baselineFR]) > 2
        if info(ipsth).PSTHrange >= options.psthRangeThresh
            info(ipsth).responsive = true;
        end
        %    end
        %end
    end
    
    if info(ipsth).excited && info(ipsth).inhibited
        info(ipsth).mixed = true;
    end
    
    z_PSTH_abs = abs(z_PSTH);
    if max(z_PSTH_abs) >= options.respThresh
        info(ipsth).zresponsive = true;
    end
    
    
    if info(ipsth).zresponsive
        onsetIdx = find(z_PSTH_abs>options.respThresh);
        onsetIdx = onsetIdx(onsetIdx>activeStart); % only during active period...
        if ~isempty(onsetIdx)
            info(ipsth).onsetLatency = binMidPoints(onsetIdx(1));
            
            % offset latency
            offsetIdx = find(z_PSTH_abs>options.respThresh,1,'last');
            
            while offsetIdx > numel(binMidPoints)
                diffbmp = diff(binMidPoints); diffbmp = diffbmp(end);
                binMidPoints(end+1) = binMidPoints(end)+diffbmp;
            end
            
            try
            info(ipsth).offsetLatency = binMidPoints(offsetIdx);
            catch
                debug = 1;
            end
        else
            info(ipsth).onsetLatency=nan;
            info(ipsth).offsetLatency=nan;
        end
    end
    
    %%% sustainedness index (how far from baseline to peak is the mean? %%%
    if info(ipsth).offsetLatency > meanDuration % use whichever is later
        susEndIdx = offsetIdx;
    else
        susEndIdx = find(binMidPoints<meanDuration, 1, 'last');
    end
    
    % sustainedness index firing rates (different from info saved)
    meanActiveFR = mean(PSTH(activeStart:susEndIdx));
    maxActiveFR = max(PSTH(activeStart:susEndIdx));
    minActiveFR = min(PSTH(activeStart:susEndIdx));
    
    if meanActiveFR >= info(ipsth).baselineMean
        info(ipsth).susIdx = (meanActiveFR - info(ipsth).baselineMean) /...
            (maxActiveFR - info(ipsth).baselineMean);
        
    elseif meanActiveFR < info(ipsth).baselineMean
        info(ipsth).susIdx = (info(ipsth).baselineMean - meanActiveFR) /...
            (info(ipsth).baselineMean - minActiveFR);
        
    end
    
    [maxAbsZ, maxAbsZidx] = max(z_PSTH_abs(activeStart:activeEnd)); % max abs z-score (idx of "peak")
    peakFR = PSTH_stim(maxAbsZidx);
    
    info(ipsth).peakSign = sign(peakFR - info(ipsth).baselineMean);
    info(ipsth).peakZ = maxAbsZ*info(ipsth).peakSign;
    info(ipsth).peakLatency = binMidPoints(baselineEnd+maxAbsZidx);
    info(ipsth).peakFR = peakFR;
    
    
    diffFR_peak = peakFR - info(ipsth).baselineMean;
    diffFR_peak90 = diffFR_peak*0.9;
    idx90 = [];
    if info(ipsth).peakSign == 1
        idx90 = find((pout(ipsth).psth > (info(ipsth).baselineMean + diffFR_peak90)), 1, 'first');
    elseif info(ipsth).peakSign == -1
        idx90 = find((pout(ipsth).psth < (info(ipsth).baselineMean + diffFR_peak90)), 1, 'first');
    end
    if isempty(idx90)
        info(ipsth).peak90latency = nan;
    else
        info(ipsth).peak90latency = binMidPoints(idx90);
    end
    
end

psthconcat = [p.psth];
minfr = floor(min(psthconcat)); % min of all psths


%% PSTH reliability

if options.getReliability
    
    % get original timings
    binMidPoints = movmean(binedges,2,'Endpoints', 'discard');
    baselineStart = find(binMidPoints>-options.preUseTime,1,'first');
    baselineEnd = find(binMidPoints<0,1, 'last');
    
    activeStart = baselineEnd+1;
    activeEnd = find(binMidPoints<trialLength+options.postUseTime, 1, 'last');
    
    onset_idx = options.onset_idx;
    offset_idx = options.offset_idx;
    
    
    maxStretch = options.maxStretch/binWidth;
    nShuffle = options.nShuffle;
    
    
    if options.nReliPerms ~= 0 % if not doing LOOCV, do half trial splits
        % split trials in  half, get means and compute distance. Shuffle
        % one of the splits and make the same comparison
        
        for ipsth = 1:numPSTH
            
            if info(ipsth).nTrials >= 10
                
                distArray = nan(options.nReliPerms,1);
                distArray_onset = distArray;
                distArray_offset = distArray;
                
                distArray_shuffle = nan(options.nReliPerms,nShuffle);
                distArray_shuffle_onset = nan(options.nReliPerms,nShuffle);
                distArray_shuffle_offset = nan(options.nReliPerms,nShuffle);
                
                swidth = pout(ipsth).swidth;
                binSwidth = swidth/binWidth;
                for iperm = 1:options.nReliPerms
                    
                    % split intervals into two
                    idx_1 = randsample(info(ipsth).nTrials, floor(info(ipsth).nTrials/2));
                    idx_2 = find(~ismember(1:info(ipsth).nTrials,idx_1))';
                    
                    meanBinSpikes = mean(vertcat(p(ipsth).t(idx_1).binSpikes),1);
                    smthBinSpikes = smoothdata(meanBinSpikes, 'gaussian', binSwidth);
                    psth_1temp = smthBinSpikes .* hz_convert; % convert to Hz fr
                    psth_1 = psth_1temp(activeStart:activeEnd);
                    % psth 2
                    meanBinSpikes = mean(vertcat(p(ipsth).t(idx_2).binSpikes),1);
                    smthBinSpikes = smoothdata(meanBinSpikes,'gaussian', binSwidth);
                    psth_2temp = smthBinSpikes .* hz_convert; % convert to Hz fr
                    psth_2 = psth_2temp(activeStart:activeEnd);
                    
                    [distArray(iperm)] = dtw(psth_1, psth_2, maxStretch);
                    
                    
                    % shuffling psth 2
                   % binSpikes = vertcat(p(ipsth).t(idx_2).binSpikes);
                   % binSpikes = binSpikes(:,activeStart:activeEnd);
                   % [nTrials, nBins] = size(binSpikes);
                   % binSpikes_shuffled = binSpikes; % pre-assign
                   
                   % get total number of spikes fired in all idx_2 trials during active period                     
                   totalUseTime = options.postUseTime+floor(trialLength);
                   totalSpikes = sum([p(ipsth).t(idx_2).tspikes]>=0 & [p(ipsth).t(idx_2).tspikes]<=totalUseTime);    
                   binedges_shuffle = binedges;
                   binedges_shuffle(binedges_shuffle<0) = [];
                    
                    for ishuffle = 1:options.nShuffle
                      %  for itrial = 1:nTrials
                      %      binSpikes_shuffled(itrial,:) = binSpikes(itrial,randperm(nBins,nBins));
                      %  end
                        
                       % get totalSpikes random spike times and bin them
                        randSpikes = sort(rand(1,totalSpikes)*totalUseTime);
                        randSpikes_binned = histcounts(randSpikes, binedges_shuffle);
                        
                        % gen psth from shuffled spikes
                        meanBinSpikes_shuffled = randSpikes_binned/numel(idx_2);
                        smthBinSpikes = smoothdata(meanBinSpikes_shuffled,'gaussian', binSwidth);
                        psth_shuffled = smthBinSpikes .* hz_convert; % convert to Hz fr
                        
                        distArray_shuffle(iperm,ishuffle) = dtw(psth_1, psth_shuffled, maxStretch);
                    end
                    
                    
                    if options.getOnsetOffsetReli
                        
                         % get total number of spikes fired in all idx_2 trials during active period                     
                        totalUseTime_onset = binWidth*numel(onset_idx);
                        totalUseTime_offset = binWidth*numel(offset_idx);
                        
                        totalSpikes_onset = sum([p(ipsth).t(idx_2).tspikes]>=binedges(baselineEnd+onset_idx(1)) & [p(ipsth).t(idx_2).tspikes]<=binedges(baselineEnd+onset_idx(end)+1));
                        totalSpikes_offset = sum([p(ipsth).t(idx_2).tspikes]>=binedges(baselineEnd+offset_idx(1)) & [p(ipsth).t(idx_2).tspikes]<=binedges(baselineEnd+offset_idx(end)+1));

                        
                        [distArray_onset(iperm)] = dtw(psth_1(options.onset_idx), psth_2(options.onset_idx), maxStretch);
                        [distArray_offset(iperm)] = dtw(psth_1(options.offset_idx), psth_2(options.offset_idx), maxStretch);
                        
                        % shuffling psth 2 ONSET
                      %  binSpikes = vertcat(p(ipsth).t(idx_2).binSpikes);
                      %  binSpikes = binSpikes(:,baselineEnd+onset_idx);
                      %  [nTrials, nBins] = size(binSpikes);
                      %  binSpikes_shuffled = binSpikes; % pre-assign
                        
                        for ishuffle = 1:options.nShuffle
                            %for itrial = 1:nTrials
                            %    binSpikes_shuffled(itrial,:) = binSpikes(itrial,randperm(nBins,nBins));
                            %end
                            
                            randSpikes = sort(rand(1,totalSpikes_onset)*totalUseTime_onset);
                            randSpikes_binned = histcounts(randSpikes, 0:binWidth:totalUseTime_onset);
                            meanBinSpikes_shuffled = randSpikes_binned/numel(idx_2);
                            
                            % gen psth from shuffled spikes
                            %meanBinSpikes_shuffled = mean(binSpikes_shuffled,1);
                            smthBinSpikes = smoothdata(meanBinSpikes_shuffled,'gaussian', binSwidth);
                            psth_shuffled = smthBinSpikes .* hz_convert; % convert to Hz fr
                            
                            distArray_shuffle_onset(iperm,ishuffle) = dtw(psth_1(onset_idx), psth_shuffled, maxStretch);
                        end
                        
                        % shuffling psth 2 OFFSET
                        %binSpikes = vertcat(p(ipsth).t(idx_2).binSpikes);
                        %binSpikes = binSpikes(:,baselineEnd+offset_idx);
                        %[nTrials, nBins] = size(binSpikes);
                        %binSpikes_shuffled = binSpikes; % pre-assign

                        
                        for ishuffle = 1:options.nShuffle
                           % for itrial = 1:nTrials
                           %     binSpikes_shuffled(itrial,:) = binSpikes(itrial,randperm(nBins,nBins));
                           % end
                            
                            randSpikes = sort(rand(1,totalSpikes_offset)*totalUseTime_offset);
                            randSpikes_binned = histcounts(randSpikes, 0:binWidth:totalUseTime_offset);
                            meanBinSpikes_shuffled = randSpikes_binned/numel(idx_2);
                            
                            % gen psth from shuffled spikes
                            %meanBinSpikes_shuffled = mean(binSpikes_shuffled,1);
                            smthBinSpikes = smoothdata(meanBinSpikes_shuffled,'gaussian', binSwidth);
                            psth_shuffled = smthBinSpikes .* hz_convert; % convert to Hz fr
                            
                            distArray_shuffle_offset(iperm,ishuffle) = dtw(psth_1(offset_idx), psth_shuffled, maxStretch);
                        end
                        
                    end
                    
                end
                % compute normalized dtw distance from distribution of distances
                info(ipsth).dtw_dist = mean(distArray) / (numel(psth_1) * info(ipsth).PSTHrange);
                info(ipsth).reliability_p = prop(distArray_shuffle(:)<=mean(distArray));
                [~, shufdist_mu, shufdist_sig] = zscore(distArray_shuffle(:));
                info(ipsth).reliability_z = z_score(mean(distArray(:)), shufdist_mu, shufdist_sig);
                
                if options.getOnsetOffsetReli
                    info(ipsth).dtw_dist_onset = mean(distArray_onset) / (numel(onset_idx) * info(ipsth).onsetRange);
                    info(ipsth).reliability_p_onset = prop(distArray_shuffle_onset(:)<=mean(distArray_onset));
                    
                    info(ipsth).dtw_dist_offset = mean(distArray_offset) / (numel(offset_idx) * info(ipsth).offsetRange);
                    info(ipsth).reliability_p_offset = prop(distArray_shuffle_offset(:)<=mean(distArray_offset));
                    
                end
                
            end
            
        end
        
    elseif options.nReliPerms == 0 % LOO-CV version
        % take 1 trial and compare to mean of the rest. shuffle that single
        % trial and make comparison to the rest to determine how reliabile
        % that single trial is.
        
        for ipsth = 1:numPSTH
            
            if info(ipsth).nTrials >= 10
                
                %distArray = nan(options.nReliPerms,1);
                %distArray_shuffle = nan(options.nReliPerms,nShuffle);
                
                distArray = nan(info(ipsth).nTrials,1);
                distArray_onset = distArray;
                distArray_offset = distArray;
                
                distArray_shuffle = nan(info(ipsth).nTrials,nShuffle);
                distArray_shuffle_onset = distArray_shuffle;
                distArray_shuffle_offset = distArray_shuffle;
                
                trial_pValues = nan(info(ipsth).nTrials,1);
                trial_pValues_onset = trial_pValues;
                trial_pValues_offset = trial_pValues;
                
                swidth = pout(ipsth).swidth; % smooth width in ms
                binSwidth = swidth/binWidth; % smooth width in bin space
                for iperm = 1:info(ipsth).nTrials
                    
                    idx_2 = iperm; % single trial
                    idx_1 = find(~ismember(1:info(ipsth).nTrials,idx_2))'; % rest of the trials
                    
                    
                    meanBinSpikes = mean(vertcat(p(ipsth).t(idx_1).binSpikes),1);
                    smthBinSpikes = smoothdata(meanBinSpikes, 'gaussian', binSwidth);
                    psth_1temp = smthBinSpikes .* hz_convert; % convert to Hz fr
                    psth_1 = psth_1temp(activeStart:activeEnd);
                    % psth 2
                    meanBinSpikes = mean(vertcat(p(ipsth).t(idx_2).binSpikes),1);
                    smthBinSpikes = smoothdata(meanBinSpikes,'gaussian', binSwidth);
                    psth_2temp = smthBinSpikes .* hz_convert; % convert to Hz fr
                    psth_2 = psth_2temp(activeStart:activeEnd);
                    
                    distArray(iperm) = dtw(psth_1, psth_2, maxStretch);
                    
                    
                     % shuffling psth 2
                   % binSpikes = vertcat(p(ipsth).t(idx_2).binSpikes);
                   % binSpikes = binSpikes(:,activeStart:activeEnd);
                   % [nTrials, nBins] = size(binSpikes);
                   % binSpikes_shuffled = binSpikes; % pre-assign
                   
                   % get total number of spikes fired in all idx_2 trials during active period                     
                   totalUseTime = options.postUseTime+floor(trialLength);
                   totalSpikes = sum([p(ipsth).t(idx_2).tspikes]>=0 & [p(ipsth).t(idx_2).tspikes]<=totalUseTime); 
                   
                   binedges_shuffle = binedges(binedges>=0);
                   
                    
                    for ishuffle = 1:options.nShuffle
                      %  for itrial = 1:nTrials
                      %      binSpikes_shuffled(itrial,:) = binSpikes(itrial,randperm(nBins,nBins));
                      %  end
                        
                       % get totalSpikes random spike times and bin them
                        randSpikes = sort(rand(1,totalSpikes)*totalUseTime);
                        randSpikes_binned = histcounts(randSpikes, binedges_shuffle);
                        
                        % gen psth from shuffled spikes
                        meanBinSpikes_shuffled = randSpikes_binned/numel(idx_2);
                        smthBinSpikes = smoothdata(meanBinSpikes_shuffled,'gaussian', binSwidth);
                        psth_shuffled = smthBinSpikes .* hz_convert; % convert to Hz fr
                        
                        distArray_shuffle(iperm,ishuffle) = dtw(psth_1, psth_shuffled, maxStretch);
                    end
                    
                    
                    if options.getOnsetOffsetReli
                        
                         % get total number of spikes fired in all idx_2 trials during active period                     
                        totalUseTime_onset = binWidth*numel(onset_idx);
                        totalUseTime_offset = binWidth*numel(offset_idx);
                        
                        totalSpikes_onset = sum([p(ipsth).t(idx_2).tspikes]>=binedges(baselineEnd+onset_idx(1)) & [p(ipsth).t(idx_2).tspikes]<=binedges(baselineEnd+onset_idx(end)+1));
                        totalSpikes_offset = sum([p(ipsth).t(idx_2).tspikes]>=binedges(baselineEnd+offset_idx(1)) & [p(ipsth).t(idx_2).tspikes]<=binedges(baselineEnd+offset_idx(end)+1));

                        
                        [distArray_onset(iperm)] = dtw(psth_1(options.onset_idx), psth_2(options.onset_idx), maxStretch);
                        [distArray_offset(iperm)] = dtw(psth_1(options.offset_idx), psth_2(options.offset_idx), maxStretch);
                        
                        % shuffling psth 2 ONSET
                      %  binSpikes = vertcat(p(ipsth).t(idx_2).binSpikes);
                      %  binSpikes = binSpikes(:,baselineEnd+onset_idx);
                      %  [nTrials, nBins] = size(binSpikes);
                      %  binSpikes_shuffled = binSpikes; % pre-assign
                        
                        for ishuffle = 1:options.nShuffle
                            %for itrial = 1:nTrials
                            %    binSpikes_shuffled(itrial,:) = binSpikes(itrial,randperm(nBins,nBins));
                            %end
                            
                            randSpikes = sort(rand(1,totalSpikes_onset)*totalUseTime_onset);
                            randSpikes_binned = histcounts(randSpikes, 0:binWidth:totalUseTime_onset);
                            meanBinSpikes_shuffled = randSpikes_binned/numel(idx_2);
                            
                            % gen psth from shuffled spikes
                            %meanBinSpikes_shuffled = mean(binSpikes_shuffled,1);
                            smthBinSpikes = smoothdata(meanBinSpikes_shuffled,'gaussian', binSwidth);
                            psth_shuffled = smthBinSpikes .* hz_convert; % convert to Hz fr
                            
                            distArray_shuffle_onset(iperm,ishuffle) = dtw(psth_1(onset_idx), psth_shuffled, maxStretch);
                        end
                        
                        % shuffling psth 2 OFFSET
                        %binSpikes = vertcat(p(ipsth).t(idx_2).binSpikes);
                        %binSpikes = binSpikes(:,baselineEnd+offset_idx);
                        %[nTrials, nBins] = size(binSpikes);
                        %binSpikes_shuffled = binSpikes; % pre-assign

                        
                        for ishuffle = 1:options.nShuffle
                           % for itrial = 1:nTrials
                           %     binSpikes_shuffled(itrial,:) = binSpikes(itrial,randperm(nBins,nBins));
                           % end
                            
                            randSpikes = sort(rand(1,totalSpikes_offset)*totalUseTime_offset);
                            randSpikes_binned = histcounts(randSpikes, 0:binWidth:totalUseTime_offset);
                            meanBinSpikes_shuffled = randSpikes_binned/numel(idx_2);
                            
                            % gen psth from shuffled spikes
                            %meanBinSpikes_shuffled = mean(binSpikes_shuffled,1);
                            smthBinSpikes = smoothdata(meanBinSpikes_shuffled,'gaussian', binSwidth);
                            psth_shuffled = smthBinSpikes .* hz_convert; % convert to Hz fr
                            
                            distArray_shuffle_offset(iperm,ishuffle) = dtw(psth_1(offset_idx), psth_shuffled, maxStretch);
                        end
                        
                    end
                    
                    [~, shufdist_mu_trial, shufdist_sig_trial] = zscore(distArray_shuffle(iperm,:));
                    trial_zscore(iperm) = z_score(distArray(iperm), shufdist_mu_trial, shufdist_sig_trial);
                    trial_pValues(iperm) = prop(distArray_shuffle(iperm,:)<=distArray(iperm));
                    trial_pValues_onset(iperm) = prop(distArray_shuffle_onset(iperm,:)<=distArray_onset(iperm));
                    trial_pValues_offset(iperm) = prop(distArray_shuffle_offset(iperm,:)<=distArray_offset(iperm));
                end
                
                [~, shufdist_mu, shufdist_sig] = zscore(distArray_shuffle(:));
                info(ipsth).reliability_z = z_score(mean(distArray(:)), shufdist_mu, shufdist_sig);
                info(ipsth).reliability_z_trials = trial_zscore;
                
                reli_pValue = prop(distArray_shuffle(:)<=mean(distArray));
                % compute normalized dtw distance from distribution of distances
                info(ipsth).dtw_dist = mean(distArray) / (numel(psth_1) * info(ipsth).PSTHrange);
                info(ipsth).reliability_p = reli_pValue;
                info(ipsth).reliability_trial_p = trial_pValues;
                
                
                if options.getOnsetOffsetReli
                    info(ipsth).dtw_dist_onset = mean(distArray_onset) / (numel(onset_idx) * info(ipsth).onsetRange);
                    info(ipsth).reliability_p_onset = prop(distArray_shuffle_onset(:)<=mean(distArray_onset));
                    info(ipsth).reliability_trial_p_onset = trial_pValues_onset;
                    
                    info(ipsth).dtw_dist_offset = mean(distArray_offset) / (numel(offset_idx) * info(ipsth).offsetRange);
                    info(ipsth).reliability_p_offset = prop(distArray_shuffle_offset(:)<=mean(distArray_offset));
                    info(ipsth).reliability_trial_p_offset = trial_pValues_offset;
                    
                end
            end
            
        end
        
    end
    
    
    
    
end


%% plotting

if options.plot
    %figure,
    hold on
    % scale PSTHs for plotting
    nTrials = size(allints,1);
    desiredSize = nTrials*options.psthSize;
    gap = ceil(nTrials*0.05); % 5% of ntrials gap...
    PSTHstart = nTrials + gap;
    
    % default plotting scale (this all needs a clean up!)
    if isempty(options.plotFRrange)
        PSTHrange = ceil(max(psthconcat))-minfr;
        scalingFactor = desiredSize/PSTHrange;
        
        
        for ipsth = 1:numPSTH
            p(ipsth).PSTHscaled = (pout(ipsth).psth-minfr)*scalingFactor;
        end
        psthconcatscaled = [p.PSTHscaled];
        maxPt = ceil(PSTHstart+max(psthconcatscaled));
        
        
        %%%% temp fix for user input firing rate range... %%%%
    else % custom FR range set
        if numel(options.plotFRrange)~=2
            error('options.plotFRrange must be [minFR, maxFR]')
        end
        minInputFR = options.plotFRrange(1);
        maxInputFR = options.plotFRrange(2);
        minfr = minInputFR;
        PSTHrange = maxInputFR-minInputFR;
        
        scalingFactor = desiredSize/PSTHrange;
        for ipsth = 1:numPSTH
            p(ipsth).PSTHscaled = (pout(ipsth).psth-minfr)*scalingFactor;
        end
        
        psthconcatscaled = [p.PSTHscaled];
        maxPt = ceil(PSTHstart+maxInputFR*scalingFactor);
    end
    
    % trial reference patch
    v = [0 0.2; trialLength 0.2; trialLength maxPt+gap; 0 maxPt+gap];
    f = [1 2 3 4];
    patch('Faces',f,'Vertices',v,...
        'EdgeColor','none','FaceColor',options.patchCol,'LineWidth',2);
    % plot the PSTHs
    for ipsth = 1:numPSTH
        try
            plot(pout(ipsth).BMP,PSTHstart+p(ipsth).PSTHscaled, '-', 'LineWidth', 1.5,...
                'Color', cols(ipsth,:))
        catch
            debug=1;
        end
    end
    %     if options.plot95CI
    %
    %         stdBinSpikes = std(vertcat(t.binSpikes),1);
    %         smthStdSpikes = smoothdata(stdBinSpikes, options.smoothType, options.smoothWidth);
    %         scaledSTD = smthStdSpikes*scalingFactor;
    %         CI95 = scaledSTD*1.96;
    %         maxPt = ceil(PSTHstart+max(PSTHscaled+CI95));
    %         plot([0 0], [0, maxPt+gap], 'r-', 'LineWidth', 1)
    %         plot([trialLength trialLength], [0, maxPt+gap], 'r-', 'LineWidth', 1)
    %         plot 95% CI
    %         plotCIs(PSTHstart, PSTHscaled, CI95, binMidPoints);
    %         plot(binMidPoints,PSTHstart+PSTHscaled+CI95, 'Color', [.8 .8 .8], 'LineWidth', 2);
    %         plot(binMidPoints,PSTHstart+PSTHscaled-CI95, 'Color', [.8 .8 .8], 'LineWidth', 2);
    %         plot scaled PSTH
    %         plot(binMidPoints,PSTHstart+PSTHscaled, 'k-', 'LineWidth', 1);
    %     else % just mean
    %    end
    
    % plot spike rasters as dots
    startPoint = 0;
    for ipsth = 1:numPSTH
        for iint = 1:size(p(ipsth).intervals,1)
            plot(p(ipsth).t(iint).tspikes, startPoint+repelem(iint, p(ipsth).t(iint).nSpikes),...
                '.', 'MarkerSize', 3, 'Color', cols(ipsth,:)*0.7)
        end
        startPoint = startPoint + size(p(ipsth).intervals,1);
    end
    
    %% set axes properties
    tickJumps = round(PSTHrange/5); % good ?
    tickJumpsAxis = tickJumps * scalingFactor;
    nTrialsVec = cellfun(@numel, {p.t});
    ax = gca;
    if options.plotNTrials
        cumsumnTrialsVec = cumsum(nTrialsVec);
        ax.YTick = [cumsumnTrialsVec, PSTHstart:tickJumpsAxis:maxPt];
        ax.YTickLabel = {nTrialsVec, minfr:tickJumps:minfr+tickJumps*numel(ax.YTick)};
    else
        ax.YTick = [PSTHstart:tickJumpsAxis:maxPt];
        ax.YTickLabel = {minfr:tickJumps:minfr+tickJumps*numel(ax.YTick)};
    end
    xlim([-options.preUseTime, trialLength+options.postUseTime])
    ylim([-1 maxPt+gap]);
    xlabel('Trial time (s)');
    ylabel('Firing Rate (Hz)')
    ylabh = get(gca,'ylabel');
    set(ylabh,'Units','normalized');
    shift = (options.psthSize/(1+options.psthSize))/2;
    set(ylabh,'position',[-0.0557 0.5+shift 0]);
end

end


function plotCIs(PSTHstart, PSTHscaled, CI95, binMidPoints)
% takes patch code from shadedErrorBar
edgeColor=[0 0 0]+(1-[0 0 0])*0.55;
faceAlpha=0.2;
patchColor=[ 0 0 0];
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
%Make pretty edges around the patch.
H.edge(1)=plot(x,lE,'-');
H.edge(2)=plot(x,uE,'-');
set([H.edge], 'color',edgeColor, ...
    'HandleVisibility','off', ...
    'Tag', 'shadedErrorBar_edge')
end




function [pout, info] = initOutputFields(ipsthVals, pout, info)

for ipsth = ipsthVals
    
    info(ipsth).nTrials = 0;
    info(ipsth).responsive = false;
    info(ipsth).zresponsive = false;
    info(ipsth).excited = false;
    info(ipsth).inhibited = false;
    info(ipsth).mixed = false;
    
    info(ipsth).baselineMean = nan;
    info(ipsth).baselineSigma = nan;
    
    info(ipsth).PSTHrange = nan;
    info(ipsth).onsetRange = nan;
    info(ipsth).offsetRange = nan;
    
    
    info(ipsth).peakSign = nan;
    info(ipsth).susIdx = nan;
    info(ipsth).onsetLatency = nan;
    info(ipsth).offsetLatency = nan;
    info(ipsth).peakLatency = nan;
    info(ipsth).peak90latency = nan;
    
    info(ipsth).swidth = nan;
    info(ipsth).dtw_dist = nan;
    info(ipsth).reliability_p = nan;
    info(ipsth).reliability_trial_p = nan;
    
    info(ipsth).dtw_dist_onset = nan;
    info(ipsth).reliability_p_onset = nan;
    info(ipsth).reliability_trial_p_onset = nan;
    
    info(ipsth).dtw_dist_offset = nan;
    info(ipsth).reliability_p_offset = nan;
    info(ipsth).reliability_trial_p_offset = nan;
    
    pout(ipsth).optw = nan;
    pout(ipsth).psth = nan;
    pout(ipsth).zPSTH = nan;
    pout(ipsth).BMP = nan;
end



end


