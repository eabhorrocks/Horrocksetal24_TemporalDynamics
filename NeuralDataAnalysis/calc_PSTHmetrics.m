function metrics = calc_PSTHmetrics(psth, zpsth, binMidPoints, baselineFR, zThresh, minRespEnd)


%% initialise metric fields
metrics.responsive = nan;
metrics.excited = nan;
metrics.suppressed = nan;
metrics.onsetLatency = nan;
metrics.offsetLatency = nan;
metrics.peakLatency = nan;

metrics.peakSign = nan;
metrics.peakAbsZ = nan;
metrics.peak90Latency = nan;
metrics.meanZ = nan;
metrics.minZ = nan;
metrics.maxZ = nan;
metrics.rangeZ = nan;
metrics.minZLat = nan;
metrics.maxZLat = nan;

metrics.maxFR = nan;
metrics.minFR = nan;
metrics.meanFR = nan;
metrics.rangeFR = nan;
metrics.susIdx = nan;


%% calculate metrics

% check responsiveness
absPSTHz = abs(zpsth);
responsive = any(absPSTHz>=zThresh);

% excited/suppressed
excited = any(zpsth>=zThresh);
suppressed = any(zpsth<=-zThresh);

% index where t=0
respStartIndex = find(binMidPoints>=0,1,'first');

%onset latency
onsetLatency = binMidPoints(find(absPSTHz(respStartIndex:end)>=zThresh,1,'first')+respStartIndex-1);
offsetLatency = binMidPoints(find(absPSTHz(respStartIndex:end)>=zThresh,1,'last')+respStartIndex-1);
respEndIndex = find(binMidPoints>=max([find(absPSTHz>=zThresh,1,'last'), minRespEnd]),1,'first');

% find z-score peak and related latencies
[peakAbsZ, maxAbsZidx] = max(absPSTHz);
peakLatency = binMidPoints(maxAbsZidx);
peakSign = sign(zpsth(maxAbsZidx));

idx90 = [];
if peakSign == 1
    idx90 = find(zpsth>=peakAbsZ*0.9, 1, 'first');
elseif peakSign == -1
    idx90 = find(zpsth<=-peakAbsZ*0.9, 1, 'first');
end
peak90Latency = binMidPoints(idx90);

% other z-score properties
meanZ = mean(zpsth);
[minZ, minZidx] = min(zpsth);
[maxZ, maxZidx] = max(zpsth);
rangeZ = maxZ-minZ;

minZLat = binMidPoints(minZidx);
maxZLat = binMidPoints(maxZidx);


% get firing rate properties
maxFR = max(psth(respStartIndex:respEndIndex));
minFR = min(psth(respStartIndex:respEndIndex));
rangeFR = maxFR-minFR;
meanFR = mean(psth(respStartIndex:respEndIndex));


% sustainedness index
if meanFR>=baselineFR
    deltaMean = meanFR-baselineFR;
    deltaPeak = maxFR-baselineFR;
elseif meanFR<baselineFR
    deltaMean = baselineFR-meanFR;
    deltaPeak = baselineFR-minFR;
else
    deltaMean = nan;
    deltaPeak = nan;
end

susIdx = deltaMean/deltaPeak;




%% assign metrics to output struct

if ~isempty(responsive), metrics.responsive = responsive; end
if ~isempty(excited), metrics.excited = excited; end
if ~isempty(suppressed), metrics.suppressed = suppressed; end
if ~isempty(onsetLatency), metrics.onsetLatency = onsetLatency; end
if ~isempty(offsetLatency), metrics.offsetLatency = offsetLatency; end
if ~isempty(peakLatency), metrics.peakLatency = peakLatency; end
if ~isempty(peak90Latency), metrics.peak90Latency = peak90Latency; end

if ~isempty(peakSign), metrics.peakSign = peakSign; end
if ~isempty(peakAbsZ), metrics.peakAbsZ = peakAbsZ; end
if ~isempty(meanZ), metrics.meanZ = meanZ; end
if ~isempty(minZ), metrics.minZ = minZ; end
if ~isempty(maxZ), metrics.maxZ = maxZ; end
if ~isempty(rangeZ), metrics.rangeZ = rangeZ; end
if ~isempty(minZLat), metrics.minZLat = minZLat; end
if ~isempty(maxZLat), metrics.maxZLat = maxZLat; end

if ~isempty(maxFR), metrics.maxFR = maxFR; end
if ~isempty(minFR), metrics.minFR = minFR; end
if ~isempty(meanFR), metrics.meanFR = meanFR; end
if ~isempty(rangeFR), metrics.rangeFR = rangeFR; end
if ~isempty(susIdx), metrics.susIdx = susIdx; end



end

