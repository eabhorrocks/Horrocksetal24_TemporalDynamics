function output = calcPopulationGeoMetrics(neuralTrajectory, weights, timeBinVector, initialSSIdx, finalSSIdx, nWins)

%% input args
% - neuralTrajectory is array of factor scores (itime, idim)
% - weights is the fraction of shared variance explained by each component
% (optional)
% - timeBinVector is the times (first edge) associated with itime indexes
% - initialSSIdx is the time indces which define the initial steady-state
% - finalSSIdx is as above, for the final steady-state
% - nWins is the number of contiguous time bins that a trajectory must be in the
% - finalSSIdx before it is considered to reach it.

%% output args
% struct with fields:
% 'maxInitialSSDist', maxmimum distance that defines the initial SS
% 'maxFinalSSDist', as above, for the final SS
% 'distFromInitialSS', distance from the initial SS for each time bin
% 'distFromFinalSS', as above, for the final SS
% 'directDistance', direct distance between initial SS and the final SS
% 'cumulativeDistanceTravelled, distnace travelled by trajectory between initial and final SS
% 'distanceRatio', the ratio of cumulative distance travelled divided by direct distance (therefore >=1)
% 'distTravelledVector', distance travelled from time bin i-1 to i, where i is the current bin
% 'exitInitialSSLatency', latency for when the trajectory leaves the initial SS
% 'enterFinalSSLatency', latency for when the trajectory reaches the final SS



%% get the start point of the trajectory

zeroIdx = initialSSIdx(end)+1;
thisTraj = neuralTrajectory(zeroIdx:finalSSIdx(end),:);

%% define the initial steady state
initialSSPoints = neuralTrajectory(initialSSIdx(1):initialSSIdx(end),:);

% calculate distances from mean point of steady state to all other points
meanInitialSSPoint = mean(initialSSPoints,1); % mean point of steady state
initialSSDists = nan(size(initialSSPoints,1),1);
for ipt = 1:size(initialSSPoints,1)
    initialSSDists(ipt) = sqrt(sum(weights.*(initialSSPoints(ipt,:)-meanInitialSSPoint).^2)); % all distances from this mean
end
maxInitialSSDist = max(initialSSDists); % max distance from mean in steady state


%% define the final steady state
finalSSPoints = neuralTrajectory(finalSSIdx(1):finalSSIdx(end),:);

% calculate distances from mean point of steady state to all other points
meanFinalSSPoint = mean(finalSSPoints,1); % mean point of steady state
finalSSDists = nan(size(finalSSPoints,1),1);
for ipt = 1:size(finalSSPoints,1)
    finalSSDists(ipt) = sqrt(sum(weights.*(finalSSPoints(ipt,:)-meanFinalSSPoint).^2)); % all distances from this mean
end
maxFinalSSDist = max(finalSSDists); % max distance from mean in steady state

%% get the distance between the trajectory and the final SS, and index of entry
distToFinalSS = nan(size(thisTraj,1),1);
for itime = 1:size(thisTraj,1)
    distToFinalSS(itime) = sqrt(sum(weights.*(thisTraj(itime,:)-meanFinalSSPoint).^2));
end

% get index of entry into final SS
belowThreshold = distToFinalSS<=maxFinalSSDist;
valididx = findContinousLogical(belowThreshold, nWins);
entryIndexFinalSS = valididx(1);

directDistance = sqrt(sum(weights.*(thisTraj(1,:)-thisTraj(entryIndexFinalSS,:)).^2));

%% get cumulative distance travelled from zero to entry point

if directDistance>0

distTravelledVector = nan(entryIndexFinalSS-1,1);

for itime = 2:entryIndexFinalSS
    prevPt = thisTraj(itime-1,:);
    thisPt = thisTraj(itime,:);
    distTravelledVector(itime-1) = sqrt(sum(weights.*(prevPt-thisPt).^2));
end

cumDist = cumsum(distTravelledVector);

cumulativeDistanceTravelled = cumDist(end);
distanceRatio = cumulativeDistanceTravelled/directDistance;


%% get the distance between the trajectory and the initial SS

distFromInitialSS = nan(size(thisTraj,1),1);
for itime = 1:size(thisTraj,1)
    distFromInitialSS(itime) = sqrt(sum(weights.*(thisTraj(itime,:)-meanInitialSSPoint).^2));
end

% get index of entry into final SS
aboveThreshold = distToFinalSS>=maxInitialSSDist;
valididx = findContinousLogical(aboveThreshold, nWins);
if ~isempty(valididx)
exitIndexInitialSS = valididx(1);
else
exitIndexInitialSS=1;
end

%% get valid exit and entry indices + time points

newBinVector = timeBinVector(zeroIdx:end); % to account for thisTraj 

exitInitialSSLatency = newBinVector(exitIndexInitialSS);
enterFinalSSLatency = newBinVector(entryIndexFinalSS);

else % if direct distance == 0
    cumulativeDistanceTravelled = nan;
    distanceRatio = nan;
    distFromInitialSS = nan;
    distToFinalSS = nan;
    distTravelledVector = nan;
    exitInitialSSLatency = nan;
    enterFinalSSLatency = nan;
end
%% assign outputs to output struct

output.maxInitialSSDist = maxInitialSSDist;
output.maxFinalSSDist = maxFinalSSDist;
output.distFromInitialSS = distFromInitialSS;
output.distFromFinalSS = distToFinalSS;
output.directDistance = directDistance;
output.cumulativeDistanceTravelled = cumulativeDistanceTravelled;
output.distanceRatio = distanceRatio;
output.distTravelledVector = distTravelledVector;
output.exitInitialSSLatency = exitInitialSSLatency;
output.enterFinalSSLatency = enterFinalSSLatency;

end



