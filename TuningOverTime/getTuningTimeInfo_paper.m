function [tunedStart, tunedFinish, nBouts, tunedDuration, ...
    tuneFlagWeak, tuneFlagMed, tuneFlagStrong] = ...
    getTuningTimeInfo_paper(R2Vector, R2pVector)



%% get tuning start and finish times + logical tuning flag

binSize = 0.01;
windowLength = 0.2;
binVector = -0.2:binSize:1.8;
binVector = round(binVector,2);
binStarts = 1:181;
binMidPoints = binVector(binStarts)+windowLength/2;

%R2Vector = [intervals.R2];

% get tuning bouts (R2>=0.1) for at least minWinLength
r2_threshold = 0.1; %r2 tuning threshold
r2_pThreshold = 0.05;
minWinLength = 5; % minimum consecuitve windows for a tuning bout

% code to find tuning bouts
aboveThreshold = (R2Vector >= r2_threshold & R2pVector<=r2_pThreshold);  %where above threshold
aboveThreshold = [false, aboveThreshold, false];  %pad with 0's at ends
edges = diff(aboveThreshold);
rising = find(edges==1);     %rising/falling edges
falling = find(edges==-1);
spanWidth = falling - rising;  %width of span of 1's (above threshold)
wideEnough = spanWidth >= minWinLength;
startPos = rising(wideEnough);    %start of each span
endPos = falling(wideEnough)-1;   %end of each span
allInSpan = cell2mat(arrayfun(@(x,y) x:1:y, startPos, endPos, 'uni', false));

if numel(startPos>0)
    tunedStart = binMidPoints(startPos(1));
    tunedFinish = binMidPoints(endPos(end));
    nBouts = numel(startPos);
    tunedDuration = numel(allInSpan);
    tuneFlagWeak = R2Vector>=0.1 & R2pVector<=r2_pThreshold;
    tuneFlagMed = R2Vector>=0.3 & R2pVector<=r2_pThreshold;
    tuneFlagStrong = R2Vector>=0.5 & R2pVector<=r2_pThreshold;
    
else
    tunedStart = nan;
    tunedFinish = nan;
    nBouts = 0;
    tunedDuration = nan;
    tuneFlagWeak = R2Vector>=0.1 & R2pVector<=r2_pThreshold;
    tuneFlagMed = R2Vector>=0.3 & R2pVector<=r2_pThreshold;
    tuneFlagStrong = R2Vector>=0.5 & R2pVector<=r2_pThreshold;
end