function intervals =...
    calcR2overtime_paper(spikeCellArray,binLengths,binStarts,opts)

%% get metrics over different time windows
% function simply loops through required intervals (defined in bin space)
% and

%% inputs:
% - spikeCellArray: each cell is a condition, with nTrial cols and nBin rows
% - binLengths: length of intervals in bins
% - binStarts: starting bin for each interval

% - opts (struct with options for cross-val R2 calc)
% { nPerms, randFlag, validMeans }


%% outputs:
% a struct array (nInts length), "intervals", with the following fields:
% - % tuning (mean), tuning (sem),
% - mean spike count over interval (mean over conditions)
% - cross-validated R2 value


%% configure desired intervals and set options
if (numel(binStarts)~=1 && numel(binLengths)~=1)
    if (numel(binStarts)~=numel(binLengths))
        error('binStarts and binLengths must be vectors of equal length or one must be a vector of length 1')
    end
end

if numel(binStarts)==1, binStarts = repelem(binStarts,numel(binLengths),1);      end
if numel(binLengths)==1, binLengths = repelem(binLengths(:),numel(binStarts),1); end

binStarts = binStarts(:);
binLengths = binLengths(:);

% cross-val R2 options
nPerms = opts.nPerms;
randFlag = opts.randFlag;
validMeans = opts.validMeans;
kval = opts.kval;
plotFlag = false;
nShuffle = opts.nShuffle;

% initialise variables
nInts = numel(binStarts);
nConds = numel(spikeCellArray);


intervals(nInts).tuning = [];
intervals(nInts).tuningErr = [];
intervals(nInts).meanCounts = [];
intervals(nInts).R2 = [];
intervals(nInts).R2_p = [];
intervals(nInts).dynamicRange = [];
intervals(nInts).meanFanoFactor = [];


% loop through interval windows and compute MI etc.
for iint = 1:numel(binStarts)
    
    %% get spike cell array for this interval
    
    % bins for this window
    sumbins = binStarts(iint):(binStarts(iint)+binLengths(iint)-1);
    % compute the spike cell array for this interval
    temp_sc_array = cell(1,nConds);
    for icond = 1:nConds
        temp_sc_array{icond} =...
            sum(spikeCellArray{icond}(sumbins,:),1)';
    end
    
    temp_sc_array = reshape(temp_sc_array, size(spikeCellArray));
    
    %% tuning curves
    intervals(iint).tuning = cellfun(@mean, temp_sc_array).*5; % *5 to convert to Hz (200ms window)
    intervals(iint).tuningErr = cellfun(@sem, temp_sc_array).*5;
    intervals(iint).meanCounts =  mean(intervals(iint).tuning(:));
    intervals(iint).dynamicRange = range(intervals(iint).tuning);
    intervals(iint).meanFanoFactor = nanmean(cellfun(@(x) var(x)/mean(x), temp_sc_array));
    
    %% cross-val R2
    [intervals(iint).R2, intervals(iint).R2_p] =...
        calc_kfold_R2(temp_sc_array, kval, nPerms, randFlag, validMeans, nShuffle);
    
end

