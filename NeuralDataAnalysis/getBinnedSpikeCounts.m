function [anUnits, cond] = getBinnedSpikeCounts(trials, units, varsOfInterest, options)

% outputs are:
% - spikeCounts (trials,units) array for each trial condition (this is a
% useful object for decoding
% - high dim tuning object for each unit
% - mean 1d tuning for each object + pref var value

% example inputs:
%trials = trials(strcmp([trials.stimulus_name],'dot_motion'));
%varsOfInterest = {'Dir', 'Speed'};
%units (struct array of units with raw spike times)

if ~exist('options','var'),             options=struct;                 end
% options concerning which parts of the function to run
if ~isfield(options,'intervalStart'),   options.intervalStart=0;        end
if ~isfield(options,'binSpacing'),      options.binSpacing=0.1;         end
if ~isfield(options,'intervalEnd'),     options.intervalEnd=1;          end
if ~isfield(options,'uniqueVals'),      options.uniqueVals = [];        end


nVars = numel(varsOfInterest);
nUnits = numel(units);

for ivar = 1:nVars
    if isempty(options.uniqueVals)
    vars(ivar).uniqueVals = nanunique([trials.(varsOfInterest{ivar})]);
    else
        vars(ivar).uniqueVals = options.uniqueVals{ivar};
    end
    vars(ivar).nVals = numel(vars(ivar).uniqueVals);
    vars(ivar).trialVals = [trials.(varsOfInterest{ivar})];
end


allVarCombs = combvec(vars(:).uniqueVals); % all unique var combinaitions
trialVarsCat = vertcat(vars(:).trialVals); % trial params vertically catenated
nVarCombs = size(allVarCombs,2);
anUnits(nUnits).cond = []; % initialised analysed units struct


% initialise cond struct array
cond(nVarCombs).varVals = [];
cond(nVarCombs).trials = [];
cond(nVarCombs).spikeCounts =[];

nBins = numel(options.intervalStart:options.binSpacing:options.intervalEnd);

for ivarcomb = 1:nVarCombs
    cond(ivarcomb).varVals = allVarCombs(:,ivarcomb)';
    cond(ivarcomb).trials = trials(ismember(trialVarsCat', allVarCombs(:,ivarcomb)', 'rows'));
    
    % initialise spike counts array with size (nBins, nReps, nUnits)
    cond(ivarcomb).spikeCounts = NaN*ones(nBins-1, numel(cond(ivarcomb).trials), nUnits); %-1 for bins vs edges
    
    if numel(cond(ivarcomb).trials)>0
    for itrial = 1:numel(cond(ivarcomb).trials)
        binEdges = (cond(ivarcomb).trials(itrial).start_time+options.intervalStart):...
            options.binSpacing:(cond(ivarcomb).trials(itrial).start_time+options.intervalEnd);
        
        for iunit = 1:nUnits
            cond(ivarcomb).spikeCounts(:,itrial,iunit) =...
                histcounts(units(iunit).spiketimes, binEdges);
            anUnits(iunit).cond(ivarcomb).spikeCounts(:,:) = cond(ivarcomb).spikeCounts(:,:,iunit);
        end
        
    end
    else % if there are no trials
        cond(ivarcomb).spikeCounts = [];
        for iunit = 1:nUnits
        anUnits(iunit).cond(ivarcomb).spikeCounts = [];
        end
    end
    
end
% put binned spikes into cell array
if numel([vars.nVals])>1
    outputShape = [vars.nVals];
else
    outputShape = [1,vars.nVals];
end
for iunit=1:nUnits
    anUnits(iunit).allSpikes = reshape({anUnits(iunit).cond.spikeCounts},outputShape);
end


