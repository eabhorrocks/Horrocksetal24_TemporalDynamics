function predcond = PID_DecodeFromSpikeCellArray(tempUnits, tuningAxes)

% temp units has the field spikeCell, cell array of spike counts which has
% the shape of the neural tuning being used

tuningShape = size(tempUnits(1).spikeCell);
numDims = sum(size(tempUnits(1).spikeCell)>1);
nReps = numel(tempUnits(1).spikeCell{1});
nCells = numel(tempUnits);
nanpred = NaN*ones(tuningShape);
nConds = numel(tempUnits(1).spikeCell);

for icond = 1:nConds
    predcond(icond).preds = nan*ones(nReps,numDims);
end

% LOOOCV
allTrialIndexes = nchoosek(1:nReps,nReps-1); % each row is trials to train on:

for irep = 1:nReps
    trainTrials = allTrialIndexes(irep,:);
    testTrials = find(~ismember(1:nReps,trainTrials));
    tempCells = [];
    
    %%%%% TRAIN DECODER -> LEARN TUNING CURVES %%%%%
    for iunit = 1:nCells
        t = [tempUnits(iunit).spikeCell{:}]; 
        t = t(trainTrials,:);
        tuning = mean(t);
        tuning = reshape(tuning, tuningShape);
        %tempCells(iunit).tuning = tuning; %interp2(opts.X,opts.Y,tuning,opts.Xq,opts.Yq);
        tempCells(iunit).tuning = tuning + eps; % for log(0)
        clear t tuning
    end
    
    PID = PoissonIndependentDecoder;
    PID.tuningAxes = tuningAxes;
    PID = PID.trainDecoder(tempCells);
       
    %%%%%  TEST DECODER -> GENERATE PREDICTION %%%%%
    for itrial = 1:numel(testTrials) % for leave one out, just 1.
        for icond = 1:nConds
            testSpikeCounts = nan*ones(1,nCells);
            for iunit = 1:nCells
                t = [tempUnits(iunit).spikeCell{:}];
                testSpikeCounts(iunit) = t(testTrials(itrial),icond);
            end

            predcond(icond).preds(testTrials(itrial),:) = PID.testDecoder(testSpikeCounts);
        end
    end
end

end