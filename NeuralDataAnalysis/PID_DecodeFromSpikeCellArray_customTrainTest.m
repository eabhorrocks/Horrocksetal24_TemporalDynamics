function predcond = PID_DecodeFromSpikeCellArray_customTrainTest(tempUnits, tuningAxes)

% temp units has the fields trainSpikeCell (learn tuning) and testSpikeCell 
% (spike counts to gen predictions from). Shapes are of the tuning shape.

tuningShape = size(tempUnits(1).trainSpikeCell);
numDims = sum(size(tempUnits(1).trainSpikeCell)>1);
nReps = numel(tempUnits(1).trainSpikeCell{1});
nCells = numel(tempUnits);
nanpred = NaN*ones(tuningShape);
nConds = numel(tempUnits(1).trainSpikeCell);

for icond = 1:nConds
    predcond(icond).preds = [];
end

%%%%% TRAIN DECODER -> LEARN TUNING CURVES %%%%%
for iunit = 1:nCells
    % tuning as mean of spike counts
    tempCells(iunit).tuning = cellfun(@mean, tempUnits(iunit).trainSpikeCell);
    %interp2(opts.X,opts.Y,tuning,opts.Xq,opts.Yq);
    tempCells(iunit).tuning = tempCells(iunit).tuning + eps; % for log(0)
end

PID = PoissonIndependentDecoder;
PID.tuningAxes = tuningAxes;
PID = PID.trainDecoder(tempCells);

%%%%%  TEST DECODER -> GENERATE PREDICTION %%%%%
for icond = 1:nConds
    for itrial = 1:numel(tempUnits(1).testSpikeCell{icond}) % for leave one out, just 1.
        testSpikeCounts = nan*ones(1,nCells);
        for iunit = 1:nCells
            testSpikeCounts(iunit) = tempUnits(iunit).testSpikeCell{icond}(itrial);
        end
        
        predcond(icond).preds(itrial) = PID.testDecoder(testSpikeCounts);
    end
end
end
