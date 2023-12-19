function [meanPerf, semPerf, pCorrect, cond, meanError, semError, confMatrix, Mdl_out] = doLDA(dataStruct, options)

%% input args
% datastruct is of format: struct(icondition).data(time,unit,trial)

%%

nConds = numel(dataStruct);

nTrials = size(dataStruct(1).data,3);
if options.kfold>=nTrials
    options.kfold=nTrials; %LOOCV
    options.nShuffle = 1; %no point repeating LOOCV
end

nReps = options.kfold;

cond(nConds).preds=[];

cvobj = cvpartition(nTrials, 'kfold', options.kfold); % generate train and test sets

for irep = 1:nReps
    trainTrials = find(training(cvobj,irep)); % logical index
    testTrials = find(test(cvobj,irep));
    
    nTrain = sum(trainTrials);
    nTest = sum(testTrials);
    
    ii = 0;
    for icond = 1:nConds
        for itrial = 1:numel(trainTrials)
            ii = ii+1;
            trainingData(ii,:) = sum(dataStruct(icond).data(:,1:end,trainTrials(itrial),1),1);
            trainingLabels(ii,1) = icond;
        end
    end
    
    ii = 0;
    for icond = 1:nConds
        for itrial = 1:numel(testTrials)
            ii = ii+1;
            testingData(ii,:) =  sum(dataStruct(icond).data(:,1:end,testTrials(itrial),1),1);
            testingLabels(ii,1) = icond;
        end
    end    
    
    %     Mdl = fitcecoc(trainingData,trainingLabels);
    %     try
    %         % fast version
    %         preds = classify(testingData,trainingData,trainingLabels,'linear');
    %     catch
    % slow version
    Mdl = fitcdiscr(trainingData,trainingLabels, 'DiscrimType', 'diaglinear');
    Mdl_out(irep).mdl = Mdl;
    preds = predict(Mdl,testingData);
    %     end
    %
    
    % get the predictions for each speed shown and add to the existing
    % predicitons
    for icond = 1:nConds
        cond(icond).preds = cat(1,cond(icond).preds, preds(testingLabels==icond));
    end
    % get pCorrect and decoding error for this kfold
    pCorrect(irep,1) = prop(preds==testingLabels);
    decodingError(irep,1) = mean(abs(preds-testingLabels));
    
end

% average perf over kfold cross-validations
meanPerf = mean(pCorrect,1); 
semPerf = sem(pCorrect,1);

meanError = mean(decodingError,1);
semError = sem(decodingError,1);

% get performance for each condition
for icond = 1:nConds
    cond(icond).perf = prop(cond(icond).preds==icond);
    cond(icond).meanError = mean(abs(cond(icond).preds-icond));
    cond(icond).semError = sem(abs(cond(icond).preds-icond));
end

% confusion matrix
confMatrix = nan([nConds]);

for icond1 = 1:nConds
    for icond2 = 1:nConds
        confMatrix(icond1,icond2) = prop(cond(icond1).preds == icond2);
    end
end


end

