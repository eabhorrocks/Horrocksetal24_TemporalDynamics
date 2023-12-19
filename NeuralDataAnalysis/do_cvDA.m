function [meanPerf, semPerf, pCorrect, cond, meanError, semError, confMatrix] = do_cvDA(dataStruct, options)

%% input args
% datastruct is of format: struct(icondition).data(time,unit,trial)

if ~exist('options','var'),                     options=struct;                         end

if ~isfield(options,'kfold'),                   options.kfold=inf;                      end % LOOCV by default
if ~isfield(options,'DiscrimType'),             options.DiscrimType='linear';           end % linear by default
% hyperparameter optimization options
if ~isfield(options,'OptimizeHyperparameters'), options.OptimizeHyperparameters='none'; end % no opt by defualt
if ~isfield(options,'Gamma'),                   options.Gamma=[];                       end

%% basic info
nConds = numel(dataStruct); % number of classes
nUnits = size(dataStruct(1).data,2);
nTrials = size(dataStruct(1).data,3); % assumes equal trial count for each class (error check here)

%% remove units with 0 variance

% for iunit = 1:nUnits
%     tempVec=[];
%     for istim = 1:nConds
%         tempVec=cat(1,tempVec,squeeze(mean(dataStruct(istim).data(:,iunit,:),1)));
%     end
%     unitVar(iunit)=var(tempVec);
% end
% 
% idx2remove = find(unitVar<eps);
% 
% for istim = 1:nConds
%     dataStruct(istim).data(:,idx2remove,:)=[];
% end
% 
% nUnits = nUnits - numel(idx2remove);

%% Main


if options.kfold>=nTrials % check for leave-one-out cv
    options.kfold=nTrials; 
end

nReps = options.kfold; 
cvobj = cvpartition(nTrials, 'kfold', options.kfold); % generate train and test sets

cond(nConds).preds=[]; % struct array to store predictions in
pCorrect = nan(options.kfold,1); % array to store pCorrect for each k-fold
decodingError = nan(options.kfold,1); % array to store decoding errors for each k-fold

for irep = 1:nReps % loop over k-folds

    % get training and testing data for this k-fold
    trainTrials = find(training(cvobj,irep)); % logical index
    testTrials = find(test(cvobj,irep));
    nTrain = numel(trainTrials);
    nTest = numel(testTrials);

    trainingData = nan(nConds*nTrain, nUnits);
    trainingLabels = nan(nConds*nTrain,1);
    testingData = nan(nConds*nTest, nUnits);
    testingLabels = nan(nConds*nTest, 1);

    ii = 0;
    for icond = 1:nConds
        for itrial = 1:numel(trainTrials)
            ii = ii+1;
            trainingData(ii,:) = mean(dataStruct(icond).data(:,1:end,trainTrials(itrial),1),1); % average over time bins
            trainingLabels(ii,1) = icond;
        end
    end
    
    ii = 0;
    for icond = 1:nConds
        for itrial = 1:numel(testTrials)
            ii = ii+1;
            testingData(ii,:) =  mean(dataStruct(icond).data(:,1:end,testTrials(itrial),1),1); % average over time bins
            testingLabels(ii,1) = icond;
        end
    end    

    % remove 0 variance predictors
    idx2remove = find(var(trainingData,1)<eps);
    trainingData(:,idx2remove)=[];
    testingData(:,idx2remove)=[];
    
    % fit the discriminant analysis model on the training data
    Mdl = fitcdiscr(trainingData,trainingLabels,'DiscrimType', ...
        options.DiscrimType,...
        'OptimizeHyperparameters',options.OptimizeHyperparameters,...
        'Gamma',options.Gamma,...
        'HyperparameterOptimizationOptions',...
        struct('ShowPlots',false,...
        'AcquisitionFunctionName','expected-improvement-plus',...
        'UseParallel',true,'Verbose',0));

    %Mdl_out(irep).mdl = Mdl;

    % generate class predictions for the testing data
    preds = predict(Mdl,testingData);

    % add these predictions to existing predictions from previous k-folds
%     for icond = 1:nConds
%         cond(icond).preds = cat(1,cond(icond).preds, preds(testingLabels==icond));
%     end

    for icond = 1:nConds
        cond(icond).preds(testTrials) = preds(testingLabels==icond);
    end

    % get pCorrect and decoding error for this kfold
    pCorrect(irep,1) = prop(preds==testingLabels);
    decodingError(irep,1) = mean(abs(preds-testingLabels));
    
end % finished all k-folds

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
confMatrix = nan(nConds);

for icond1 = 1:nConds
    for icond2 = 1:nConds
        confMatrix(icond1,icond2) = prop(cond(icond1).preds == icond2);
    end
end


end

