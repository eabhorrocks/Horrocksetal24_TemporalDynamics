function [modelParams, R2] =...
    fitDescFunAndValidate_spikeCounts(spikeCountCellArray,xAxisVals, k, randFlag,...
    fitMeans, validMeans, customFun, x0, lb, ub, nSolvers)

%spikeCountCellArray = units(4).dmSC_trial_speed;
%xAxisVals = 1:7;
if numel(xAxisVals)~=numel(spikeCountCellArray)
    error('xAxisVals must have same number of elements as spikeCountCellArray')
end
%% cross-validated model fitting

%randFlag = 1;
%k = 5; % i.e. train on 80% of data, test on remaining 20%
nConds = numel(spikeCountCellArray);


%% get train and test data
ytestData = cell(size(spikeCountCellArray));
ytrainData = ytestData;

if randFlag % randomly pick train and test trials
    for icond = 1:nConds
        nReps = size(spikeCountCellArray{icond},1);
        nTest = round(nReps/k);
        shuf = randperm(nReps);
        ytestData{icond} = spikeCountCellArray{icond}(shuf(1:nTest));
        ytrainData{icond} = spikeCountCellArray{icond}(shuf(nTest+1:end));
    end
else % sequential segmenting of train/test trials
    for icond = 1:nConds
        nReps = size(spikeCountCellArray{icond},1);
        nTest = round(nReps/k);
        ytestData{icond} = spikeCountCellArray{icond}(1:nTest);
        ytrainData{icond} = spikeCountCellArray{icond}(nTest+1:end);
    end
end

if fitMeans
    ytestData = cellfun(@(x) num2cell(mean(x)), ytestData); % cell fo mean spikes
    yTrainVals = cellfun(@mean, ytrainData); % vector of spike counts
    xTrainVals = xAxisVals; % vector of stim vals
else % fit all the data
    xTrainVals = [];
    for icond = 1:nConds
        xTrainVals = [xTrainVals; repelem(xAxisVals(icond),numel(ytrainData{icond}),1)];
    end
    yTrainVals = vertcat(ytrainData{:});
end

%% fit function on training data...
ms = MultiStart();
problem = createOptimProblem('lsqcurvefit','x0',x0,'objective',customFun,...
    'lb', lb, 'ub', ub, 'xdata',xTrainVals,'ydata',yTrainVals);
% problem = createOptimProblem('fmincon','x0',params0,'objective',fun,...
% 'lb',lb,'ub',ub,'xdata',xvals,'ydata',yvals);
[modelParams,errormulti] = run(ms,problem,nSolvers);


%% evaluate function on test data

yfit = customFun(modelParams,xAxisVals);
y_train_const = mean(cellfun(@mean, ytrainData));

if validMeans
    ytestData = num2cell(cellfun(@mean, ytestData));
end

SSres_cond = nan([1, nConds]);
for icond = 1:nConds
    SSres_cond(icond) = sum((yfit(icond)-ytestData{icond}).^2);
end
SSres = sum(SSres_cond);
SStot = sum((y_train_const-vertcat(ytestData{:})).^2);

ratio = SSres/SStot;
if ratio > 1
    R2 = -1 + (1/ratio);
else
    R2 = 1-ratio;
end

% pause
% hold off
% yfit2 = customFun(modelParams,xAxisVals(1):1:xAxisVals(end));
% plot(xAxisVals, cellfun(@mean, testData), 'ko'), hold on
% plot(xAxisVals(1):1:xAxisVals(end),yfit2,'r-')
% title(['R^2 = ' num2str(R2,2), ', Params: ' num2str(modelParams,3)])

end

