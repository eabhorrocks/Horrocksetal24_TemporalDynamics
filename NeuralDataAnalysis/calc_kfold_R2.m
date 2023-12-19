function [R2, pval] = calc_kfold_R2(spikeCountCellArray, kfold, nPerms, randFlag, validMeans, nShuffle)

% function calculates cross-validated R2 with option of multiple perms

% add check for columns vectors in cell array!

nConds = numel(spikeCountCellArray);
nReps = cellfun(@numel,spikeCountCellArray);
nTestTrials = arrayfun(@(x) floor(x/kfold), nReps); % test on 1/k
nTrainTrials = arrayfun(@(x) floor((x/kfold))*(kfold-1), nReps); % train on (k-1)/k
allTrials = arrayfun(@(x) 1:x, nReps, 'UniformOutput', false);

meanR2 = nan*ones(1,nPerms); % initialise nPerm R2 values

for iperm = 1:nPerms % repeat nPerm times
    
    R2 = nan*ones(1,kfold); % initalise kfold R2 value for this perm
    
    if randFlag % assign test trials randomly
        testTrials = cellfun(@(x) randperm(numel(x)), spikeCountCellArray', 'UniformOutput', false);
    else % or just sequentially
        testTrials = cellfun(@(x) 1:numel(x), spikeCountCellArray', 'UniformOutput', false);
    end
    
    for ik = 1:kfold % do kfold cross-validation
        for icond = 1:nConds
            testTrialsRep = testTrials{icond}((ik*nTestTrials(icond)-(nTestTrials(icond)-1)):(ik*nTestTrials(icond)));
            % component of test curve
            testData{icond} = spikeCountCellArray{icond}(testTrialsRep);
            
            % get component of mean tuning curvbe
            trainTrialsRep = allTrials{icond}(~ismember(allTrials{icond},testTrialsRep));
            trainedModel(icond) = mean(spikeCountCellArray{icond}(trainTrialsRep));
            
            % calc residual
            %resid(icond) = trainData(icond)-testData(icond);
        end
        
        if validMeans
            testData = num2cell(cellfun(@mean, testData)); % keep as a cell to stay compatible...
        end
        
        % null model
        y_train_const = mean(trainedModel);
        
        % get sum of squared residuals for trained tuning model and null mean model
        SSres_cond = nan([1, nConds]);
        for icond = 1:nConds
            SSres_cond(icond) = sum((trainedModel(icond)-testData{icond}).^2);
        end
        SSres = sum(SSres_cond);
        SStot = sum((y_train_const-vertcat(testData{:})).^2);
        
        ratio = SSres/SStot;
        if ratio > 1
            R2(ik) = -1 + (1/ratio); % adjust for -ve R2 values...
        else
            R2(ik) = 1-ratio;
        end
        
        % * optional plotting code goes here
        R2(ik);
        
%                 subplot(121)
%                 hold off
%                 plot([testData{:}], 'k--')
%                 hold on
%                 plot([1 7], [y_train_const y_train_const], 'r-') % null model
%                 plot(trainedModel, 'b-') % training curve model
%                 legend({'test', 'null', 'model'})
%                 box off
%                 defaultAxesProperties(gca, true)
%                 title('Data')
%         
%         
%                 subplot(122)
%                 hold off,
%                 plot(SSres_cond, 'b-o')
%                 hold on
%                 plot(((y_train_const-vertcat(testData{:})).^2), 'r-o')
%                 legend({'SSmodel', 'SSnull'})
%                 box off
%                 defaultAxesProperties(gca, true)
%                 title('SS residuals')
        
    end
    % get mean R2 for k fold cross-validations for this permutation
    meanR2(iperm) = mean(R2);
    
end

% take mean all permutations for final result
R2 = mean(meanR2);



%% if calculating p value, shuffle spike counts and compute R2
if nShuffle > 0    
    shuffleR2vals = nan(1,nShuffle);
    
    for ishuf = 1:nShuffle
        
        % shuffle spike count cell array
        allSpikes = cell2mat(spikeCountCellArray);
        allSpikesShuffled = allSpikes(randperm(numel(allSpikes)));
        
        for icond = 1:numel(spikeCountCellArray)
            shuffleSpikeCell{icond} = allSpikesShuffled(1:numel(spikeCountCellArray{icond}));
            allSpikesShuffled(1:numel(spikeCountCellArray{icond})) = [];
        end
        
        R2shuf = nan*ones(1,kfold); % initalise kfold R2 value for this perm

            % get test trials (already shuffled)
            testTrials = cellfun(@(x) 1:numel(x), shuffleSpikeCell', 'UniformOutput', false);
        
        for ik = 1:kfold % do kfold cross-validation
            for icond = 1:nConds
                testTrialsRep = testTrials{icond}((ik*nTestTrials(icond)-(nTestTrials(icond)-1)):(ik*nTestTrials(icond)));
                % component of test curve
                testData{icond} = shuffleSpikeCell{icond}(testTrialsRep);
                
                % get component of mean tuning curvbe
                trainTrialsRep = allTrials{icond}(~ismember(allTrials{icond},testTrialsRep));
                trainedModel(icond) = mean(shuffleSpikeCell{icond}(trainTrialsRep));
                
                % calc residual
                %resid(icond) = trainData(icond)-testData(icond);
            end
            
            if validMeans
                testData = num2cell(cellfun(@mean, testData)); % keep as a cell to stay compatible...
            end
            
            % null model
            y_train_const = mean(trainedModel);
            
            % get sum of squared residuals for trained tuning model and null mean model
            SSres_cond = nan([1, nConds]);
            for icond = 1:nConds
                SSres_cond(icond) = sum((trainedModel(icond)-testData{icond}).^2);
            end
            SSres = sum(SSres_cond);
            SStot = sum((y_train_const-vertcat(testData{:})).^2);
            
            ratio = SSres/SStot;
            if ratio > 1
                R2shuf(ik) = -1 + (1/ratio); % adjust for -ve R2 values...
            else
                R2shuf(ik) = 1-ratio;
            end
            
            
        end
        
        shuffleR2vals(ishuf) = mean(R2shuf);
    end
    
    % now calculate p-value by comparing R2 value to shuffled values
     pval = 1 - prop(R2>shuffleR2vals);
            
end