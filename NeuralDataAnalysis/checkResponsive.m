function [pval, stats, distRatio] = checkResponsive(stimTrials, blankResp, distType)
%% info
% function determine whether PSTH is responsive
% logic: if a unit is responsive to a stimulus individual trial responses
% should be more similar to the average stimulus response (excluding that
% trial)  than to the average blank response (response to blank stimulus or 
% constant spontaneous rate).

% input args: 
% stimTrials array of size nTime x nTtrials
% blankResp: either a vector of length nTime or a single value
% distType: the distnace metric to use, e..g 'euclidean', can use a custom
% function

% output args:
% pval from signrank test on distance measures
% distRatio of distance meaasures for trial/stim to trial/spon

% example usage
% stimTrials = smoothdata(unitex.trials_stat{ispeed}.*100,1,'gaussian', 17.5);
% blankResp = sponSmoothed;
% blankResp = mean(smoothdata(unitex.blankTrials_stat{1}.*100,1,'gaussian', 17.5));
% distType = 'euclidean';
% [pval, stats, distRatio] = checkResponsive(stimTrials, blankResp, distType)

%% main function
% if blankResp is a single value firing rate, make it a flat response
if numel(blankResp)==1
    blankResp = blankResp*ones(size(stimTrials,1),1);
end

if numel(blankResp)~=size(stimTrials,1)
    error('blankResp must be either a single vale firing rate or a vector with nTime bins equal to stimTrials')
end

% loop through trials and calculate distance to meanStim and blankResp
for itrial = 1:size(stimTrials,2)
    trainTrials = stimTrials;
    trainTrials(:,itrial)=[];
    meanTrain = mean(trainTrials,2);
    testTrial = stimTrials(:,itrial);    
    
    stimDist(itrial) = pdist([meanTrain'; testTrial'],distType); % distance between training stim and testing stim
    sponDist(itrial) = pdist([blankResp'; testTrial'], distType); % distance between blank resp and test stim
end
 
% alt hypothesis is that stimDist should be smaller than sponDist
[pval, ~, stats] = signrank(stimDist, sponDist, 'tail', 'left'); 
distRatio = median(stimDist)/median(sponDist);


%% plots
% figure, hold on
% shadedErrorBar(1:200, mean(stimTrials,2), sem(stimTrials,2))
% plot(1:200, blankResp, 'r')
% figure, hold on
% mVal = max(vertcat(stimDist(:), sponDist(:)));
% plot([0 mVal], [0 mVal])
% plot(stimDist, sponDist, 'ko')
% xlabel('Stim distance'), ylabel('Spon Distance')
% title(pval)
% pause=1;
% close all

end
