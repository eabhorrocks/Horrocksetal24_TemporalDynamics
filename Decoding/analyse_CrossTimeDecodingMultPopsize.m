%% analyse cross-time decoding
sessionTags = {'M22027', '20220517';...
    'M22029', '20220607';...
    'M22031', '20220524';
    'M22032', '20220621';
    'M22033', '20220706'};

dataDir = 'C:\Users\edward.horrocks\Documents\Code\V1Dynamics\Data\basic_111022';
dataDir = 'E:\V1Data\Data\v1_fromC24';

for isession = 1:size(sessionTags,1)
    fname = [sessionTags{isession,1},'_', sessionTags{isession,2},'_crossTimeDecoding_popSize_June24.mat']; %PSTH_noEye %psth3rds

    load(fullfile(dataDir,fname))
    [units.session] = deal(isession);

    s(isession) = session;



end

%%

popSizeVector = [10 20 40 80];
nWorkers=18;


for ipop = 1:4

    allStat=[];
    allRun=[];

    for isession = 1:5
        s(isession).popSize(ipop).allStat = [];
        s(isession).popSize(ipop).allRun = [];

        % get decoding performance of all reps
        for irep =1:numel(s(isession).popSize(ipop).rep)
            thisStatPerf = mean(cat(3,s(isession).popSize(ipop).rep(irep).stat.perm.perf),3);
            thisRunPerf = mean(cat(3,s(isession).popSize(ipop).rep(irep).run.perm.perf),3);

            s(isession).popSize(ipop).allStat = cat(3,s(isession).popSize(ipop).allStat,thisStatPerf);
            s(isession).popSize(ipop).allRun = cat(3,s(isession).popSize(ipop).allRun,thisRunPerf);
        end
    end

    timeIdx2Use =3:12;
    trainIdx=4;

    % calculate fractional performance of different train test windows
    for isession = 1:5
        s(isession).popSize(ipop).statDiff=[];
        s(isession).popSize(ipop).runDiff=[];

        if ~isempty(s(isession).popSize(ipop).rep)

            for irep =1:numel(s(isession).popSize(ipop).rep)

                thisStatPerf = s(isession).popSize(ipop).allStat(:,:,irep);
                thisRunPerf = s(isession).popSize(ipop).allRun(:,:,irep);

                diagPerf_stat = diag(thisStatPerf);
                diagPerf_stat = diagPerf_stat-1/6;
                diagPerf_stat(diagPerf_stat<0)=0;
                diagPerf_stat = diagPerf_stat(timeIdx2Use);

                diagPerf_run = diag(thisRunPerf);
                diagPerf_run = diagPerf_run-1/6;
                diagPerf_run(diagPerf_run<0)=0;
                diagPerf_run = diagPerf_run(timeIdx2Use);


                for itime = 1:timeIdx2Use(end)
                    tempPerf = thisStatPerf(itime,:)-1/6; tempPerf(tempPerf<0)=0;
                    tempPerf = tempPerf(timeIdx2Use);
                    thisDiffStatPerf(itime) = sum(tempPerf)/sum(diagPerf_stat);

                    tempPerf = thisRunPerf(itime,:)-1/6; tempPerf(tempPerf<0)=0;
                    tempPerf = tempPerf(timeIdx2Use);
                    thisDiffRunPerf(itime) = sum(tempPerf)/sum(diagPerf_run);

                end

                s(isession).popSize(ipop).statDiff=cat(1,s(isession).popSize(ipop).statDiff,thisDiffStatPerf);
                s(isession).popSize(ipop).runDiff=cat(1,s(isession).popSize(ipop).runDiff,thisDiffRunPerf);

            end

        end


    end
end

%% plot normalised performance for all reps

figure
for ipop = 1:4

    allStatDiff=[];
    allRunDiff=[];

    for isession = 1:size(sessionTags,1)

        allStatDiff = cat(1,allStatDiff, s(isession).popSize(ipop).statDiff);
        allRunDiff = cat(1,allRunDiff, s(isession).popSize(ipop).runDiff);

        allStatDiff(allStatDiff>1) = 1;
        allRunDiff(allRunDiff>1) = 1;
    end

    trainIdx = 4;
    vals = {allStatDiff(:,trainIdx), allRunDiff(:,trainIdx)};
     hold on
    distributionPlot(vals,'histOpt',1,'globalNorm',3,'showMM',6,'divFactor',3,...
        'xValues',[ipop-0.2, ipop+0.2],'color',{'k','r'})
    plot([0 4],[0 0],'k:')
    ylabel('Relative Performanace')
    xlabel('pop size')
    n(ipop)=numel(vals{1})
end
ax =gca; ax.XTick = 1:4;
defaultAxesProperties(gca, true)
% figure(999)
% subplot(1,4,ipop)
% bar([1 2], [mean(allStatDiff(:,trainIdx)), mean(allRun(:,trainIdx))])
% for isesh = 1:5
%     plot([1 2], [allStat(isesh,trainIdx), allRun(isesh,trainIdx)])
%     pause
% end
% 




%%



% figure(999), hold on
%
% vals = {allStatDiff(:,trainIdx), allRunDiff(:,trainIdx)};
% temp{ipop} = [vals{1},vals{2}];
% hold on
% distributionPlot(vals,'histOpt',1,'globalNorm',3,'showMM',6,'divFactor',2,...
%     'xValues',[ipop-0.2, ipop+0.2],'color',{'k','r'})
% plot([0 4],[0 0],'k:')
% ylabel('Relative Performanace')
% xlabel('pop size')

% statVals = allStatDiff(:,trainIdx);
% runVals = allRunDiff(:,trainIdx);
% 
% figure(999)
% subplot(1,4,ipop)
% bar([1 2], [mean(allStatDiff(:,trainIdx)), mean(allRun(:,trainIdx))])
% for isesh = 1:5
%     plot([1 2], [allStat(isesh,trainIdx), allRun(isesh,trainIdx)])
%     pause
% end
% 




%%


